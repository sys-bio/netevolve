// M3: integrators (Dormand-Prince RK45, Rosenbrock RODAS3), guards and policy.
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "integrator.hpp"
#include "model.hpp"
#include "network.hpp"

using namespace ne;

namespace {

SimResult run(const Network& net, double t_end, int n_points, const SolverOptions& o) {
    const Model m = make_model(net);
    return simulate(m, net.init, t_end, n_points, o);
}

SolverOptions with(Method method, double rtol = 1e-6, double atol = 1e-9) {
    SolverOptions o;
    o.method = method;
    o.rtol = rtol;
    o.atol = atol;
    return o;
}

std::string name(Method m) {
    switch (m) {
        case Method::Auto: return "Auto";
        case Method::RK45: return "RK45";
        case Method::Rosenbrock: return "Rosenbrock";
        case Method::CVODE: return "CVODE";
    }
    return "?";
}

// Mean period of species `sp` from upward crossings of `level` after t_from,
// with linear interpolation between samples.
double mean_period(const SimResult& r, int sp, double level, double t_from) {
    std::vector<double> crossings;
    for (int i = 1; i < r.n_samples(); ++i) {
        if (r.t[i - 1] < t_from) continue;
        const double a = r.at(i - 1, sp) - level, b = r.at(i, sp) - level;
        if (a < 0 && b >= 0) crossings.push_back(r.t[i - 1] + (r.t[i] - r.t[i - 1]) * (-a) / (b - a));
    }
    REQUIRE(crossings.size() >= 3);
    return (crossings.back() - crossings.front()) / static_cast<double>(crossings.size() - 1);
}

const Network kDecay = parse_antimony("const X0;\nJ0: S0 -> X0; k0*S0;\nk0 = 0.5;\nS0 = 1;\nX0 = 0;\n");

// Robertson (1966), written as mass action:
//   y1 -> y2 (0.04);  y2 + y3 -> y1 + y3 (1e4);  2 y2 -> y2 + y3 (3e7)
const Network kRobertson = parse_antimony(
    "J0: S0 -> S1; k0*S0;\n"
    "J1: S1 + S2 -> S0 + S2; k1*S1*S2;\n"
    "J2: S1 + S1 -> S1 + S2; k2*S1*S1;\n"
    "k0 = 0.04; k1 = 1e4; k2 = 3e7;\n"
    "S0 = 1; S1 = 0; S2 = 0;\n");

}  // namespace

TEST_CASE("LU factorisation solves small dense systems") {
    Rng rng(11);
    for (int trial = 0; trial < 200; ++trial) {
        const int n = rng.uniform_int(1, 20);
        std::vector<double> A(static_cast<std::size_t>(n) * n), x(n), b(n, 0.0);
        for (double& a : A) a = rng.uniform(-1.0, 1.0);
        for (int i = 0; i < n; ++i) A[static_cast<std::size_t>(i) * n + i] += (i % 2 ? -1.0 : 1.0) * 0.1;
        for (double& v : x) v = rng.uniform(-5.0, 5.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) b[i] += A[static_cast<std::size_t>(i) * n + j] * x[j];
        std::vector<double> LU = A;
        std::vector<int> piv(n);
        REQUIRE(lu_factor(LU, n, piv));
        lu_solve(LU, n, piv, b);
        // Residual check: A * x_solved should reproduce A * x.
        double worst = 0.0;
        for (int i = 0; i < n; ++i) {
            double r = 0.0, rhs_i = 0.0;
            for (int j = 0; j < n; ++j) {
                r += A[static_cast<std::size_t>(i) * n + j] * b[j];
                rhs_i += A[static_cast<std::size_t>(i) * n + j] * x[j];
            }
            worst = std::max(worst, std::abs(r - rhs_i));
        }
        REQUIRE(worst < 1e-9);
    }
    std::vector<double> singular = {1, 2, 2, 4};
    std::vector<int> piv(2);
    CHECK_FALSE(lu_factor(singular, 2, piv));
}

TEST_CASE("dy/dt = -k y matches the analytic solution at t = 10 (both solvers)") {
    for (Method method : {Method::RK45, Method::Rosenbrock}) {
        CAPTURE(name(method));
        const SimResult r = run(kDecay, 10.0, 101, with(method, 1e-9, 1e-14));
        REQUIRE(r.ok());
        REQUIRE(r.n_samples() == 101);
        CHECK(r.t.back() == 10.0);
        const double exact = std::exp(-0.5 * 10.0);
        CHECK(std::abs(r.at(100, 0) - exact) / exact < 1e-6);
        // Dense output at the interior sample times is accurate too.
        double worst = 0.0;
        for (int i = 0; i < r.n_samples(); ++i) {
            const double e = std::exp(-0.5 * r.t[i]);
            worst = std::max(worst, std::abs(r.at(i, 0) - e) / e);
        }
        CHECK(worst < 1e-6);
    }
}

TEST_CASE("Lotka-Volterra: conserved quantity drifts < 1e-5 relative over t = 50") {
    // dx/dt = a x - b x y, dy/dt = b x y - c y with a = k0*X0 = 1, b = 0.5, c = 0.75.
    // Conserved: V = b x - c ln x + b y - a ln y.
    const Network net = parse_antimony(
        "const X0, X1;\n"
        "J0: X0 + S0 -> S0 + S0; k0*X0*S0;\n"
        "J1: S0 + S1 -> S1 + S1; k1*S0*S1;\n"
        "J2: S1 -> X1; k2*S1;\n"
        "k0 = 1; k1 = 0.5; k2 = 0.75;\n"
        "S0 = 2; S1 = 1;\n"
        "X0 = 1; X1 = 0;\n");
    const double a = 1.0, b = 0.5, c = 0.75;
    auto V = [&](double x, double y) { return b * x - c * std::log(x) + b * y - a * std::log(y); };

    for (Method method : {Method::RK45, Method::Auto}) {
        CAPTURE(name(method));
        const SimResult r = run(net, 50.0, 5001, with(method, 1e-10, 1e-12));
        REQUIRE(r.ok());
        const double v0 = V(r.at(0, 0), r.at(0, 1));
        double worst = 0.0;
        for (int i = 0; i < r.n_samples(); ++i) worst = std::max(worst, std::abs(V(r.at(i, 0), r.at(i, 1)) - v0) / std::abs(v0));
        MESSAGE(name(method) << ": max relative drift of V = " << worst << ", steps = " << r.steps);
        CHECK(worst < 1e-5);
    }
}

TEST_CASE("Robertson stiff problem: Rosenbrock reaches t = 1e4 within the step budget") {
    SolverOptions o = with(Method::Rosenbrock, 1e-6, 1e-12);
    const SimResult r = run(kRobertson, 1e4, 1001, o);
    REQUIRE_MESSAGE(r.ok(), to_string(r.status));
    CHECK(r.steps < o.max_steps);
    MESSAGE("Robertson: " << r.steps << " accepted, " << r.rejected << " rejected steps");

    double worst = 0.0;
    for (int i = 0; i < r.n_samples(); ++i) worst = std::max(worst, std::abs(r.at(i, 0) + r.at(i, 1) + r.at(i, 2) - 1.0));
    CHECK(worst < 1e-6);

    // Reference at t = 1e4 from SciPy solve_ivp (Radau and LSODA agree to 1e-10),
    // rtol = 1e-12, atol = 1e-16.
    CHECK(std::abs(r.at(1000, 0) - 1.073004285e-01) / 1.073004285e-01 < 1e-3);
    CHECK(std::abs(r.at(1000, 1) - 4.800166973e-07) / 4.800166973e-07 < 1e-3);
    CHECK(std::abs(r.at(1000, 2) - 8.926990914e-01) / 8.926990914e-01 < 1e-3);
}

TEST_CASE("Auto policy switches to the stiff solver on Robertson") {
    const SimResult r = run(kRobertson, 1e4, 101, with(Method::Auto, 1e-6, 1e-12));
    REQUIRE_MESSAGE(r.ok(), to_string(r.status));
    CHECK(r.used_stiff);
    double worst = 0.0;
    for (int i = 0; i < r.n_samples(); ++i) worst = std::max(worst, std::abs(r.at(i, 0) + r.at(i, 1) + r.at(i, 2) - 1.0));
    CHECK(worst < 1e-6);
}

TEST_CASE("Brusselator (A = 1, B = 3) settles to the reference period within 2%") {
    // Reference period 7.1569199212, from SciPy solve_ivp with DOP853 and with
    // Radau (rtol = 1e-12, atol = 1e-14, both agree to 1e-10): the mean
    // interval between upward crossings of x = 1 over t in [200, 400],
    // starting from (x, y) = (1, 1).
    constexpr double kReferencePeriod = 7.1569199212;
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    for (Method method : {Method::RK45, Method::Rosenbrock, Method::Auto}) {
        CAPTURE(name(method));
        const SimResult r = run(net, 200.0, 20001, with(method));
        REQUIRE_MESSAGE(r.ok(), to_string(r.status));
        const double period = mean_period(r, 0, 1.0, 100.0);
        MESSAGE(name(method) << ": period " << period << ", steps " << r.steps);
        CHECK(std::abs(period - kReferencePeriod) / kReferencePeriod < 0.02);
    }
}

TEST_CASE("an explosive network triggers the max_conc guard instead of hanging") {
    const Network net = parse_antimony("J0: S0 -> S0 + S0; k0*S0;\nk0 = 50;\nS0 = 1;\n");
    for (Method method : {Method::RK45, Method::Rosenbrock, Method::Auto}) {
        CAPTURE(name(method));
        const auto start = std::chrono::steady_clock::now();
        const SimResult r = run(net, 100.0, 1000, with(method));
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        CHECK(r.status == SimStatus::MaxConc);
        CHECK(r.n_samples() < 1000);
        CHECK(secs < 1.0);
    }
}

TEST_CASE("step budget and wall-clock guards stop a simulation") {
    SUBCASE("step budget") {
        SolverOptions o = with(Method::RK45, 1e-6, 1e-12);
        o.max_steps = 1000;
        const SimResult r = run(kRobertson, 1e4, 101, o);
        CHECK(r.status == SimStatus::StepBudget);
        CHECK(r.steps == 1000);
    }
    SUBCASE("wall clock") {
        SolverOptions o = with(Method::RK45, 1e-6, 1e-12);
        o.max_steps = 100000000;
        o.max_wall_seconds = 0.0;
        const SimResult r = run(kRobertson, 1e4, 101, o);
        CHECK(r.status == SimStatus::WallClock);
    }
    SUBCASE("concentration floor") {
        SolverOptions o = with(Method::RK45);
        o.min_conc = 0.5;  // the decay crosses 0.5 at t = 2 ln 2
        const SimResult r = run(kDecay, 10.0, 101, o);
        CHECK(r.status == SimStatus::Negative);
        CHECK(r.t.back() < 2 * std::log(2.0));
    }
}

TEST_CASE("stiff quasi-steady species: Rosenbrock samples do not overshoot") {
    // Found during M5. dS3/dt = k2*X1 - k3*S1*S3 relaxes at rate ~k3*S1 ~ 4e4
    // onto its quasi-steady value k2*X1/(k3*S1), which declines smoothly as S1
    // grows linearly. Cubic Hermite dense output turned this into spurious
    // oscillations of ~1e-6 (fitness 1.0 for a non-oscillating network).
    const Network net = parse_antimony(
        "const X0, X1;\n"
        "J0: X1 -> S1; k0*X1;\n"
        "J1: X0 -> S0 + X0; k1*X0;\n"
        "J2: X1 -> S3; k2*X1;\n"
        "J3: S1 + S3 -> S0 + S2; k3*S1*S3;\n"
        "J4: S0 + X1 -> X0; k4*S0*X1;\n"
        "J5: S0 + S3 -> S3 + X0; k5*S0*S3;\n"
        "k0 = 8.252664554556926; k1 = 0.023360102852485538; k2 = 0.09976968186966347; k3 = 50.44104264115237; "
        "k4 = 0.2525921537647396; k5 = 0.672908797472974;\n"
        "S0 = 1; S1 = 1; S2 = 1; S3 = 1;\n"
        "X0 = 1; X1 = 1;\n");
    for (Method method : {Method::Rosenbrock, Method::Auto}) {
        CAPTURE(name(method));
        const SimResult r = run(net, 100.0, 1000, with(method));
        REQUIRE_MESSAGE(r.ok(), to_string(r.status));
        // Over the second half, S3 tracks k2/(k3*S1) closely and decreases.
        for (int i = 500; i < r.n_samples(); ++i) {
            const double qss = 0.09976968186966347 / (50.44104264115237 * r.at(i, 1));
            REQUIRE_MESSAGE(std::abs(r.at(i, 3) - qss) < 0.01 * qss, "t = " << r.t[i] << ", S3 = " << r.at(i, 3)
                                                                             << ", quasi-steady " << qss);
            if (i > 500) REQUIRE(r.at(i, 3) <= r.at(i - 1, 3) * (1 + 1e-6));
        }
    }
}

TEST_CASE("simulation is deterministic") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const SimResult a = run(net, 50.0, 501, with(Method::Auto));
    const SimResult b = run(net, 50.0, 501, with(Method::Auto));
    CHECK(a.y == b.y);
    CHECK(a.steps == b.steps);
}
