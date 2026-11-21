// Optional CVODE solver (built with NETEVOLVE_WITH_CVODE=ON). Without it, only
// checks that asking for CVODE fails cleanly.
#include <chrono>
#include <cmath>
#include <vector>

#include "config.hpp"
#include "doctest/doctest.h"
#include "fitness.hpp"
#include "integrator.hpp"
#include "model.hpp"
#include "network.hpp"

using namespace ne;

#ifndef NETEVOLVE_HAVE_CVODE

TEST_CASE("CVODE: absent from this build, and asking for it fails cleanly") {
    CHECK_FALSE(cvode_available());
    Config c;
    c.solver = "cvode";
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    const Network net = parse_antimony("const X0;\nJ0: S0 -> X0; k0*S0;\nk0 = 0.5;\nS0 = 1;\nX0 = 0;\n");
    SolverOptions o;
    o.method = Method::CVODE;
    CHECK_THROWS(simulate(make_model(net), net.init, 1.0, 11, o));
}

#else

namespace {

SimResult run(const Network& net, double t_end, int n_points, double rtol = 1e-6, double atol = 1e-9) {
    SolverOptions o;
    o.method = Method::CVODE;
    o.rtol = rtol;
    o.atol = atol;
    return simulate(make_model(net), net.init, t_end, n_points, o);
}

}  // namespace

TEST_CASE("CVODE: available, and selectable from the config") {
    CHECK(cvode_available());
    Config c;
    c.solver = "cvode";
    CHECK_NOTHROW(validate_config(c));
    CHECK(fitness_params(c).solver.method == Method::CVODE);
    CHECK(fitness_params(c).recheck_solver.method == Method::CVODE);
}

TEST_CASE("CVODE: dy/dt = -k y matches the analytic solution at every sample") {
    const Network net = parse_antimony("const X0;\nJ0: S0 -> X0; k0*S0;\nk0 = 0.5;\nS0 = 1;\nX0 = 0;\n");
    const SimResult r = run(net, 10.0, 101, 1e-10, 1e-14);
    REQUIRE_MESSAGE(r.ok(), to_string(r.status));
    REQUIRE(r.n_samples() == 101);
    CHECK(r.t.back() == 10.0);
    double worst = 0.0;
    for (int i = 0; i < r.n_samples(); ++i) {
        const double e = std::exp(-0.5 * r.t[i]);
        worst = std::max(worst, std::abs(r.at(i, 0) - e) / e);
    }
    MESSAGE("CVODE decay: worst relative error " << worst << ", steps " << r.steps);
    CHECK(worst < 1e-6);
}

TEST_CASE("CVODE: Robertson stiff problem reaches t = 1e4 with mass conserved") {
    const Network net = parse_antimony(
        "J0: S0 -> S1; k0*S0;\n"
        "J1: S1 + S2 -> S0 + S2; k1*S1*S2;\n"
        "J2: S1 + S1 -> S1 + S2; k2*S1*S1;\n"
        "k0 = 0.04; k1 = 1e4; k2 = 3e7;\n"
        "S0 = 1; S1 = 0; S2 = 0;\n");
    const SimResult r = run(net, 1e4, 1001, 1e-6, 1e-12);
    REQUIRE_MESSAGE(r.ok(), to_string(r.status));
    MESSAGE("CVODE Robertson: " << r.steps << " steps");
    double worst = 0.0;
    for (int i = 0; i < r.n_samples(); ++i) worst = std::max(worst, std::abs(r.at(i, 0) + r.at(i, 1) + r.at(i, 2) - 1.0));
    CHECK(worst < 1e-6);
    // Same SciPy reference as the Rosenbrock test.
    CHECK(std::abs(r.at(1000, 0) - 1.073004285e-01) / 1.073004285e-01 < 1e-3);
    CHECK(std::abs(r.at(1000, 1) - 4.800166973e-07) / 4.800166973e-07 < 1e-3);
    CHECK(std::abs(r.at(1000, 2) - 8.926990914e-01) / 8.926990914e-01 < 1e-3);
}

TEST_CASE("CVODE: Brusselator agrees with RK45 and is a genuine oscillator") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const SimResult c = run(net, 100.0, 1001, 1e-10, 1e-14);
    SolverOptions o;
    o.method = Method::RK45;
    o.rtol = 1e-10;
    o.atol = 1e-14;
    const SimResult k = simulate(make_model(net), net.init, 100.0, 1001, o);
    REQUIRE(c.ok());
    REQUIRE(k.ok());
    double worst = 0.0;
    for (std::size_t i = 0; i < c.y.size(); ++i) worst = std::max(worst, std::abs(c.y[i] - k.y[i]));
    MESSAGE("CVODE vs RK45 on the Brusselator: max abs difference " << worst);
    CHECK(worst < 1e-5);

    Config cfg;
    cfg.solver = "cvode";
    CHECK(oscillation_fitness(net, fitness_params(cfg)) >= 0.95);
}

TEST_CASE("CVODE: guards and limits") {
    SUBCASE("max_conc") {
        const Network net = parse_antimony("J0: S0 -> S0 + S0; k0*S0;\nk0 = 50;\nS0 = 1;\n");
        const auto start = std::chrono::steady_clock::now();
        const SimResult r = run(net, 100.0, 1000);
        CHECK(r.status == SimStatus::MaxConc);
        CHECK(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() < 1.0);
    }
    SUBCASE("step budget") {
        const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
        SolverOptions o;
        o.method = Method::CVODE;
        o.max_steps = 50;
        const SimResult r = simulate(make_model(net), net.init, 100.0, 101, o);
        CHECK(r.status == SimStatus::StepBudget);
        CHECK(r.steps == 50);
    }
    SUBCASE("concentration floor") {
        const Network net = parse_antimony("const X0;\nJ0: S0 -> X0; k0*S0;\nk0 = 0.5;\nS0 = 1;\nX0 = 0;\n");
        SolverOptions o;
        o.method = Method::CVODE;
        o.min_conc = 0.5;
        const SimResult r = simulate(make_model(net), net.init, 10.0, 101, o);
        CHECK(r.status == SimStatus::Negative);
    }
}

TEST_CASE("CVODE: deterministic, and the integrator artifact is rejected") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const SimResult a = run(net, 50.0, 501), b = run(net, 50.0, 501);
    CHECK(a.y == b.y);
    CHECK(a.steps == b.steps);

    const Network bad = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/integrator_artifact.ant");
    Config cfg;
    cfg.solver = "cvode";
    const double f = oscillation_fitness(bad, fitness_params(cfg));
    MESSAGE("integrator_artifact.ant with CVODE: fitness " << f);
    CHECK(f < 0.5);
}

TEST_CASE("CVODE: the deep-dip oscillator is rejected by the concentration floor") {
    // examples/deep_dip.ant scores 1.0 with CVODE when nothing guards against
    // its species dipping to about 1e-22; the floor (default 1e-12) rejects it.
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/deep_dip.ant");
    Config cfg;
    cfg.solver = "cvode";
    const FitnessReport r = evaluate_oscillation(net, fitness_params(cfg));
    MESSAGE("deep_dip.ant: floor species " << r.floor_species << ", minimum " << r.floor_min);
    CHECK(r.floor_species >= 0);
    CHECK(network_fitness(net, fitness_params(cfg)) == 0.0);

    cfg.osc_conc_floor = 0;  // without the floor it gets through
    CHECK(network_fitness(net, fitness_params(cfg)) >= 0.95);
}

#endif
