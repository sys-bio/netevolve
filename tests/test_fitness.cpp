// M4: extrema detection and the peak-and-trough oscillation fitness.
#include <cmath>
#include <functional>
#include <vector>

#include "config.hpp"
#include "doctest/doctest.h"
#include "fitness.hpp"
#include "network.hpp"

using namespace ne;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Samples g(t) at n uniform points over [0, t_end], including both ends,
// matching the sample times used by simulate().
std::vector<double> sample(const std::function<double(double)>& g, double t_end, int n) {
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) y[i] = g(t_end * i / (n - 1));
    return y;
}

const OscillationParams kDefaults{};  // transient 0.5, n_target 5, min_rel_amplitude 0.01

}  // namespace

TEST_CASE("find_extrema: alternation and threshold") {
    // Up, down, a small ripple, then down further.
    const std::vector<double> y = {0, 1, 0, 0.5, 0.45, 0.48, 0, 1, 0};
    const auto ext = find_extrema(y, 0.1);
    REQUIRE(ext.size() == 5);
    CHECK(ext[0].peak);
    CHECK(ext[0].value == 1.0);
    CHECK_FALSE(ext[1].peak);
    CHECK(ext[2].peak);
    CHECK(ext[2].value == 0.5);  // the 0.45 / 0.48 ripple is merged into this peak
    CHECK_FALSE(ext[3].peak);
    CHECK(ext[3].value == 0.0);
    CHECK(ext[4].peak);
    for (std::size_t i = 1; i < ext.size(); ++i) CHECK(ext[i].peak != ext[i - 1].peak);
}

TEST_CASE("find_extrema: ripple at the window edges is not an extremum") {
    // Trough at index 2, peak at index 5, then a small wobble while still
    // falling at the end, plus a small wobble at the very start.
    const std::vector<double> y = {0.50, 0.52, 0.0, 0.5, 0.9, 1.0, 0.6, 0.3, 0.31, 0.2};
    const auto ext = find_extrema(y, 0.1);
    REQUIRE(ext.size() == 2);
    CHECK_FALSE(ext[0].peak);
    CHECK(ext[0].index == 2);
    CHECK(ext[1].peak);
    CHECK(ext[1].index == 5);

    // An end trough is kept when the signal afterwards clearly rises.
    const std::vector<double> z = {0.0, 1.0, 0.0, 0.5};
    CHECK(find_extrema(z, 0.1).size() == 2);
}

TEST_CASE("find_extrema: flat-topped peaks count once") {
    const std::vector<double> y = {0, 1, 1, 1, 0, 0, 1, 0};
    const auto ext = find_extrema(y, 0.1);
    REQUIRE(ext.size() == 3);
    CHECK(ext[0].index == 1);
    CHECK(ext[1].index == 4);
    CHECK(ext[2].index == 6);
}

TEST_CASE("sin(t) over 10 periods: n >= n_target, d ~ 0, fitness ~ 1") {
    const auto y = sample([](double t) { return std::sin(t); }, 20 * kPi, 2000);
    const SignalScore s = signal_score(y, kDefaults);
    CHECK(s.n_pairs >= kDefaults.n_target);
    CHECK(s.damping < 1e-3);
    CHECK(s.score > 0.99);
}

TEST_CASE("e^{-0.3 t} sin(t): fitness well below 1") {
    const auto y = sample([](double t) { return std::exp(-0.3 * t) * std::sin(t); }, 20 * kPi, 2000);
    const SignalScore s = signal_score(y, kDefaults);
    MESSAGE("damped sine: pairs " << s.n_pairs << ", damping " << s.damping << ", score " << s.score);
    CHECK(s.score < 0.5);
}

TEST_CASE("constant and monotone signals score 0") {
    SUBCASE("constant") {
        const std::vector<double> y(1000, 3.7);
        const SignalScore s = signal_score(y, kDefaults);
        CHECK(s.n_pairs == 0);
        CHECK(s.score == 0.0);
    }
    SUBCASE("increasing ramp") {
        const auto y = sample([](double t) { return 1 + 2 * t; }, 100, 1000);
        CHECK(signal_score(y, kDefaults).score == 0.0);
    }
    SUBCASE("decreasing ramp") {
        const auto y = sample([](double t) { return 500 - 2 * t; }, 100, 1000);
        CHECK(signal_score(y, kDefaults).score == 0.0);
    }
}

TEST_CASE("sin(t) with 0.1% noise: extrema are not double-counted") {
    const auto clean = sample([](double t) { return std::sin(t); }, 20 * kPi, 2000);
    const SignalScore s_clean = signal_score(clean, kDefaults);

    Rng rng(8);
    for (int trial = 0; trial < 20; ++trial) {
        std::vector<double> noisy = clean;
        for (double& v : noisy) v += 1e-3 * rng.uniform(-1.0, 1.0);  // +/- 0.1% of the unit amplitude
        const SignalScore s = signal_score(noisy, kDefaults);
        CAPTURE(trial);
        CHECK(s.n_pairs == s_clean.n_pairs);
        CHECK(s.score > 0.99);
    }
}

TEST_CASE("densely sampled sin(t) with 0.1% noise: many raw ripples, none counted") {
    // At 20,000 samples over 10 periods the signal moves at most ~0.003 per
    // sample, comparable to the noise, so the noise creates many spurious
    // local extrema. The amplitude filter must remove all of them.
    const auto clean = sample([](double t) { return std::sin(t); }, 20 * kPi, 20000);
    const int clean_pairs = signal_score(clean, kDefaults).n_pairs;
    Rng rng(9);
    for (int trial = 0; trial < 10; ++trial) {
        std::vector<double> noisy = clean;
        for (double& v : noisy) v += 1e-3 * rng.uniform(-1.0, 1.0);
        CAPTURE(trial);
        const SignalScore s = signal_score(noisy, kDefaults);
        CHECK(s.n_pairs == clean_pairs);
        CHECK(s.score > 0.99);

        const std::vector<double> w(noisy.begin() + static_cast<long>(noisy.size() / 2), noisy.end());
        const std::size_t raw = find_extrema(w, 0.0).size();
        const std::size_t filtered = find_extrema(w, 0.01 * 0.6366).size();  // 0.01 * mean|sin|
        CHECK(filtered == 2 * static_cast<std::size_t>(clean_pairs));
        CHECK(raw > 2 * filtered);
    }
}

TEST_CASE("examples/brusselator.ant: fitness >= 0.95") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const FitnessReport r = evaluate_oscillation(net, fitness_params(Config{}));
    REQUIRE(r.status == SimStatus::Ok);
    for (std::size_t i = 0; i < r.species.size(); ++i)
        MESSAGE("S" << i << ": pairs " << r.species[i].n_pairs << ", damping " << r.species[i].damping << ", score "
                    << r.species[i].score);
    CHECK(r.fitness >= 0.95);
}

TEST_CASE("examples/damped.ant (stable focus): fitness < 0.5") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/damped.ant");
    const FitnessReport r = evaluate_oscillation(net, fitness_params(Config{}));
    REQUIRE(r.status == SimStatus::Ok);
    MESSAGE("damped.ant fitness " << r.fitness);
    CHECK(r.fitness < 0.5);
}

TEST_CASE("a network that blows up has fitness exactly 0") {
    const Network net = parse_antimony("const X0;\nJ0: S0 -> S0 + S0; k0*S0;\nJ1: S0 -> X0; k1*S0;\nk0 = 5; k1 = 0.1;\nS0 = 1;\nX0 = 0;\n");
    const FitnessReport r = evaluate_oscillation(net, fitness_params(Config{}));
    CHECK(r.status == SimStatus::MaxConc);
    CHECK(r.fitness == 0.0);
    CHECK(oscillation_fitness(net, fitness_params(Config{})) == 0.0);
}

TEST_CASE("integrator noise around zero is not an oscillation") {
    // Found during M5: a random network in which nothing produces S2, so S2
    // decays to ~0 and then jitters at the 1e-9 level. Before the absolute
    // amplitude floor this scored a perfect 1.0.
    const Network net = parse_antimony(
        "const X0, X1;\n"
        "J0: S3 + X1 -> S0 + S0; k0*S3*X1;\n"
        "J1: S0 + S2 -> S3; k1*S0*S2;\n"
        "J2: S2 + S3 -> S1 + X1; k2*S2*S3;\n"
        "J3: S2 + X1 -> X0 + X1; k3*S2*X1;\n"
        "J4: S2 -> S1 + X0; k4*S2;\n"
        "J5: X0 + X0 -> S0 + X0; k5*X0*X0;\n"
        "k0 = 0.17457805103950375; k1 = 0.1185337363411992; k2 = 2.1633076150978705; k3 = 0.030354211136767922; "
        "k4 = 0.012140261307771828; k5 = 0.10919957822863836;\n"
        "S0 = 1; S1 = 1; S2 = 1; S3 = 1;\n"
        "X0 = 1; X1 = 1;\n");
    const FitnessReport r = evaluate_oscillation(net, fitness_params(Config{}));
    REQUIRE(r.status == SimStatus::Ok);
    CHECK(r.species[2].n_pairs == 0);
    CHECK(r.fitness == 0.0);

    // The same holds for a synthetic signal: tiny noise around a decayed value.
    Rng rng(4);
    std::vector<double> y(1000);
    for (double& v : y) v = 1e-9 * rng.uniform(-1.0, 1.0);
    CHECK(signal_score(y, kDefaults).score == 0.0);
}

TEST_CASE("fitness respects n_target and transient_fraction") {
    const auto y = sample([](double t) { return 2 + std::sin(t); }, 20 * kPi, 2000);
    OscillationParams p;
    p.n_target = 20;  // only ~5 pairs are available after the transient
    const SignalScore s = signal_score(y, p);
    CHECK(s.score == doctest::Approx(s.n_pairs / 20.0));
    p.transient_fraction = 0.0;  // all 10 periods count
    CHECK(signal_score(y, p).n_pairs >= 9);
}

TEST_CASE("tight-tolerance re-check: an integrator artifact loses its fitness") {
    // The artifact fools the in-house solvers at the default tolerance, so pin
    // solver = auto. (CVODE already rejects it at the default tolerance; see
    // test_cvode.cpp.) The concentration floor also rejects it on its own (see
    // the floor tests below), so turn the floor off to test the re-check alone.
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/integrator_artifact.ant");
    Config cfg;
    cfg.solver = "auto";
    cfg.osc_conc_floor = 0;
    const FitnessParams fp = fitness_params(cfg);
    const FitnessReport raw = evaluate_oscillation(net, fp);  // default tolerance only
    MESSAGE("default-tolerance fitness " << raw.fitness);
    REQUIRE(raw.fitness >= 0.95);  // looks like a perfect oscillator ...
    FitnessParams tight = fp;
    tight.solver = fp.recheck_solver;
    MESSAGE("tight-tolerance fitness " << evaluate_oscillation(net, tight).fitness);
    CHECK(oscillation_fitness(net, fp) < 0.5);  // ... but the re-check rejects it

    FitnessParams off = fp;
    off.recheck_threshold = 1.5;  // re-check off: the raw score comes through
    CHECK(oscillation_fitness(net, off) == raw.fitness);
}

TEST_CASE("tight-tolerance re-check: genuine and low-scoring networks are unaffected") {
    const FitnessParams fp = fitness_params(Config{});
    const Network bruss = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const double raw = evaluate_oscillation(bruss, fp).fitness;
    REQUIRE(raw >= fp.recheck_threshold);  // so the re-check runs
    CHECK(oscillation_fitness(bruss, fp) == doctest::Approx(raw).epsilon(0.01));

    const Network damped = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/damped.ant");
    const double raw_d = evaluate_oscillation(damped, fp).fitness;
    REQUIRE(raw_d < fp.recheck_threshold);  // below the threshold: no re-check
    CHECK(oscillation_fitness(damped, fp) == raw_d);
}

TEST_CASE("tight-tolerance re-check: config keys reach the fitness parameters") {
    Config c;
    c.recheck_threshold = 0.5;
    c.recheck_rtol = 1e-9;
    c.recheck_atol = 1e-12;
    c.max_conc = 1e5;
    const FitnessParams fp = fitness_params(c);
    CHECK(fp.recheck_threshold == 0.5);
    CHECK(fp.recheck_solver.rtol == 1e-9);
    CHECK(fp.recheck_solver.atol == 1e-12);
    CHECK(fp.recheck_solver.max_conc == 1e5);
}

TEST_CASE("solver and max_steps config keys reach both the normal and the re-check solver") {
    Config c;
    // Default: CVODE when the build has it, otherwise the in-house Auto policy.
    CHECK(fitness_params(c).solver.method == (cvode_available() ? Method::CVODE : Method::Auto));
    CHECK(fitness_params(c).solver.max_steps == 50000);
    c.solver = "rk45";
    c.max_steps = 500000;
    FitnessParams fp = fitness_params(c);
    CHECK(fp.solver.method == Method::RK45);
    CHECK(fp.recheck_solver.method == Method::RK45);
    CHECK(fp.solver.max_steps == 500000);
    c.solver = "rosenbrock";
    fp = fitness_params(c);
    CHECK(fp.solver.method == Method::Rosenbrock);
    CHECK(fp.recheck_solver.method == Method::Rosenbrock);

    // The Brusselator is a genuine oscillator with every solver choice.
    const Network bruss = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    for (const char* s : {"auto", "rk45", "rosenbrock"}) {
        c.solver = s;
        CHECK(oscillation_fitness(bruss, fitness_params(c)) >= 0.95);
    }
}

// --- Concentration floor ----------------------------------------------------------

namespace {

// A synthetic SimResult: species columns given as functions of t on [0, 100].
SimResult synthetic(const std::vector<std::function<double(double)>>& cols) {
    SimResult r;
    r.n_state = static_cast<int>(cols.size());
    for (int i = 0; i <= 1000; ++i) {
        const double t = 0.1 * i;
        r.t.push_back(t);
        for (const auto& f : cols) r.y.push_back(f(t));
    }
    return r;
}

}  // namespace

TEST_CASE("concentration floor: an oscillating species dipping below it zeroes the fitness") {
    OscillationParams p;  // default floor 1e-12
    const auto healthy = [](double t) { return 2.0 + std::sin(t); };
    const auto deep = [](double t) { return 1e-15 + std::pow(0.5 * (1.0 + std::sin(t)), 8); };  // min 1e-15
    const auto decaying = [](double t) { return std::exp(-t); };  // reaches 4e-44, never oscillates

    // A healthy oscillator next to a species that just decays to ~0: no effect.
    const FitnessReport ok = oscillation_report(synthetic({healthy, decaying}), p);
    CHECK(ok.fitness > 0.99);
    CHECK(ok.floor_species == -1);

    // An oscillating species that dips to 1e-15: the whole network scores 0,
    // even though the other species is a clean oscillator.
    const FitnessReport bad = oscillation_report(synthetic({healthy, deep}), p);
    REQUIRE(bad.species.size() == 2);
    CHECK(bad.species[1].n_pairs > 0);
    CHECK(bad.floor_species == 1);
    CHECK(bad.floor_min < 1e-12);
    CHECK(bad.fitness == 0.0);

    // Floor off: the clean oscillator's score comes through again.
    p.conc_floor = 0;
    CHECK(oscillation_report(synthetic({healthy, deep}), p).fitness > 0.99);
}

TEST_CASE("concentration floor: config key, and the integrator artifact is caught by it alone") {
    Config c;
    CHECK(fitness_params(c).osc.conc_floor == 1e-12);
    c.osc_conc_floor = 1e-9;
    CHECK(fitness_params(c).osc.conc_floor == 1e-9);
    c.osc_conc_floor = -1;
    CHECK_THROWS_AS(validate_config(c), ConfigError);

    // integrator_artifact.ant fools solver = auto at the default tolerance, but
    // its oscillating species dip below 1e-12, so the floor rejects it even
    // with the re-check turned off.
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/integrator_artifact.ant");
    Config cfg;
    cfg.solver = "auto";
    cfg.recheck_threshold = 1.5;  // re-check off
    const FitnessReport r = evaluate_oscillation(net, fitness_params(cfg));
    CHECK(r.floor_species >= 0);
    CHECK(r.fitness == 0.0);
    cfg.osc_conc_floor = 0;
    CHECK(evaluate_oscillation(net, fitness_params(cfg)).fitness >= 0.95);  // what it scored before the floor
}

TEST_CASE("concentration floor: the genuine Brusselator is unaffected") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const FitnessReport r = evaluate_oscillation(net, fitness_params(Config{}));
    CHECK(r.floor_species == -1);
    CHECK(r.fitness >= 0.95);
}
