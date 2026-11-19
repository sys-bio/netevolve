// M4: extrema detection and the peak-and-trough oscillation fitness.
#include <cmath>
#include <functional>
#include <vector>

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

TEST_CASE("fitness respects n_target and transient_fraction") {
    const auto y = sample([](double t) { return 2 + std::sin(t); }, 20 * kPi, 2000);
    OscillationParams p;
    p.n_target = 20;  // only ~5 pairs are available after the transient
    const SignalScore s = signal_score(y, p);
    CHECK(s.score == doctest::Approx(s.n_pairs / 20.0));
    p.transient_fraction = 0.0;  // all 10 periods count
    CHECK(signal_score(y, p).n_pairs >= 9);
}
