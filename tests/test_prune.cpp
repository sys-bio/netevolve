// M7: prune and the long-run check.
#include <string>

#include "doctest/doctest.h"
#include "evolve.hpp"
#include "fitness.hpp"
#include "network.hpp"

using namespace ne;

namespace {

const FitnessParams kFp = fitness_params(Config{});

FitnessFn osc() {
    return [](const Network& n) { return oscillation_fitness(n, kFp); };
}

// The Brusselator of examples/brusselator.ant with B = X1 and the initial
// state (x0, y0) set by the caller.
Network brusselator(const std::string& B, const std::string& x0, const std::string& y0) {
    return parse_antimony(
        "const X0, X1, X2;\n"
        "J0: X0 -> S0; k0*X0;\n"
        "J1: 2 S0 + S1 -> 3 S0; k1*S0*S0*S1;\n"
        "J2: X1 + S0 -> S1 + X2; k2*X1*S0;\n"
        "J3: S0 -> X2; k3*S0;\n"
        "k0 = 1; k1 = 1; k2 = 1; k3 = 1;\n"
        "S0 = " + x0 + "; S1 = " + y0 + ";\n"
        "X0 = 1; X1 = " + B + "; X2 = 0;\n");
}

}  // namespace

TEST_CASE("prune removes dispensable reactions and keeps the core oscillator") {
    // The Brusselator plus two reactions with negligible rates, placed
    // between and after the core reactions.
    const Network padded = parse_antimony(
        "const X0, X1, X2;\n"
        "J0: X0 -> S0; k0*X0;\n"
        "J1: 2 S0 + S1 -> 3 S0; k1*S0*S0*S1;\n"
        "J2: S1 -> X2; k2*S1;\n"
        "J3: X1 + S0 -> S1 + X2; k3*X1*S0;\n"
        "J4: S0 -> X2; k4*S0;\n"
        "J5: S0 + S1 -> S1 + S1; k5*S0*S1;\n"
        "k0 = 1; k1 = 1; k2 = 1e-6; k3 = 1; k4 = 1; k5 = 1e-6;\n"
        "S0 = 1; S1 = 1;\n"
        "X0 = 1; X1 = 3; X2 = 0;\n");
    const Network core = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");

    const PruneResult pr = prune(padded, osc(), 1);
    CHECK(pr.removed == 2);
    REQUIRE(pr.net.reactions.size() == 4);
    for (std::size_t j = 0; j < 4; ++j) {
        CAPTURE(j);
        CHECK(pr.net.reactions[j] == core.reactions[j]);
    }
    CHECK(pr.original_fitness >= 0.95);
    CHECK(pr.fitness >= 0.95 * pr.original_fitness);
    CHECK(pr.fitness == oscillation_fitness(pr.net, kFp));
}

TEST_CASE("prune keeps every reaction of a minimal oscillator") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const PruneResult pr = prune(net, osc(), 1);
    CHECK(pr.removed == 0);
    CHECK(pr.net == net);
}

TEST_CASE("prune never goes below r_min") {
    // Fitness 0, so every removal is accepted until r_min stops it.
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/damped.ant");
    CHECK(prune(net, osc(), 3).net.reactions.size() == 3);
    CHECK(prune(net, osc(), 4).net.reactions.size() == 4);
    CHECK(prune(net, osc(), 1).net.reactions.size() == 1);
}

TEST_CASE("prune with a toy fitness: the keep threshold is relative to the original") {
    // Fitness = sum of k / 100. Removing a reaction is kept only while the
    // total stays >= 95% of the original total.
    Network net = parse_antimony(
        "const X0;\n"
        "J0: X0 -> S0; k0*X0;\nJ1: S0 -> S1; k1*S0;\nJ2: S1 -> X0; k2*S1;\nJ3: S1 -> S0; k3*S1;\n"
        "k0 = 50; k1 = 46; k2 = 3; k3 = 1;\nS0 = 1; S1 = 1;\nX0 = 1;\n");
    auto sum_k = [](const Network& n) {
        double s = 0;
        for (const Reaction& r : n.reactions) s += r.k;
        return s / 100;
    };
    const PruneResult pr = prune(net, sum_k, 1);
    // Original 1.0. Removing J2 (k = 3) leaves 0.97, then J3 (k = 1) leaves 0.96;
    // removing anything else would drop below 0.95.
    CHECK(pr.removed == 2);
    REQUIRE(pr.net.reactions.size() == 2);
    CHECK(pr.net.reactions[0].k == 50);
    CHECK(pr.net.reactions[1].k == 46);
    CHECK(pr.fitness == doctest::Approx(0.96));
}

TEST_CASE("long-run check: sustained oscillation passes") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const LongRunReport lr = long_run_check(net, kFp);
    MESSAGE("brusselator: amplitude " << lr.amp_ref << " -> " << lr.amp_last << " (" << 100 * lr.rel_change << "%)");
    CHECK(lr.status == SimStatus::Ok);
    CHECK(lr.passed);
    CHECK(lr.amp_ref > 1.0);
}

TEST_CASE("long-run check: slow damping that fitness misses is caught") {
    // Just past the Hopf bifurcation (B < 1 + A^2 = 2) the focus is stable,
    // with linear decay rate (2 - B) / 2 = 0.0005. Starting close to the
    // fixed point (1, B) keeps the orbit small, so the decay stays close to
    // that slow linear rate: the fitness window [50, 100] shows little damping
    // (fitness 0.97), but by t = 500 the amplitude has fallen by about 25%.
    // (Starting far away, e.g. at (1, 1), the cubic term dominates and the
    // damping is already visible in the fitness: about 0.77.)
    const Network net = brusselator("1.999", "1.02", "1.999");
    const double f = oscillation_fitness(net, kFp);
    const LongRunReport lr = long_run_check(net, kFp);
    MESSAGE("B = 1.999 from (1.02, 1.999): fitness " << f << ", amplitude " << lr.amp_ref << " -> " << lr.amp_last
                                                     << " (" << 100 * lr.rel_change << "%)");
    CHECK(f >= 0.95);  // looks like a good oscillator to the fitness function
    CHECK(lr.status == SimStatus::Ok);
    CHECK_FALSE(lr.passed);
    CHECK(lr.rel_change < -0.10);
}

TEST_CASE("long-run check: non-oscillating and failing networks do not pass") {
    SUBCASE("stable focus") {
        const LongRunReport lr = long_run_check(load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/damped.ant"), kFp);
        CHECK_FALSE(lr.passed);
        CHECK(lr.amp_ref == 0.0);
    }
    SUBCASE("blow-up") {
        const Network net = parse_antimony("const X0;\nJ0: S0 -> S0 + S0; k0*S0;\nJ1: S0 -> X0; k1*S0;\nk0 = 5; k1 = 0.1;\nS0 = 1;\nX0 = 0;\n");
        const LongRunReport lr = long_run_check(net, kFp);
        CHECK_FALSE(lr.passed);
        CHECK(lr.status == SimStatus::MaxConc);
    }
}

TEST_CASE("window_amplitude measures mean peak-to-trough difference") {
    std::vector<double> y;
    for (int i = 0; i < 2000; ++i) y.push_back(5 + 2 * std::sin(i * 0.05));
    CHECK(window_amplitude(y, OscillationParams{}) == doctest::Approx(4.0).epsilon(0.001));
    CHECK(window_amplitude(std::vector<double>(100, 1.0), OscillationParams{}) == 0.0);
}
