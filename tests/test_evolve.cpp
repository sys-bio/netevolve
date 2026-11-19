// M5: mutation operators, selection with elitism, the generation loop and
// the evolve command outputs.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "config.hpp"
#include "doctest/doctest.h"
#include "evolve.hpp"
#include "fitness.hpp"
#include "network.hpp"

using namespace ne;
namespace fs = std::filesystem;

namespace {

std::string read_file(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    REQUIRE_MESSAGE(in.good(), "cannot read " << p.string());
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

fs::path fresh_dir(const std::string& name) {
    const fs::path d = fs::temp_directory_path() / ("netevolve_test_" + name);
    fs::remove_all(d);
    return d;
}

// A random network with a random shape and reaction count in [r_min, r_max].
Network random_test_network(Rng& rng, const Config& cfg) {
    GenParams p = gen_params(cfg);
    p.n_floating = rng.uniform_int(2, 6);
    p.n_boundary = rng.uniform_int(0, 3);
    p.r_initial = rng.uniform_int(cfg.r_min, cfg.r_max);
    return random_network(rng, p);
}

// Number of molecules on each side of each reaction (the reaction "type").
std::vector<std::pair<int, int>> shape(const Network& net) {
    std::vector<std::pair<int, int>> out;
    for (const Reaction& r : net.reactions) {
        int a = 0, b = 0;
        for (const Term& t : r.reactants) a += t.coeff;
        for (const Term& t : r.products) b += t.coeff;
        out.emplace_back(a, b);
    }
    return out;
}

// A small, fast configuration for runs that use the real oscillation fitness.
Config small_config(const std::string& out) {
    Config c;
    c.population = 24;
    c.generations = 8;
    c.output_dir = out;
    c.seed = 7;
    return c;
}

}  // namespace

TEST_CASE("every mutation operator, applied 10,000 times, yields valid networks") {
    const Config cfg;
    for (MutationKind kind : {MutationKind::Rate, MutationKind::Add, MutationKind::Delete, MutationKind::Rewire}) {
        CAPTURE(static_cast<int>(kind));
        Rng rng(1000 + static_cast<int>(kind));
        int applied = 0;
        for (int i = 0; i < 10000; ++i) {
            Network net = random_test_network(rng, cfg);
            const Network before = net;
            const bool ok = apply_mutation(net, kind, rng, cfg);
            applied += ok;
            if (!ok) REQUIRE(net == before);  // a failed mutation leaves the network unchanged
            const std::string err = check_validity(net, cfg.r_min, cfg.r_max);
            REQUIRE_MESSAGE(err.empty(), "iteration " << i << ": " << err << "\n" << write_antimony(net));
            for (const Reaction& r : net.reactions) {
                REQUIRE(r.k >= cfg.k_min);
                REQUIRE(r.k <= cfg.k_max);
            }
        }
        MESSAGE("operator " << static_cast<int>(kind) << " applied " << applied << " / 10000 times");
        CHECK(applied > 5000);  // each operator is genuinely exercised
    }
}

TEST_CASE("make_child, applied 10,000 times, yields valid networks") {
    const Config cfg;
    Rng rng(77);
    int changed = 0;
    for (int i = 0; i < 10000; ++i) {
        const Network parent = random_test_network(rng, cfg);
        const Network child = make_child(parent, rng, cfg);
        REQUIRE(check_validity(child, cfg.r_min, cfg.r_max).empty());
        changed += !(child == parent);
    }
    CHECK(changed > 9900);
}

TEST_CASE("operator specifics") {
    const Config cfg;
    Rng rng(5);
    SUBCASE("Add and Delete change the reaction count by one, Rate and Rewire keep it") {
        for (int i = 0; i < 1000; ++i) {
            const Network net = random_test_network(rng, cfg);
            const std::size_t R = net.reactions.size();
            Network a = net, d = net, r = net, w = net;
            if (apply_mutation(a, MutationKind::Add, rng, cfg)) REQUIRE(a.reactions.size() == R + 1);
            if (apply_mutation(d, MutationKind::Delete, rng, cfg)) REQUIRE(d.reactions.size() == R - 1);
            if (apply_mutation(r, MutationKind::Rate, rng, cfg)) REQUIRE(r.reactions.size() == R);
            if (apply_mutation(w, MutationKind::Rewire, rng, cfg)) {
                REQUIRE(w.reactions.size() == R);
                REQUIRE(shape(w) == shape(net));  // rewiring keeps every reaction's type
            }
        }
    }
    SUBCASE("Add is refused at r_max and Delete at r_min") {
        GenParams p = gen_params(cfg);
        p.r_initial = cfg.r_max;
        Network full = random_network(rng, p);
        CHECK_FALSE(apply_mutation(full, MutationKind::Add, rng, cfg));
        p.r_initial = cfg.r_min;
        Network minimal = random_network(rng, p);
        CHECK_FALSE(apply_mutation(minimal, MutationKind::Delete, rng, cfg));
    }
}

TEST_CASE("rate mutation keeps k in [k_min, k_max]") {
    Config cfg;
    cfg.mut_sigma = 3.0;  // large steps, so the bounds are hit often
    Rng rng(9);
    Network net = random_network(rng, gen_params(cfg));
    bool hit_min = false, hit_max = false;
    for (int i = 0; i < 10000; ++i) {
        REQUIRE(apply_mutation(net, MutationKind::Rate, rng, cfg));
        for (const Reaction& r : net.reactions) {
            REQUIRE(r.k >= cfg.k_min);
            REQUIRE(r.k <= cfg.k_max);
            hit_min |= r.k == cfg.k_min;
            hit_max |= r.k == cfg.k_max;
        }
    }
    CHECK(hit_min);
    CHECK(hit_max);
}

TEST_CASE("mutation kinds are chosen with the configured probabilities") {
    Config cfg;
    Rng rng(3);
    int counts[4] = {0, 0, 0, 0};
    const int N = 200000;
    for (int i = 0; i < N; ++i) ++counts[static_cast<int>(choose_mutation(rng, cfg))];
    CHECK(counts[0] / double(N) == doctest::Approx(0.70).epsilon(0.02));
    CHECK(counts[1] / double(N) == doctest::Approx(0.15).epsilon(0.03));
    CHECK(counts[2] / double(N) == doctest::Approx(0.10).epsilon(0.03));
    CHECK(counts[3] / double(N) == doctest::Approx(0.05).epsilon(0.05));
    cfg.p_rewire = 0;
    for (int i = 0; i < 10000; ++i) REQUIRE(choose_mutation(rng, cfg) != MutationKind::Rewire);
}

TEST_CASE("toy fitness (exactly 7 reactions) is solved within 50 generations") {
    Config cfg;
    cfg.population = 50;
    cfg.generations = 50;
    cfg.r_initial = 3;  // needs four successful Add mutations
    cfg.target_fitness = 1.0;
    auto toy = [](const Network& n) {
        const int R = static_cast<int>(n.reactions.size());
        return 1.0 / (1.0 + std::abs(R - 7));
    };
    const EvolutionResult res = evolve(cfg, toy);
    CHECK(res.best_fitness == 1.0);
    CHECK(res.best.reactions.size() == 7);
    CHECK(res.log.back().generation <= 50);
    CHECK(res.log.back().best == 1.0);  // stopped because the target was reached
    MESSAGE("toy problem solved at generation " << res.best_generation);
}

TEST_CASE("elitism: best fitness is non-decreasing across generations") {
    SUBCASE("oscillation fitness") {
        Config cfg = small_config("unused");
        cfg.generations = 15;
        const FitnessParams fp = fitness_params(cfg);
        const EvolutionResult res = evolve(cfg, [&](const Network& n) { return oscillation_fitness(n, fp); });
        for (std::size_t g = 1; g < res.log.size(); ++g) REQUIRE(res.log[g].best >= res.log[g - 1].best);
        MESSAGE("oscillation run: best " << res.log.front().best << " -> " << res.log.back().best);
    }
    SUBCASE("noisy-looking toy fitness") {
        // A fitness with many local optima, so the population churns.
        Config cfg;
        cfg.population = 40;
        cfg.generations = 60;
        cfg.target_fitness = 2.0;  // never reached
        auto bumpy = [](const Network& n) {
            double s = 0;
            for (const Reaction& r : n.reactions) s += std::sin(7.0 * std::log(r.k));
            return 0.5 + 0.5 * std::tanh(s / 4);
        };
        const EvolutionResult res = evolve(cfg, bumpy);
        REQUIRE(res.log.size() == 61);
        for (std::size_t g = 1; g < res.log.size(); ++g) REQUIRE(res.log[g].best >= res.log[g - 1].best);
    }
}

TEST_CASE("determinism: two evolve runs with the same seed give byte-identical outputs") {
    const fs::path a = fresh_dir("det_a"), b = fresh_dir("det_b");
    run_evolve_command(small_config(a.string()));
    run_evolve_command(small_config(b.string()));
    const std::string log_a = read_file(a / "log.csv");
    CHECK(log_a == read_file(b / "log.csv"));
    CHECK(read_file(a / "best_final.ant") == read_file(b / "best_final.ant"));
    CHECK(read_file(a / "best_final.csv") == read_file(b / "best_final.csv"));

    // A different seed gives a different run. (Compare the best network: with
    // a small population every oscillation fitness may be 0, which makes the
    // logs of different seeds identical too.)
    Config other = small_config(fresh_dir("det_c").string());
    other.seed = 8;
    run_evolve_command(other);
    CHECK(read_file(a / "best_final.ant") != read_file(fs::path(other.output_dir) / "best_final.ant"));
    fs::remove_all(a);
    fs::remove_all(b);
    fs::remove_all(other.output_dir);
}

TEST_CASE("determinism: log rows are identical for the same seed when fitness varies") {
    // Same check as above with a fitness that is non-zero for every network,
    // so the logged statistics are not trivially all zero.
    Config cfg;
    cfg.population = 60;
    cfg.generations = 40;
    cfg.target_fitness = 2.0;
    auto f = [](const Network& n) {
        double s = 0;
        for (const Reaction& r : n.reactions) s += std::sin(3.0 * std::log(r.k)) + 0.1 * r.reactants.size();
        return 0.5 + 0.5 * std::tanh(s / 5);
    };
    auto rows = [&](std::uint64_t seed) {
        Config c = cfg;
        c.seed = seed;
        std::string out;
        for (const GenerationStats& s : evolve(c, f).log) out += log_row(s) + "\n";
        return out;
    };
    const std::string first = rows(11);
    CHECK(first == rows(11));
    CHECK(first != rows(12));
    CHECK(first.find(",0,0,0,") == std::string::npos);  // statistics are not all zero
}

TEST_CASE("M6: log.csv is byte-identical for threads = 1 and threads = 8") {
    // Seed 27 with population 100 reaches a non-zero oscillation fitness
    // (0.257) by generation 4, so the logs carry real numbers.
    auto run_with = [](int threads, const std::string& name) {
        Config c;
        c.population = 100;
        c.generations = 4;
        c.seed = 27;
        c.threads = threads;
        c.output_dir = fresh_dir(name).string();
        run_evolve_command(c);
        return fs::path(c.output_dir);
    };
    const fs::path one = run_with(1, "threads_1"), eight = run_with(8, "threads_8"), all = run_with(0, "threads_0");
    const std::string log1 = read_file(one / "log.csv");
    CHECK(log1 == read_file(eight / "log.csv"));
    CHECK(log1 == read_file(all / "log.csv"));
    for (const char* f : {"best_final.ant", "best_final.csv"}) {
        CAPTURE(f);
        CHECK(read_file(one / f) == read_file(eight / f));
        CHECK(read_file(one / f) == read_file(all / f));
    }
    // The run is not trivial: some generation has a non-zero best fitness.
    const Network best = load_antimony_file((one / "best_final.ant").string());
    CHECK(oscillation_fitness(best, fitness_params(Config{})) > 0.0);
    for (const fs::path& d : {one, eight, all}) fs::remove_all(d);
}

TEST_CASE("M6: evolve() results do not depend on the thread count") {
    Config cfg;
    cfg.population = 120;
    cfg.generations = 30;
    cfg.target_fitness = 2.0;
    auto f = [](const Network& n) {
        double s = 0;
        for (const Reaction& r : n.reactions) s += std::sin(3.0 * std::log(r.k)) + 0.1 * r.reactants.size();
        return 0.5 + 0.5 * std::tanh(s / 5);
    };
    std::string reference;
    for (int threads : {1, 2, 3, 8, 17}) {
        CAPTURE(threads);
        Config c = cfg;
        c.threads = threads;
        const EvolutionResult res = evolve(c, f);
        std::string rows;
        for (const GenerationStats& s : res.log) rows += log_row(s) + "\n";
        rows += write_antimony(res.best);
        if (reference.empty()) reference = rows;
        CHECK(rows == reference);
        // 120 in generation 0, then 108 new children (12 elites) in each of 30 generations.
        CHECK(res.n_evaluations == 120 + 30 * 108);
    }
}

TEST_CASE("M6: an exception in a worker thread reaches the caller") {
    Config cfg;
    cfg.population = 50;
    cfg.generations = 3;
    cfg.threads = 4;
    auto bad = [](const Network& n) -> double {
        if (n.reactions.size() == 6 && n.reactions[0].k > 50) throw std::runtime_error("fitness failed");
        return 0.0;
    };
    CHECK_THROWS_WITH_AS(evolve(cfg, bad), "fitness failed", std::runtime_error);
}

TEST_CASE("evolve command writes the documented output files") {
    const fs::path d = fresh_dir("outputs");
    const Config cfg = small_config(d.string());
    const EvolutionResult res = run_evolve_command(cfg);

    // log.csv: header plus one row per generation.
    std::istringstream log(read_file(d / "log.csv"));
    std::string line;
    std::getline(log, line);
    CHECK(line == "generation,best,mean,median,best_reactions");
    int rows = 0;
    while (std::getline(log, line)) ++rows;
    CHECK(rows == static_cast<int>(res.log.size()));

    // run.cfg is the fully resolved config, including the seed.
    const Config reread = load_config_file((d / "run.cfg").string());
    CHECK(write_config(reread) == write_config(cfg));
    CHECK(reread.seed == 7);

    // best_final.ant loads and is the best network; best_genNNNN.ant exists for the best generation.
    const Network best = load_antimony_file((d / "best_final.ant").string());
    CHECK(best == res.best);
    char name[32];
    std::snprintf(name, sizeof name, "best_gen%04d.ant", res.best_generation);
    CHECK(fs::exists(d / name));
    CHECK(fs::exists(d / "best_gen0000.ant"));
    CHECK(read_file(d / "best_final.csv").rfind("time,S0,S1,S2,S3\n", 0) == 0);
    fs::remove_all(d);
}
