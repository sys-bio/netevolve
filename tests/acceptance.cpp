// M7 acceptance test (CTest label "slow"; excluded from the default run).
//
// With the default config and seeds 1-5: at least 3 of the 5 runs must
// reach fitness >= 0.95 within 500 generations, and the pruned result must
// pass the long-run check. Outputs for each seed are kept in
// <output base>/seedN for inspection, including pruned.ant.
//
// Runs only when the environment variable NETEVOLVE_RUN_SLOW is set to 1;
// otherwise it exits with code 77, which CTest reports as "skipped".
//
// Usage: netevolve_acceptance CONFIG [OUTPUT_BASE]

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "config.hpp"
#include "integrator.hpp"
#include "evolve.hpp"
#include "fitness.hpp"
#include "network.hpp"

int main(int argc, char** argv) {
    const char* run = std::getenv("NETEVOLVE_RUN_SLOW");
    if (!run || std::string(run) != "1") {
        std::cout << "acceptance test skipped (set NETEVOLVE_RUN_SLOW=1 to run it)\n";
        return 77;
    }
    if (argc < 2) {
        std::cerr << "usage: netevolve_acceptance CONFIG [OUTPUT_BASE]\n";
        return 2;
    }
    const std::string base = argc >= 3 ? argv[2] : "acceptance_out";

    try {
        const ne::Config defaults = ne::load_config_file(argv[1]);
        int passed = 0;
        std::printf("%-5s %-9s %-6s %-10s %-9s %-8s %-28s %s\n", "seed", "best", "gen", "evals", "pruned", "R", "long-run amplitude", "result");
        for (std::uint64_t seed = 1; seed <= 5; ++seed) {
            ne::Config cfg = defaults;
            cfg.seed = seed;
            cfg.output_dir = (std::filesystem::path(base) / ("seed" + std::to_string(seed))).string();
            std::filesystem::remove_all(cfg.output_dir);

            const ne::EvolutionResult res = ne::run_evolve_command(cfg);
            const ne::FitnessParams fp = ne::fitness_params(cfg);
            const ne::PruneResult pr = ne::prune(
                res.best, [&](const ne::Network& n) { return ne::oscillation_fitness(n, fp); }, cfg.r_min);
            const ne::LongRunReport lr = ne::long_run_check(pr.net, fp);

            char header[128];
            std::snprintf(header, sizeof header, "netevolve: seed %llu, pruned, fitness %.6g",
                          static_cast<unsigned long long>(seed), pr.fitness);
            std::ofstream(std::filesystem::path(cfg.output_dir) / "pruned.ant", std::ios::binary)
                << ne::write_antimony(pr.net, header);

            const bool reached = res.best_fitness >= 0.95;
            const bool ok = reached && lr.passed;
            passed += ok;
            char amp[64];
            if (lr.status != ne::SimStatus::Ok)
                std::snprintf(amp, sizeof amp, "sim: %s", ne::to_string(lr.status));
            else
                std::snprintf(amp, sizeof amp, "%.4g -> %.4g (%+.1f%%)", lr.amp_ref, lr.amp_last, 100 * lr.rel_change);
            std::printf("%-5llu %-9.4f %-6d %-10ld %-9.4f %zu->%-5zu %-28s %s\n", static_cast<unsigned long long>(seed),
                        res.best_fitness, res.best_generation, res.n_evaluations, pr.fitness, res.best.reactions.size(),
                        pr.net.reactions.size(), amp, ok ? "PASS" : (reached ? "FAIL (long-run)" : "FAIL (fitness)"));
            std::fflush(stdout);
        }
        std::printf("%d of 5 runs passed (need >= 3)\n", passed);
        return passed >= 3 ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
