// Mutation operators, selection with elitism, and the generation loop
// (SPEC §1.3), plus the file outputs of the `evolve` command (SPEC §2).
#pragma once

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

#include "config.hpp"
#include "network.hpp"
#include "rng.hpp"

namespace ne {

enum class MutationKind { Rate, Add, Delete, Rewire };

// Applies one mutation of the given kind in place. Returns false, leaving
// `net` unchanged, if the operator is not applicable (Add at r_max, Delete at
// r_min) or its random draw produced an invalid network.
bool apply_mutation(Network& net, MutationKind kind, Rng& rng, const Config& cfg);

// Chooses a mutation kind with probabilities p_rate : p_add : p_delete : p_rewire.
MutationKind choose_mutation(Rng& rng, const Config& cfg);

// A child of `parent` with exactly one mutation. Tries up to kMutationAttempts
// times (choosing the kind afresh each time); if every attempt fails the child
// is an unchanged copy of the parent.
constexpr int kMutationAttempts = 10;
Network make_child(const Network& parent, Rng& rng, const Config& cfg);

// Random generation parameters taken from the config.
GenParams gen_params(const Config& cfg);

struct GenerationStats {
    int generation = 0;
    double best = 0.0;
    double mean = 0.0;
    double median = 0.0;
    int best_reactions = 0;
};

struct EvolutionResult {
    Network best;
    double best_fitness = 0.0;
    int best_generation = 0;
    long n_evaluations = 0;  // fitness evaluations performed (cached elites excluded)
    std::vector<GenerationStats> log;
};

// Must be safe to call concurrently from several threads (cfg.threads > 1).
using FitnessFn = std::function<double(const Network&)>;

// Number of worker threads for `threads` (0 = all hardware threads).
int resolve_threads(int threads);

struct EvolutionCallbacks {
    // Called after each generation is evaluated.
    std::function<void(const GenerationStats&)> on_generation;
    // Called whenever the best fitness strictly improves.
    std::function<void(int generation, const Network& best, double fitness)> on_improvement;
};

// Runs the evolutionary algorithm. Generation 0 is the random initial
// population; each of the following cfg.generations generations is produced
// by elitism plus mutation. Stops early once best >= cfg.target_fitness.
// All random draws come from one Rng seeded with cfg.seed.
EvolutionResult evolve(const Config& cfg, const FitnessFn& fitness, const EvolutionCallbacks& cb = {});

// The `evolve` command: validates cfg, runs evolve() with the oscillation
// fitness, and writes log.csv, best_genNNNN.ant, best_final.ant,
// best_final.csv and run.cfg into cfg.output_dir. Progress goes to `log`
// if non-null.
EvolutionResult run_evolve_command(const Config& cfg, std::ostream* log = nullptr);

// One log.csv data row, without the trailing newline.
std::string log_row(const GenerationStats& s);

// --- Pruning (SPEC M7) ------------------------------------------------------

struct PruneResult {
    Network net;                  // the pruned network
    double original_fitness = 0;  // fitness of the input network
    double fitness = 0;           // fitness of the pruned network
    int removed = 0;              // reactions removed
};

// Removes each reaction in turn, keeping the removal if the fitness stays
// >= keep_fraction * (fitness of the original network), and repeats full
// passes until no reaction can be removed. Never goes below r_min reactions.
// Deterministic: reactions are tried in order, first to last.
PruneResult prune(const Network& net, const FitnessFn& fitness, int r_min, double keep_fraction = 0.95);

}  // namespace ne
