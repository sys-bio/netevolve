#include "evolve.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <ostream>
#include <stdexcept>

#include "fitness.hpp"
#include "integrator.hpp"
#include "model.hpp"

namespace ne {

// ---------------------------------------------------------------------------
// Mutation
// ---------------------------------------------------------------------------

namespace {

bool mutate_rate(Network& net, Rng& rng, const Config& cfg) {
    if (net.reactions.empty()) return false;
    Reaction& r = net.reactions[rng.uniform_int(0, static_cast<int>(net.reactions.size()) - 1)];
    r.k = std::clamp(r.k * std::exp(cfg.mut_sigma * rng.normal()), cfg.k_min, cfg.k_max);
    return true;
}

bool mutate_add(Network& net, Rng& rng, const Config& cfg) {
    if (static_cast<int>(net.reactions.size()) >= cfg.r_max) return false;
    Reaction r = random_reaction(rng, net.n_floating, net.n_boundary, cfg.k_min, cfg.k_max);
    if (!reaction_acceptable(net, r)) return false;
    net.reactions.push_back(std::move(r));
    return true;
}

bool mutate_delete(Network& net, Rng& rng, const Config& cfg) {
    if (static_cast<int>(net.reactions.size()) <= cfg.r_min) return false;
    net.reactions.erase(net.reactions.begin() + rng.uniform_int(0, static_cast<int>(net.reactions.size()) - 1));
    return true;
}

// Expands a canonical side into one entry per molecule.
std::vector<int> expand(const std::vector<Term>& terms) {
    std::vector<int> out;
    for (const Term& t : terms)
        for (int c = 0; c < t.coeff; ++c) out.push_back(t.species);
    return out;
}

bool mutate_rewire(Network& net, Rng& rng, const Config&) {
    const int n_species = net.n_floating + net.n_boundary;
    if (net.reactions.empty() || n_species < 2) return false;
    const int j = rng.uniform_int(0, static_cast<int>(net.reactions.size()) - 1);

    std::vector<int> react = expand(net.reactions[j].reactants);
    std::vector<int> prod = expand(net.reactions[j].products);
    const int slot = rng.uniform_int(0, static_cast<int>(react.size() + prod.size()) - 1);
    int& target = slot < static_cast<int>(react.size()) ? react[slot] : prod[slot - react.size()];
    // A different species, uniformly.
    int replacement = rng.uniform_int(0, n_species - 2);
    if (replacement >= target) ++replacement;
    target = replacement;

    Reaction r;
    r.reactants = make_terms(std::move(react));
    r.products = make_terms(std::move(prod));
    r.k = net.reactions[j].k;

    // Acceptable against every other reaction.
    Network others = net;
    others.reactions.erase(others.reactions.begin() + j);
    if (!reaction_acceptable(others, r)) return false;
    net.reactions[j] = std::move(r);
    return true;
}

}  // namespace

bool apply_mutation(Network& net, MutationKind kind, Rng& rng, const Config& cfg) {
    switch (kind) {
        case MutationKind::Rate: return mutate_rate(net, rng, cfg);
        case MutationKind::Add: return mutate_add(net, rng, cfg);
        case MutationKind::Delete: return mutate_delete(net, rng, cfg);
        case MutationKind::Rewire: return mutate_rewire(net, rng, cfg);
    }
    return false;
}

MutationKind choose_mutation(Rng& rng, const Config& cfg) {
    const double total = cfg.p_rate + cfg.p_add + cfg.p_delete + cfg.p_rewire;
    double u = rng.uniform01() * total;
    if ((u -= cfg.p_rate) < 0) return MutationKind::Rate;
    if ((u -= cfg.p_add) < 0) return MutationKind::Add;
    if ((u -= cfg.p_delete) < 0) return MutationKind::Delete;
    if (cfg.p_rewire > 0) return MutationKind::Rewire;
    // Rounding pushed u past the end: return the last kind with p > 0.
    if (cfg.p_delete > 0) return MutationKind::Delete;
    if (cfg.p_add > 0) return MutationKind::Add;
    return MutationKind::Rate;
}

Network make_child(const Network& parent, Rng& rng, const Config& cfg) {
    for (int attempt = 0; attempt < kMutationAttempts; ++attempt) {
        Network child = parent;
        if (apply_mutation(child, choose_mutation(rng, cfg), rng, cfg) &&
            check_validity(child, cfg.r_min, cfg.r_max).empty())
            return child;
    }
    return parent;
}

GenParams gen_params(const Config& cfg) {
    GenParams p;
    p.n_floating = cfg.n_floating;
    p.n_boundary = cfg.n_boundary;
    p.r_initial = cfg.r_initial;
    p.k_min = cfg.k_min;
    p.k_max = cfg.k_max;
    return p;
}

// ---------------------------------------------------------------------------
// Generation loop
// ---------------------------------------------------------------------------

namespace {

struct Individual {
    Network net;
    double fitness = 0.0;
    bool evaluated = false;
};

GenerationStats stats_of(int gen, const std::vector<Individual>& sorted) {
    GenerationStats s;
    s.generation = gen;
    const std::size_t n = sorted.size();
    s.best = sorted.front().fitness;
    s.best_reactions = static_cast<int>(sorted.front().net.reactions.size());
    double sum = 0.0;
    for (const Individual& ind : sorted) sum += ind.fitness;
    s.mean = sum / static_cast<double>(n);
    // `sorted` is in descending order of fitness.
    s.median = (n % 2) ? sorted[n / 2].fitness : 0.5 * (sorted[n / 2 - 1].fitness + sorted[n / 2].fitness);
    return s;
}

}  // namespace

EvolutionResult evolve(const Config& cfg, const FitnessFn& fitness, const EvolutionCallbacks& cb) {
    validate_config(cfg);
    Rng rng(cfg.seed);
    const int P = cfg.population;
    int n_elite = static_cast<int>(std::lround(cfg.elite_fraction * P));
    if (cfg.elite_fraction > 0) n_elite = std::max(n_elite, 1);
    n_elite = std::min(n_elite, P - 1);
    const int n_parents = std::clamp(static_cast<int>(std::lround(cfg.parent_fraction * P)), 1, P);

    std::vector<Individual> pop(P);
    const GenParams gp = gen_params(cfg);
    for (Individual& ind : pop) ind.net = random_network(rng, gp);

    EvolutionResult result;
    result.best_fitness = -1.0;
    for (int gen = 0;; ++gen) {
        // Evaluate (elites carried over keep their fitness: evaluation is pure).
        for (Individual& ind : pop) {
            if (!ind.evaluated) {
                ind.fitness = fitness(ind.net);
                ind.evaluated = true;
            }
        }
        // Sort by fitness, best first; ties keep their order.
        std::stable_sort(pop.begin(), pop.end(),
                         [](const Individual& a, const Individual& b) { return a.fitness > b.fitness; });

        const GenerationStats s = stats_of(gen, pop);
        result.log.push_back(s);
        if (cb.on_generation) cb.on_generation(s);
        if (pop.front().fitness > result.best_fitness) {
            result.best = pop.front().net;
            result.best_fitness = pop.front().fitness;
            result.best_generation = gen;
            if (cb.on_improvement) cb.on_improvement(gen, result.best, result.best_fitness);
        }
        if (result.best_fitness >= cfg.target_fitness || gen >= cfg.generations) break;

        // Next generation: elites unchanged, the rest mutated children of
        // parents drawn uniformly from the top n_parents.
        std::vector<Individual> next;
        next.reserve(P);
        for (int i = 0; i < n_elite; ++i) next.push_back(pop[i]);
        while (static_cast<int>(next.size()) < P) {
            const Individual& parent = pop[rng.uniform_int(0, n_parents - 1)];
            next.push_back({make_child(parent.net, rng, cfg), 0.0, false});
        }
        pop = std::move(next);
    }
    return result;
}

// ---------------------------------------------------------------------------
// The `evolve` command
// ---------------------------------------------------------------------------

namespace {

std::string fmt(double x) {
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, res.ptr);
}

void write_text(const std::filesystem::path& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write '" + path.string() + "'");
    out << text;
    if (!out) throw std::runtime_error("error writing '" + path.string() + "'");
}

std::string ant_header(int gen, double fitness) {
    return "netevolve: generation " + std::to_string(gen) + ", fitness " + fmt(fitness);
}

}  // namespace

std::string log_row(const GenerationStats& s) {
    return std::to_string(s.generation) + "," + fmt(s.best) + "," + fmt(s.mean) + "," + fmt(s.median) + "," +
           std::to_string(s.best_reactions);
}

EvolutionResult run_evolve_command(const Config& cfg, std::ostream* log) {
    validate_config(cfg);
    namespace fs = std::filesystem;
    const fs::path dir(cfg.output_dir);
    fs::create_directories(dir);
    write_text(dir / "run.cfg", write_config(cfg));

    std::ofstream log_csv(dir / "log.csv", std::ios::binary);
    if (!log_csv) throw std::runtime_error("cannot write '" + (dir / "log.csv").string() + "'");
    log_csv << "generation,best,mean,median,best_reactions\n";

    const FitnessParams fp = fitness_params(cfg);
    EvolutionCallbacks cb;
    cb.on_generation = [&](const GenerationStats& s) {
        log_csv << log_row(s) << "\n";
        log_csv.flush();
        if (log) *log << "gen " << s.generation << "  best " << s.best << "  mean " << s.mean << "  reactions "
                      << s.best_reactions << "\n";
    };
    cb.on_improvement = [&](int gen, const Network& best, double f) {
        char name[32];
        std::snprintf(name, sizeof name, "best_gen%04d.ant", gen);
        write_text(dir / name, write_antimony(best, ant_header(gen, f)));
    };

    EvolutionResult res = evolve(cfg, [&](const Network& n) { return oscillation_fitness(n, fp); }, cb);
    log_csv.close();

    write_text(dir / "best_final.ant", write_antimony(res.best, ant_header(res.best_generation, res.best_fitness)));
    const SimResult sim = simulate(make_model(res.best), res.best.init, fp.t_end, fp.n_points, fp.solver);
    write_text(dir / "best_final.csv", simulation_csv(sim, res.best));
    return res;
}

}  // namespace ne
