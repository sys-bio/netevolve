// netevolve command-line entry point.
#include <charconv>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

#include "config.hpp"
#include "evolve.hpp"
#include "fitness.hpp"
#include "integrator.hpp"
#include "model.hpp"
#include "network.hpp"

namespace {

void print_usage(std::ostream& os) {
    os << "usage:\n"
          "  netevolve evolve   --config FILE [--seed N] [--out DIR]\n"
          "  netevolve simulate --model FILE.ant [--t-end T] [--points N] [--csv OUT.csv]\n"
          "  netevolve fitness  --model FILE.ant [--config FILE]\n"
          "  netevolve prune    --model FILE.ant [--config FILE] [--out FILE.ant]\n"
          "  netevolve --version\n";
}

// Command-line misuse: reported with the usage text, exit code 2.
struct UsageError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Parses "--flag value" pairs after the command. Every flag takes a value.
std::map<std::string, std::string> parse_flags(int argc, char** argv, const std::set<std::string>& allowed) {
    std::map<std::string, std::string> flags;
    for (int i = 2; i < argc; ++i) {
        const std::string flag = argv[i];
        if (!allowed.count(flag)) throw UsageError("unknown option '" + flag + "' for '" + argv[1] + "'");
        if (i + 1 >= argc) throw UsageError("option '" + flag + "' needs a value");
        if (!flags.emplace(flag, argv[++i]).second) throw UsageError("option '" + flag + "' given twice");
    }
    return flags;
}

std::string require_flag(const std::map<std::string, std::string>& f, const std::string& name) {
    auto it = f.find(name);
    if (it == f.end()) throw UsageError("missing required option '" + name + "'");
    return it->second;
}

template <class T>
T flag_number(const std::string& flag, const std::string& text) {
    T v{};
    auto res = std::from_chars(text.data(), text.data() + text.size(), v);
    if (res.ec != std::errc() || res.ptr != text.data() + text.size())
        throw UsageError("invalid value '" + text + "' for '" + flag + "'");
    return v;
}

void write_file(const std::string& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write '" + path + "'");
    out << content;
    if (!out) throw std::runtime_error("error writing '" + path + "'");
}

void print_long_run(std::ostream& os, const ne::Network& net, const ne::FitnessParams& fp) {
    const ne::LongRunReport lr = ne::long_run_check(net, fp);
    os << "long-run check (to " << 5 * fp.t_end << "): ";
    if (lr.status != ne::SimStatus::Ok) {
        os << "FAILED, simulation: " << ne::to_string(lr.status) << "\n";
        return;
    }
    if (lr.species < 0 || lr.amp_ref <= 0) {
        os << "FAILED, no oscillation at t_end\n";
        return;
    }
    os << (lr.passed ? "passed" : "FAILED") << ", " << ne::species_name(net, lr.species) << " amplitude "
       << lr.amp_ref << " -> " << lr.amp_last << " (" << 100 * lr.rel_change << "%)\n";
}

int cmd_simulate(int argc, char** argv) {
    const auto f = parse_flags(argc, argv, {"--model", "--t-end", "--points", "--csv"});
    const ne::Network net = ne::load_antimony_file(require_flag(f, "--model"));
    const ne::Config defaults;
    const double t_end = f.count("--t-end") ? flag_number<double>("--t-end", f.at("--t-end")) : defaults.t_end;
    const int points = f.count("--points") ? flag_number<int>("--points", f.at("--points")) : defaults.n_points;
    if (!(t_end > 0)) throw UsageError("--t-end must be > 0");
    if (points < 2) throw UsageError("--points must be >= 2");

    ne::SolverOptions opts;
    opts.max_conc = defaults.max_conc;
    const ne::SimResult r = ne::simulate(ne::make_model(net), net.init, t_end, points, opts);
    const std::string csv = ne::simulation_csv(r, net);
    if (f.count("--csv")) write_file(f.at("--csv"), csv);
    else std::cout << csv;

    if (!r.ok()) {
        std::cerr << "warning: simulation stopped at t = " << (r.t.empty() ? 0.0 : r.t.back()) << ": "
                  << ne::to_string(r.status) << "\n";
        return 1;
    }
    return 0;
}

int cmd_fitness(int argc, char** argv) {
    const auto f = parse_flags(argc, argv, {"--model", "--config"});
    const ne::Network net = ne::load_antimony_file(require_flag(f, "--model"));
    ne::Config cfg;
    if (f.count("--config")) cfg = ne::load_config_file(f.at("--config"));
    ne::validate_config(cfg);

    const ne::FitnessReport rep = ne::evaluate_oscillation(net, ne::fitness_params(cfg));
    std::cout << "status: " << ne::to_string(rep.status) << "\n";
    for (std::size_t s = 0; s < rep.species.size(); ++s) {
        const ne::SignalScore& sc = rep.species[s];
        std::cout << ne::species_name(net, static_cast<int>(s)) << ": pairs " << sc.n_pairs << ", damping "
                  << sc.damping << ", score " << sc.score << "\n";
    }
    std::cout << "fitness " << rep.fitness << "\n";
    return 0;
}

int cmd_evolve(int argc, char** argv) {
    const auto f = parse_flags(argc, argv, {"--config", "--seed", "--out"});
    ne::Config cfg = ne::load_config_file(require_flag(f, "--config"));
    if (f.count("--seed")) cfg.seed = flag_number<std::uint64_t>("--seed", f.at("--seed"));
    if (f.count("--out")) {
        if (f.at("--out").empty()) throw UsageError("--out must not be empty");
        cfg.output_dir = f.at("--out");
    }
    ne::validate_config(cfg);
    std::cout << "netevolve evolve: seed " << cfg.seed << ", population " << cfg.population << ", generations "
              << cfg.generations << ", threads " << ne::resolve_threads(cfg.threads) << ", output "
              << cfg.output_dir << "\n";
    const auto start = std::chrono::steady_clock::now();
    const ne::EvolutionResult res = ne::run_evolve_command(cfg, &std::cout);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "best fitness " << res.best_fitness << " (generation " << res.best_generation << ", "
              << res.best.reactions.size() << " reactions)\n";
    std::cout << res.n_evaluations << " evaluations in " << secs << " s ("
              << (secs > 0 ? res.n_evaluations / secs : 0.0) << " per second)\n";
    print_long_run(std::cout, res.best, ne::fitness_params(cfg));
    return 0;
}

int cmd_prune(int argc, char** argv) {
    const auto f = parse_flags(argc, argv, {"--model", "--config", "--out"});
    const ne::Network net = ne::load_antimony_file(require_flag(f, "--model"));
    ne::Config cfg;
    if (f.count("--config")) cfg = ne::load_config_file(f.at("--config"));
    ne::validate_config(cfg);
    const ne::FitnessParams fp = ne::fitness_params(cfg);

    const ne::PruneResult pr =
        ne::prune(net, [&](const ne::Network& n) { return ne::oscillation_fitness(n, fp); }, cfg.r_min);
    char header[128];
    std::snprintf(header, sizeof header, "netevolve: pruned, fitness %.6g (original %.6g)", pr.fitness,
                  pr.original_fitness);
    const std::string text = ne::write_antimony(pr.net, header);

    // The model goes to stdout unless --out is given; the summary then goes to stderr.
    const bool to_file = f.count("--out") > 0;
    if (to_file) write_file(f.at("--out"), text);
    else std::cout << text;
    std::ostream& info = to_file ? std::cout : std::cerr;
    info << "pruned " << pr.removed << " of " << net.reactions.size() << " reactions: fitness "
         << pr.original_fitness << " -> " << pr.fitness << "\n";
    print_long_run(info, pr.net, fp);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(std::cerr);
        return 2;
    }
    const std::string_view cmd = argv[1];
    try {
        if (cmd == "--version") {
            std::cout << "netevolve " << NETEVOLVE_VERSION << "\n";
            return 0;
        }
        if (cmd == "--help" || cmd == "-h") {
            print_usage(std::cout);
            return 0;
        }
        if (cmd == "simulate") return cmd_simulate(argc, argv);
        if (cmd == "fitness") return cmd_fitness(argc, argv);
        if (cmd == "evolve") return cmd_evolve(argc, argv);
        if (cmd == "prune") return cmd_prune(argc, argv);
        throw UsageError("unknown command '" + std::string(cmd) + "'");
    } catch (const UsageError& e) {
        std::cerr << "error: " << e.what() << "\n";
        print_usage(std::cerr);
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
