// Configuration: a plain "key = value" file with '#' comments (SPEC §2).
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

namespace ne {

struct Config {
    std::uint64_t seed = 42;
    int population = 200;
    int generations = 500;
    double elite_fraction = 0.1;
    double parent_fraction = 0.3;
    int n_floating = 4;
    int n_boundary = 2;
    int r_min = 3;
    int r_max = 12;
    int r_initial = 6;
    double k_min = 0.01;
    double k_max = 100;
    double mut_sigma = 0.5;
    double p_rate = 0.7;
    double p_add = 0.15;
    double p_delete = 0.10;
    double p_rewire = 0.05;
    double t_end = 100;
    int n_points = 1000;
    double transient_fraction = 0.5;
    int n_target = 5;
    double min_rel_amplitude = 0.01;
    double max_conc = 1e6;
    double target_fitness = 0.99;
    // Tight-tolerance re-check: a network scoring >= recheck_threshold is
    // re-simulated at recheck_rtol / recheck_atol and gets the lower of the
    // two scores, so fitness that only exists because of integrator error is
    // not selected. A threshold above 1 turns the re-check off.
    double recheck_threshold = 0.9;
    double recheck_rtol = 1e-10;
    double recheck_atol = 1e-14;
    // ODE solver: "cvode" (SUNDIALS CVODE, BDF), "auto" (RK45, switching to
    // Rosenbrock on stiffness or step budget), "rk45" or "rosenbrock". Used
    // for fitness and the re-check. The default is "cvode"; builds configured
    // with NETEVOLVE_WITH_CVODE=OFF have no CVODE and default to "auto".
#ifdef NETEVOLVE_HAVE_CVODE
    std::string solver = "cvode";
#else
    std::string solver = "auto";
#endif
    long max_steps = 50000;  // step budget per solver, normal-tolerance simulations
    // Oscillation objective: if a species that oscillates in the scored window
    // drops below this concentration there, the network scores 0. Such dips
    // (e.g. to 1e-22) are physically meaningless and numerically unreliable.
    // 0 turns the floor off.
    double osc_conc_floor = 1e-12;
    // Objective: "oscillation" (SPEC §1.4) or "timeseries" (match the
    // time/concentration pairs in target_file; see fitness.hpp).
    std::string objective = "oscillation";
    // Target file for objective = timeseries. A relative path is resolved
    // against the directory of the config file it was read from.
    std::string target_file;
    std::string output_dir = "out";
    int threads = 0;  // 0 = use all cores
};

struct ConfigError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Sets one key from its text value. Throws ConfigError for an unknown key or
// a value that doesn't parse.
void set_config_value(Config& cfg, const std::string& key, const std::string& value);

// Parses config text, starting from the defaults. Errors carry the line number.
Config parse_config(const std::string& text);
Config load_config_file(const std::string& path);

// Checks value ranges and relations between keys. Throws ConfigError.
void validate_config(const Config& cfg);

// Writes every key in "key = value" form; parse_config(write_config(c)) == c.
std::string write_config(const Config& cfg);

}  // namespace ne
