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
