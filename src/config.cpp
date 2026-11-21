#include "config.hpp"

#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace ne {

namespace {

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r");
    return s.substr(b, e - b + 1);
}

template <class T>
T parse_num(const std::string& key, const std::string& value) {
    T v{};
    const char* first = value.data();
    const char* last = first + value.size();
    auto res = std::from_chars(first, last, v);
    if (res.ec != std::errc() || res.ptr != last)
        throw ConfigError("invalid value '" + value + "' for '" + key + "'");
    return v;
}

std::string fmt(double x) {
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, res.ptr);
}

}  // namespace

void set_config_value(Config& c, const std::string& key, const std::string& value) {
    auto I = [&](int& f) { f = parse_num<int>(key, value); };
    auto D = [&](double& f) { f = parse_num<double>(key, value); };
    if (key == "seed") c.seed = parse_num<std::uint64_t>(key, value);
    else if (key == "population") I(c.population);
    else if (key == "generations") I(c.generations);
    else if (key == "elite_fraction") D(c.elite_fraction);
    else if (key == "parent_fraction") D(c.parent_fraction);
    else if (key == "n_floating") I(c.n_floating);
    else if (key == "n_boundary") I(c.n_boundary);
    else if (key == "r_min") I(c.r_min);
    else if (key == "r_max") I(c.r_max);
    else if (key == "r_initial") I(c.r_initial);
    else if (key == "k_min") D(c.k_min);
    else if (key == "k_max") D(c.k_max);
    else if (key == "mut_sigma") D(c.mut_sigma);
    else if (key == "p_rate") D(c.p_rate);
    else if (key == "p_add") D(c.p_add);
    else if (key == "p_delete") D(c.p_delete);
    else if (key == "p_rewire") D(c.p_rewire);
    else if (key == "t_end") D(c.t_end);
    else if (key == "n_points") I(c.n_points);
    else if (key == "transient_fraction") D(c.transient_fraction);
    else if (key == "n_target") I(c.n_target);
    else if (key == "min_rel_amplitude") D(c.min_rel_amplitude);
    else if (key == "max_conc") D(c.max_conc);
    else if (key == "target_fitness") D(c.target_fitness);
    else if (key == "recheck_threshold") D(c.recheck_threshold);
    else if (key == "recheck_rtol") D(c.recheck_rtol);
    else if (key == "recheck_atol") D(c.recheck_atol);
    else if (key == "solver") c.solver = value;
    else if (key == "max_steps") c.max_steps = parse_num<long>(key, value);
    else if (key == "osc_conc_floor") D(c.osc_conc_floor);
    else if (key == "objective") c.objective = value;
    else if (key == "target_file") c.target_file = value;
    else if (key == "output_dir") {
        if (value.empty()) throw ConfigError("output_dir must not be empty");
        c.output_dir = value;
    } else if (key == "threads") I(c.threads);
    else throw ConfigError("unknown config key '" + key + "'");
}

Config parse_config(const std::string& text) {
    Config cfg;
    std::istringstream in(text);
    std::string raw;
    int line = 0;
    while (std::getline(in, raw)) {
        ++line;
        const auto hash = raw.find('#');
        const std::string s = trim(hash == std::string::npos ? raw : raw.substr(0, hash));
        if (s.empty()) continue;
        const auto eq = s.find('=');
        if (eq == std::string::npos) throw ConfigError("line " + std::to_string(line) + ": expected 'key = value'");
        const std::string key = trim(s.substr(0, eq));
        const std::string value = trim(s.substr(eq + 1));
        try {
            set_config_value(cfg, key, value);
        } catch (const ConfigError& e) {
            throw ConfigError("line " + std::to_string(line) + ": " + e.what());
        }
    }
    return cfg;
}

Config load_config_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ConfigError("cannot open config file '" + path + "'");
    std::ostringstream ss;
    ss << in.rdbuf();
    Config c;
    try {
        c = parse_config(ss.str());
    } catch (const ConfigError& e) {
        throw ConfigError(path + ": " + e.what());
    }
    // A relative target_file is relative to this config file's directory.
    if (!c.target_file.empty()) {
        const std::filesystem::path tf(c.target_file);
        if (tf.is_relative())
            c.target_file = (std::filesystem::absolute(path).parent_path() / tf).lexically_normal().string();
    }
    return c;
}

void validate_config(const Config& c) {
    auto require = [](bool ok, const std::string& msg) {
        if (!ok) throw ConfigError("invalid config: " + msg);
    };
    require(c.population >= 2, "population must be >= 2");
    require(c.generations >= 0, "generations must be >= 0");
    require(c.elite_fraction >= 0 && c.elite_fraction < 1, "elite_fraction must be in [0, 1)");
    require(c.parent_fraction > 0 && c.parent_fraction <= 1, "parent_fraction must be in (0, 1]");
    require(c.n_floating >= 1, "n_floating must be >= 1");
    require(c.n_boundary >= 0, "n_boundary must be >= 0");
    require(c.r_min >= 1, "r_min must be >= 1");
    require(c.r_min <= c.r_initial && c.r_initial <= c.r_max, "need r_min <= r_initial <= r_max");
    require(c.k_min > 0 && c.k_min < c.k_max && std::isfinite(c.k_max), "need 0 < k_min < k_max");
    require(c.mut_sigma >= 0, "mut_sigma must be >= 0");
    require(c.p_rate >= 0 && c.p_add >= 0 && c.p_delete >= 0 && c.p_rewire >= 0, "mutation probabilities must be >= 0");
    require(c.p_rate + c.p_add + c.p_delete + c.p_rewire > 0, "mutation probabilities must not all be 0");
    require(c.t_end > 0 && std::isfinite(c.t_end), "t_end must be > 0");
    require(c.n_points >= 3, "n_points must be >= 3");
    require(c.transient_fraction >= 0 && c.transient_fraction < 1, "transient_fraction must be in [0, 1)");
    require(c.n_target >= 1, "n_target must be >= 1");
    require(c.min_rel_amplitude >= 0, "min_rel_amplitude must be >= 0");
    require(c.max_conc > 0, "max_conc must be > 0");
    require(c.recheck_threshold >= 0, "recheck_threshold must be >= 0");
    require(c.recheck_rtol > 0 && c.recheck_atol > 0, "recheck_rtol and recheck_atol must be > 0");
    require(c.solver == "auto" || c.solver == "rk45" || c.solver == "rosenbrock" || c.solver == "cvode",
            "solver must be auto, rk45, rosenbrock or cvode");
#ifndef NETEVOLVE_HAVE_CVODE
    require(c.solver != "cvode", "solver = cvode, but this build has no CVODE (configured with NETEVOLVE_WITH_CVODE=OFF)");
#endif
    require(c.max_steps >= 1, "max_steps must be >= 1");
    require(c.osc_conc_floor >= 0, "osc_conc_floor must be >= 0");
    require(c.objective == "oscillation" || c.objective == "timeseries", "objective must be oscillation or timeseries");
    require(c.objective != "timeseries" || !c.target_file.empty(), "objective = timeseries needs target_file");
    require(c.threads >= 0, "threads must be >= 0");
}

std::string write_config(const Config& c) {
    std::ostringstream os;
    os << "seed = " << c.seed << "\n"
       << "population = " << c.population << "\n"
       << "generations = " << c.generations << "\n"
       << "elite_fraction = " << fmt(c.elite_fraction) << "\n"
       << "parent_fraction = " << fmt(c.parent_fraction) << "\n"
       << "n_floating = " << c.n_floating << "\n"
       << "n_boundary = " << c.n_boundary << "\n"
       << "r_min = " << c.r_min << "\n"
       << "r_max = " << c.r_max << "\n"
       << "r_initial = " << c.r_initial << "\n"
       << "k_min = " << fmt(c.k_min) << "\n"
       << "k_max = " << fmt(c.k_max) << "\n"
       << "mut_sigma = " << fmt(c.mut_sigma) << "\n"
       << "p_rate = " << fmt(c.p_rate) << "\n"
       << "p_add = " << fmt(c.p_add) << "\n"
       << "p_delete = " << fmt(c.p_delete) << "\n"
       << "p_rewire = " << fmt(c.p_rewire) << "\n"
       << "t_end = " << fmt(c.t_end) << "\n"
       << "n_points = " << c.n_points << "\n"
       << "transient_fraction = " << fmt(c.transient_fraction) << "\n"
       << "n_target = " << c.n_target << "\n"
       << "min_rel_amplitude = " << fmt(c.min_rel_amplitude) << "\n"
       << "max_conc = " << fmt(c.max_conc) << "\n"
       << "target_fitness = " << fmt(c.target_fitness) << "\n"
       << "recheck_threshold = " << fmt(c.recheck_threshold) << "\n"
       << "recheck_rtol = " << fmt(c.recheck_rtol) << "\n"
       << "recheck_atol = " << fmt(c.recheck_atol) << "\n"
       << "solver = " << c.solver << "\n"
       << "max_steps = " << c.max_steps << "\n"
       << "osc_conc_floor = " << fmt(c.osc_conc_floor) << "\n"
       << "objective = " << c.objective << "\n"
       << "target_file = " << c.target_file << "\n"
       << "output_dir = " << c.output_dir << "\n"
       << "threads = " << c.threads << "\n";
    return os.str();
}

}  // namespace ne
