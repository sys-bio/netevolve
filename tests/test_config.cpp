// Configuration parsing (SPEC §2).
#include <string>

#include "config.hpp"
#include "doctest/doctest.h"

using namespace ne;

namespace {

void check_config_error(const std::string& text, const std::string& fragment) {
    try {
        parse_config(text);
        FAIL("expected a ConfigError for: " << text);
    } catch (const ConfigError& e) {
        CHECK_MESSAGE(std::string(e.what()).find(fragment) != std::string::npos, "message: " << e.what());
    }
}

}  // namespace

TEST_CASE("examples/default.cfg matches the built-in defaults") {
    const Config c = load_config_file(NETEVOLVE_EXAMPLES_DIR "/default.cfg");
    CHECK(write_config(c) == write_config(Config{}));
    CHECK_NOTHROW(validate_config(c));
}

TEST_CASE("config parsing: values, comments and blank lines") {
    const Config c = parse_config("# comment\n\nseed = 7\npopulation=50   # inline comment\n  k_max = 1e3 \noutput_dir = runs/a\n");
    CHECK(c.seed == 7);
    CHECK(c.population == 50);
    CHECK(c.k_max == 1000.0);
    CHECK(c.output_dir == "runs/a");
    CHECK(c.generations == Config{}.generations);  // untouched keys keep defaults
}

TEST_CASE("config round trip through write_config") {
    Config c;
    c.seed = 123456789012345ULL;
    c.mut_sigma = 0.123456789;
    c.output_dir = "x";
    CHECK(write_config(parse_config(write_config(c))) == write_config(c));
}

TEST_CASE("config errors are reported clearly") {
    check_config_error("seed = 1\npopulaton = 5\n", "line 2: unknown config key 'populaton'");
    check_config_error("population = many\n", "line 1: invalid value 'many' for 'population'");
    check_config_error("population = 5.5\n", "invalid value");
    check_config_error("seed 5\n", "line 1: expected 'key = value'");
}

TEST_CASE("config validation catches inconsistent values") {
    Config c;
    c.r_initial = 20;  // above r_max
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.k_min = 0;
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.p_rate = c.p_add = c.p_delete = c.p_rewire = 0;
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.recheck_threshold = -1;
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.recheck_rtol = 0;
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.solver = "lsoda";  // not a known solver
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.max_steps = 0;
    CHECK_THROWS_AS(validate_config(c), ConfigError);
    c = Config{};
    c.solver = "rk45";
    CHECK_NOTHROW(validate_config(c));
}
