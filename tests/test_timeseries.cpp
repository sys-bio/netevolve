// Time-series objective: target file parsing, sampling at given times, the
// NRMSE score, config handling and an end-to-end evolve run.
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "config.hpp"
#include "doctest/doctest.h"
#include "evolve.hpp"
#include "fitness.hpp"
#include "integrator.hpp"
#include "model.hpp"
#include "network.hpp"

using namespace ne;
namespace fs = std::filesystem;

namespace {

const Network kDecay = parse_antimony("const X0;\nJ0: S0 -> X0; k0*S0;\nk0 = 0.5;\nS0 = 1;\nX0 = 0;\n");

FitnessParams target_params(const TimeSeriesTarget& tg) {
    FitnessParams p = fitness_params(Config{});
    p.objective = Objective::TimeSeries;
    p.target = std::make_shared<const TimeSeriesTarget>(tg);
    return p;
}

void write(const fs::path& p, const std::string& text) {
    std::ofstream out(p, std::ios::binary);
    out << text;
}

}  // namespace

TEST_CASE("target file: formats, comments and a header line") {
    const TimeSeriesTarget a = parse_target("# comment\ntime concentration\n1 0.5\n2\t3.25\n\n3, 1e-2  # trailing\n");
    REQUIRE(a.t.size() == 3);
    CHECK(a.t == std::vector<double>{1, 2, 3});
    CHECK(a.value == std::vector<double>{0.5, 3.25, 0.01});
    const TimeSeriesTarget b = parse_target("0 1\r\n5 2\r\n");  // CRLF, starts at t = 0
    CHECK(b.t == std::vector<double>{0, 5});
}

TEST_CASE("target file: errors") {
    CHECK_THROWS(parse_target(""));                    // no data
    CHECK_THROWS(parse_target("# only a comment\n"));  // no data
    CHECK_THROWS(parse_target("1 2\n1 3\n"));          // not strictly increasing
    CHECK_THROWS(parse_target("2 2\n1 3\n"));          // decreasing
    CHECK_THROWS(parse_target("-1 2\n1 3\n"));         // negative time
    CHECK_THROWS(parse_target("1 2 3\n"));             // three columns
    CHECK_THROWS(parse_target("1 2\nfoo bar\n"));      // text after data
    CHECK_THROWS(parse_target("head\nhead2\n1 2\n"));  // two header lines
    CHECK_THROWS(parse_target("0 2\n"));               // last time must be > 0
    CHECK_THROWS(parse_target("1 nan\n"));             // non-finite value
    CHECK_THROWS(load_target_file("no/such/target.txt"));
}

TEST_CASE("simulate_at samples exactly at the given times, with and without t = 0") {
    const Model m = make_model(kDecay);
    for (Method method : {Method::Auto, Method::RK45, Method::Rosenbrock}) {
        SolverOptions o;
        o.method = method;
        o.rtol = 1e-10;
        o.atol = 1e-14;
        const std::vector<double> times = {0.5, 1.3, 7.0};
        const SimResult r = simulate_at(m, kDecay.init, times, o);
        REQUIRE(r.ok());
        REQUIRE(r.n_samples() == 3);
        for (int i = 0; i < 3; ++i) {
            CHECK(r.t[i] == times[static_cast<std::size_t>(i)]);
            CHECK(std::abs(r.at(i, 0) - std::exp(-0.5 * times[static_cast<std::size_t>(i)])) < 1e-8);
        }
        const std::vector<double> with0 = {0.0, 2.0};
        const SimResult r0 = simulate_at(m, kDecay.init, with0, o);
        REQUIRE(r0.n_samples() == 2);
        CHECK(r0.at(0, 0) == 1.0);
    }
    const std::vector<double> bad = {1.0, 1.0};
    CHECK_THROWS(simulate_at(m, kDecay.init, bad, {}));
}

TEST_CASE("time-series score: NRMSE and 1 / (1 + NRMSE)") {
    // Target = 2 x the exact decay at t = 1, 2, 3. The simulated S0 is the
    // exact decay, so the error at every point is half the target value.
    TimeSeriesTarget tg;
    for (double t : {1.0, 2.0, 3.0}) {
        tg.t.push_back(t);
        tg.value.push_back(2 * std::exp(-0.5 * t));
    }
    double sse = 0, mean_abs = 0;
    for (double v : tg.value) {
        sse += 0.25 * v * v;
        mean_abs += v / 3;
    }
    const double expected = std::sqrt(sse / 3) / mean_abs;
    const TargetReport r = evaluate_target(kDecay, target_params(tg));
    REQUIRE(r.status == SimStatus::Ok);
    REQUIRE(r.nrmse.size() == 1);
    CHECK(r.nrmse[0] == doctest::Approx(expected).epsilon(1e-5));
    CHECK(r.fitness == doctest::Approx(1.0 / (1.0 + expected)).epsilon(1e-5));
    CHECK(r.species == 0);

    // An exact target scores ~1.
    TimeSeriesTarget exact = tg;
    for (std::size_t i = 0; i < exact.value.size(); ++i) exact.value[i] = std::exp(-0.5 * exact.t[i]);
    CHECK(evaluate_target(kDecay, target_params(exact)).fitness > 0.9999);
}

TEST_CASE("time-series objective: Brusselator peaks and troughs (examples/target_peaks.txt)") {
    const Config cfg = load_config_file(NETEVOLVE_EXAMPLES_DIR "/timeseries.cfg");
    validate_config(cfg);
    const FitnessParams p = fitness_params(cfg);
    REQUIRE(p.objective == Objective::TimeSeries);
    REQUIRE(p.target->t.size() == 11);

    // The network the peaks came from matches them with its S0 ...
    const Network bruss = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const TargetReport r = evaluate_target(bruss, p);
    CHECK(r.species == 0);
    CHECK(r.fitness > 0.999);
    CHECK(network_fitness(bruss, p) > 0.999);  // the re-check runs (score >= 0.9) and agrees

    // ... and a damped version of it doesn't.
    const Network damped = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/damped.ant");
    CHECK(network_fitness(damped, p) < 0.7);
}

TEST_CASE("time-series config: validation, and target_file relative to the config file") {
    Config c;
    c.objective = "timeseries";
    CHECK_THROWS_AS(validate_config(c), ConfigError);  // no target_file
    c.objective = "fit";
    c.target_file = "x.txt";
    CHECK_THROWS_AS(validate_config(c), ConfigError);  // unknown objective

    const fs::path d = fs::temp_directory_path() / "netevolve_test_tsconfig";
    fs::remove_all(d);
    fs::create_directories(d / "sub");
    write(d / "sub" / "peaks.txt", "1 2\n2 1\n");
    write(d / "sub" / "run.cfg", "objective = timeseries\ntarget_file = peaks.txt\n");
    const Config loaded = load_config_file((d / "sub" / "run.cfg").string());
    CHECK(fs::path(loaded.target_file).is_absolute());
    CHECK(fs::equivalent(loaded.target_file, d / "sub" / "peaks.txt"));
    CHECK(fitness_params(loaded).target->t.size() == 2);
    // write_config keeps the resolved path, so a run.cfg in an output directory still works.
    CHECK(parse_config(write_config(loaded)).target_file == loaded.target_file);
    fs::remove_all(d);
}

TEST_CASE("time-series objective: evolve and --numtrials run end to end") {
    const fs::path d = fs::temp_directory_path() / "netevolve_test_tsevolve";
    fs::remove_all(d);
    Config cfg = load_config_file(NETEVOLVE_EXAMPLES_DIR "/timeseries.cfg");
    cfg.population = 24;
    cfg.generations = 6;
    cfg.seed = 5;
    cfg.output_dir = (d / "single").string();
    const EvolutionResult res = run_evolve_command(cfg);
    CHECK(res.best_fitness > 0.0);
    CHECK(res.best_fitness <= 1.0);
    CHECK(fs::exists(d / "single" / "best_final.ant"));
    // The best network's recorded fitness is its time-series fitness.
    CHECK(network_fitness(res.best, fitness_params(cfg)) == res.best_fitness);

    cfg.output_dir = (d / "trials").string();
    const auto trials = run_trials_command(cfg, 2);
    REQUIRE(trials.size() == 2);
    {
        std::ifstream in(d / "trials" / "summary.csv");
        std::string header, row;
        std::getline(in, header);
        std::getline(in, row);
        CHECK(row.substr(row.rfind(',') + 1) == "n/a");  // no long-run check for this objective
    }
    fs::remove_all(d);
}
