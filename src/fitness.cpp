#include "fitness.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "model.hpp"

namespace ne {

std::vector<Extremum> find_extrema(std::span<const double> y, double threshold) {
    // 1. Raw candidates: y[i-1] < y[i] >= y[i+1] (peak), y[i-1] > y[i] <= y[i+1] (trough).
    std::vector<Extremum> raw;
    for (std::size_t i = 1; i + 1 < y.size(); ++i) {
        if (y[i - 1] < y[i] && y[i] >= y[i + 1]) raw.push_back({static_cast<int>(i), y[i], true});
        else if (y[i - 1] > y[i] && y[i] <= y[i + 1]) raw.push_back({static_cast<int>(i), y[i], false});
    }

    // 2. Hysteresis: keep a pending extremum; a candidate of the same kind
    //    replaces it if more extreme, a candidate of the opposite kind
    //    commits it only if the two differ by more than `threshold`.
    std::vector<Extremum> out;
    if (raw.empty()) return out;
    Extremum pending = raw[0];
    for (std::size_t i = 1; i < raw.size(); ++i) {
        const Extremum& c = raw[i];
        if (c.peak == pending.peak) {
            if (c.peak ? c.value > pending.value : c.value < pending.value) pending = c;
        } else if (std::abs(c.value - pending.value) > threshold) {
            out.push_back(pending);
            pending = c;
        }
    }
    // 3. Edges. Inside the signal, every committed extremum is confirmed on
    //    both sides by a move larger than `threshold`. At the ends, check the
    //    signal itself: the last pending extremum counts only if the signal
    //    afterwards moves away from it by more than `threshold`, and the first
    //    committed one only if the signal before it does. Otherwise they are
    //    just ripple near the window boundary.
    auto moves_away = [&](const Extremum& e, std::size_t from, std::size_t to) {
        for (std::size_t i = from; i < to; ++i)
            if (e.peak ? e.value - y[i] > threshold : y[i] - e.value > threshold) return true;
        return false;
    };
    if (moves_away(pending, static_cast<std::size_t>(pending.index) + 1, y.size())) out.push_back(pending);
    if (!out.empty() && !moves_away(out.front(), 0, static_cast<std::size_t>(out.front().index)))
        out.erase(out.begin());
    return out;
}

SignalScore signal_score(std::span<const double> y, const OscillationParams& p) {
    SignalScore s;
    const std::size_t start = static_cast<std::size_t>(p.transient_fraction * static_cast<double>(y.size()));
    if (start >= y.size()) return s;
    const std::span<const double> w = y.subspan(start);

    double mean_abs = 0.0;
    for (double v : w) mean_abs += std::abs(v);
    mean_abs /= static_cast<double>(w.size());

    const double threshold = std::max(p.min_rel_amplitude * mean_abs, p.min_abs_amplitude);
    const std::vector<Extremum> ext = find_extrema(w, threshold);
    s.n_pairs = static_cast<int>(ext.size() / 2);
    if (s.n_pairs == 0) return s;

    const double a_first = std::abs(ext[0].value - ext[1].value);
    const std::size_t last = 2 * static_cast<std::size_t>(s.n_pairs) - 2;
    const double a_last = std::abs(ext[last].value - ext[last + 1].value);
    s.damping = std::clamp((a_first - a_last) / a_first, 0.0, 1.0);
    s.score = std::min(s.n_pairs, p.n_target) / static_cast<double>(p.n_target) * (1.0 - s.damping);
    return s;
}

FitnessParams fitness_params(const Config& c) {
    FitnessParams p;
    p.t_end = c.t_end;
    p.n_points = c.n_points;
    p.osc.transient_fraction = c.transient_fraction;
    p.osc.n_target = c.n_target;
    p.osc.min_rel_amplitude = c.min_rel_amplitude;
    p.solver.max_conc = c.max_conc;
    p.solver.max_steps = c.max_steps;
    p.solver.method = c.solver == "rk45"         ? Method::RK45
                      : c.solver == "rosenbrock" ? Method::Rosenbrock
                      : c.solver == "cvode"      ? Method::CVODE
                                                 : Method::Auto;
    p.recheck_solver.method = p.solver.method;
    p.recheck_threshold = c.recheck_threshold;
    p.recheck_solver.rtol = c.recheck_rtol;
    p.recheck_solver.atol = c.recheck_atol;
    p.recheck_solver.max_conc = c.max_conc;
    p.osc.conc_floor = c.osc_conc_floor;
    if (c.objective == "timeseries") {
        p.objective = Objective::TimeSeries;
        p.target = std::make_shared<const TimeSeriesTarget>(load_target_file(c.target_file));
    }
    return p;
}

FitnessReport oscillation_report(const SimResult& sim, const OscillationParams& p) {
    FitnessReport r;
    r.status = sim.status;
    if (!sim.ok()) return r;
    const int F = sim.n_state, n = sim.n_samples();
    const int start = static_cast<int>(p.transient_fraction * static_cast<double>(n));  // as in signal_score
    std::vector<double> traj(n);
    for (int sp = 0; sp < F; ++sp) {
        for (int i = 0; i < n; ++i) traj[i] = sim.at(i, sp);
        const SignalScore s = signal_score(traj, p);
        r.species.push_back(s);
        r.fitness = std::max(r.fitness, s.score);
        // Concentration floor: an oscillating species must stay above it.
        if (p.conc_floor > 0 && s.n_pairs > 0 && start < n) {
            const double lo = *std::min_element(traj.begin() + start, traj.end());
            if (lo < p.conc_floor && r.floor_species < 0) {
                r.floor_species = sp;
                r.floor_min = lo;
            }
        }
    }
    if (r.floor_species >= 0) r.fitness = 0.0;
    return r;
}

FitnessReport evaluate_oscillation(const Network& net, const FitnessParams& p) {
    const Model m = make_model(net);
    const SimResult sim = simulate(m, net.init, p.t_end, p.n_points, p.solver);
    return oscillation_report(sim, p.osc);
}

double oscillation_fitness(const Network& net, const FitnessParams& p) {
    const double f = evaluate_oscillation(net, p).fitness;
    if (f < p.recheck_threshold) return f;
    FitnessParams tight = p;
    tight.solver = p.recheck_solver;
    return std::min(f, evaluate_oscillation(net, tight).fitness);
}

// ---------------------------------------------------------------------------
// Time-series objective
// ---------------------------------------------------------------------------

namespace {

// Splits a data line on whitespace and commas.
std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : line) {
        if (ch == ',' || ch == ' ' || ch == '\t' || ch == '\r') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur += ch;
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

bool to_double(const std::string& s, double& v) {
    const auto res = std::from_chars(s.data(), s.data() + s.size(), v);
    return res.ec == std::errc() && res.ptr == s.data() + s.size() && std::isfinite(v);
}

}  // namespace

TimeSeriesTarget parse_target(const std::string& text) {
    TimeSeriesTarget tg;
    std::istringstream in(text);
    std::string raw;
    int line = 0;
    bool header_allowed = true;
    while (std::getline(in, raw)) {
        ++line;
        const auto hash = raw.find('#');
        const std::vector<std::string> f = fields(hash == std::string::npos ? raw : raw.substr(0, hash));
        if (f.empty()) continue;
        double t = 0, v = 0;
        const bool numeric = f.size() == 2 && to_double(f[0], t) && to_double(f[1], v);
        if (!numeric) {
            if (header_allowed && !(f.size() >= 1 && to_double(f[0], t))) {  // one text header line
                header_allowed = false;
                continue;
            }
            throw std::runtime_error("target line " + std::to_string(line) +
                                     ": expected two numbers, time and concentration");
        }
        header_allowed = false;
        if (t < 0) throw std::runtime_error("target line " + std::to_string(line) + ": time must be >= 0");
        if (!tg.t.empty() && !(t > tg.t.back()))
            throw std::runtime_error("target line " + std::to_string(line) + ": times must be strictly increasing");
        tg.t.push_back(t);
        tg.value.push_back(v);
    }
    if (tg.t.empty()) throw std::runtime_error("target has no data points");
    if (!(tg.t.back() > 0)) throw std::runtime_error("target: the last time must be > 0");
    return tg;
}

TimeSeriesTarget load_target_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open target file '" + path + "'");
    std::ostringstream ss;
    ss << in.rdbuf();
    try {
        return parse_target(ss.str());
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(path + ": " + e.what());
    }
}

TargetReport evaluate_target(const Network& net, const FitnessParams& p) {
    if (!p.target) throw std::invalid_argument("evaluate_target: no target loaded");
    const TimeSeriesTarget& tg = *p.target;
    TargetReport r;
    const SimResult sim = simulate_at(make_model(net), net.init, tg.t, p.solver);
    r.status = sim.status;
    if (!sim.ok()) return r;

    const int n = static_cast<int>(tg.t.size());
    double scale = 0.0;
    for (double v : tg.value) scale += std::abs(v);
    scale /= n;
    if (!(scale > 0)) scale = 1.0;  // all-zero target: plain RMS error

    for (int sp = 0; sp < sim.n_state; ++sp) {
        double sse = 0.0;
        for (int i = 0; i < n; ++i) {
            const double d = sim.at(i, sp) - tg.value[static_cast<std::size_t>(i)];
            sse += d * d;
        }
        const double e = std::sqrt(sse / n) / scale;
        r.nrmse.push_back(e);
        const double score = 1.0 / (1.0 + e);
        if (score > r.fitness) {
            r.fitness = score;
            r.species = sp;
        }
    }
    return r;
}

double network_fitness(const Network& net, const FitnessParams& p) {
    if (p.objective == Objective::Oscillation) return oscillation_fitness(net, p);
    const double f = evaluate_target(net, p).fitness;
    if (f < p.recheck_threshold) return f;
    FitnessParams tight = p;
    tight.solver = p.recheck_solver;
    return std::min(f, evaluate_target(net, tight).fitness);
}

namespace {

double amplitude_threshold(std::span<const double> w, const OscillationParams& p) {
    double mean_abs = 0.0;
    for (double v : w) mean_abs += std::abs(v);
    mean_abs /= static_cast<double>(std::max<std::size_t>(w.size(), 1));
    return std::max(p.min_rel_amplitude * mean_abs, p.min_abs_amplitude);
}

}  // namespace

double window_amplitude(std::span<const double> w, const OscillationParams& p) {
    const std::vector<Extremum> ext = find_extrema(w, amplitude_threshold(w, p));
    if (ext.size() < 2) return 0.0;
    double sum = 0.0;
    for (std::size_t i = 1; i < ext.size(); ++i) sum += std::abs(ext[i].value - ext[i - 1].value);
    return sum / static_cast<double>(ext.size() - 1);
}

LongRunReport long_run_check(const Network& net, const FitnessParams& p, double tolerance, double factor) {
    LongRunReport rep;

    // Species to check: the one with the best score over the normal run.
    const FitnessReport fr = evaluate_oscillation(net, p);
    rep.status = fr.status;
    if (fr.status != SimStatus::Ok || fr.species.empty()) return rep;
    rep.species = 0;
    for (int s = 1; s < static_cast<int>(fr.species.size()); ++s)
        if (fr.species[s].score > fr.species[rep.species].score) rep.species = s;

    // Same sample spacing as the fitness simulation, `factor` times as long.
    // The step budget and wall-clock limit are per simulation and sized for
    // t_end, so they scale with the simulated time too.
    const int intervals = static_cast<int>(std::lround(factor * (p.n_points - 1)));
    const double t_long = p.t_end * intervals / (p.n_points - 1);
    SolverOptions opts = p.solver;
    opts.max_steps = static_cast<long>(std::ceil(factor * static_cast<double>(opts.max_steps)));
    opts.max_wall_seconds *= factor;
    const SimResult sim = simulate(make_model(net), net.init, t_long, intervals + 1, opts);
    rep.status = sim.status;
    if (!sim.ok()) return rep;

    // Windows of (1 - transient_fraction) of the fitness run, in samples.
    const int window = static_cast<int>(std::lround((1.0 - p.osc.transient_fraction) * (p.n_points - 1))) + 1;
    auto window_ending_at = [&](int last_index) {
        std::vector<double> w(window);
        for (int i = 0; i < window; ++i) w[i] = sim.at(last_index - window + 1 + i, rep.species);
        return w;
    };
    rep.amp_ref = window_amplitude(window_ending_at(p.n_points - 1), p.osc);
    rep.amp_last = window_amplitude(window_ending_at(intervals), p.osc);
    if (rep.amp_ref <= 0.0) return rep;
    rep.rel_change = (rep.amp_last - rep.amp_ref) / rep.amp_ref;
    rep.passed = rep.amp_last > 0.0 && std::abs(rep.rel_change) <= tolerance;
    return rep;
}

}  // namespace ne
