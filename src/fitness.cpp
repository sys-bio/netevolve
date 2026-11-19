#include "fitness.hpp"

#include <algorithm>
#include <cmath>

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
    return p;
}

FitnessReport oscillation_report(const SimResult& sim, const OscillationParams& p) {
    FitnessReport r;
    r.status = sim.status;
    if (!sim.ok()) return r;
    const int F = sim.n_state, n = sim.n_samples();
    std::vector<double> traj(n);
    for (int sp = 0; sp < F; ++sp) {
        for (int i = 0; i < n; ++i) traj[i] = sim.at(i, sp);
        const SignalScore s = signal_score(traj, p);
        r.species.push_back(s);
        r.fitness = std::max(r.fitness, s.score);
    }
    return r;
}

FitnessReport evaluate_oscillation(const Network& net, const FitnessParams& p) {
    const Model m = make_model(net);
    const SimResult sim = simulate(m, net.init, p.t_end, p.n_points, p.solver);
    return oscillation_report(sim, p.osc);
}

double oscillation_fitness(const Network& net, const FitnessParams& p) {
    return evaluate_oscillation(net, p).fitness;
}

}  // namespace ne
