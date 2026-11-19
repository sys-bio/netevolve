// Objective functions. v1 has one: the peak-and-trough oscillation fitness
// of SPEC §1.4. Further objectives are meant to be added as separate
// functions alongside it.
#pragma once

#include <span>
#include <vector>

#include "config.hpp"
#include "integrator.hpp"
#include "network.hpp"

namespace ne {

struct OscillationParams {
    double transient_fraction = 0.5;
    int n_target = 5;
    double min_rel_amplitude = 0.01;
    // Absolute floor on the peak-trough difference. Without it, a species
    // that has decayed to ~0 has a ~0 relative threshold, and integrator
    // noise at the atol level (~1e-9) scores as a perfect oscillation.
    // Not in the spec; see PROGRESS.md (M5).
    double min_abs_amplitude = 1e-6;
};

struct Extremum {
    int index;     // sample index within the analysed signal
    double value;
    bool peak;     // true = peak, false = trough
};

// Alternating peaks and troughs of `y`, with every consecutive pair differing
// by more than `threshold`. Smaller ripples are merged away.
std::vector<Extremum> find_extrema(std::span<const double> y, double threshold);

struct SignalScore {
    int n_pairs = 0;       // complete peak-trough pairs
    double damping = 0.0;  // d in [0, 1]
    double score = 0.0;    // f in [0, 1]
};

// Scores one species trajectory. `y` is the full trajectory; the first
// transient_fraction of its samples is discarded here.
SignalScore signal_score(std::span<const double> y, const OscillationParams& p);

struct FitnessParams {
    double t_end = 100;
    int n_points = 1000;
    OscillationParams osc;
    SolverOptions solver;
};

FitnessParams fitness_params(const Config& cfg);

struct FitnessReport {
    double fitness = 0.0;
    SimStatus status = SimStatus::Ok;
    std::vector<SignalScore> species;  // per floating species; empty on failure
};

// Oscillation fitness of a simulation result: 0 unless it succeeded,
// otherwise the best species score.
FitnessReport oscillation_report(const SimResult& sim, const OscillationParams& p);

// Simulates `net` and returns its oscillation fitness with details.
FitnessReport evaluate_oscillation(const Network& net, const FitnessParams& p);

// Just the fitness value, in [0, 1].
double oscillation_fitness(const Network& net, const FitnessParams& p);

}  // namespace ne
