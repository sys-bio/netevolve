// Objective functions: the peak-and-trough oscillation fitness of SPEC §1.4,
// and (added after M7) matching a target time series read from a file.
#pragma once

#include <memory>
#include <span>
#include <string>
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
    // If a species that oscillates in the scored window (at least one
    // peak-trough pair) drops below this there, the network scores 0.
    // 0 = off. See Config::osc_conc_floor.
    double conc_floor = 1e-12;
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

// --- Time-series target -------------------------------------------------------
//
// A text file of (time, concentration) pairs, one per line, separated by
// spaces, tabs or a comma. '#' starts a comment; blank lines are skipped; one
// non-numeric header line is allowed before the data. Times must be >= 0,
// strictly increasing, and the last one > 0. Typically alternating peak and
// trough values at the times they should occur.

struct TimeSeriesTarget {
    std::vector<double> t;
    std::vector<double> value;
};

TimeSeriesTarget parse_target(const std::string& text);  // throws std::runtime_error
TimeSeriesTarget load_target_file(const std::string& path);

enum class Objective { Oscillation, TimeSeries };

struct FitnessParams {
    double t_end = 100;
    int n_points = 1000;
    OscillationParams osc;
    SolverOptions solver;
    Objective objective = Objective::Oscillation;
    std::shared_ptr<const TimeSeriesTarget> target;  // set for Objective::TimeSeries
    // Tight-tolerance re-check used by oscillation_fitness (see Config).
    double recheck_threshold = 0.9;
    SolverOptions recheck_solver = tight_solver();

    static SolverOptions tight_solver() {
        SolverOptions s;
        s.rtol = 1e-10;
        s.atol = 1e-14;
        s.max_steps = 5000000;
        s.max_wall_seconds = 30.0;
        return s;
    }
};

// Maps the config to fitness parameters. For objective = timeseries this
// reads target_file (throws if it is missing or malformed).
FitnessParams fitness_params(const Config& cfg);

struct FitnessReport {
    double fitness = 0.0;
    SimStatus status = SimStatus::Ok;
    std::vector<SignalScore> species;  // per floating species; empty on failure
    // Set if the concentration floor zeroed the fitness: the oscillating
    // species that dipped below it, and its minimum in the scored window.
    int floor_species = -1;
    double floor_min = 0.0;
};

// Oscillation fitness of a simulation result: 0 unless it succeeded,
// otherwise the best species score.
FitnessReport oscillation_report(const SimResult& sim, const OscillationParams& p);

// Simulates `net` and returns its oscillation fitness with details.
FitnessReport evaluate_oscillation(const Network& net, const FitnessParams& p);

// The fitness value used for selection and pruning, in [0, 1]. If the normal
// simulation scores >= p.recheck_threshold, the network is simulated again
// with p.recheck_solver (tight tolerances) and the lower score is returned,
// so fitness that is an artifact of integrator error is not selected. A
// failed re-check simulation scores 0.
double oscillation_fitness(const Network& net, const FitnessParams& p);

// --- Time-series objective -----------------------------------------------------
//
// The network is simulated from t = 0 and sampled at the target times. For
// each floating species, the error is the RMS difference from the target
// values divided by the mean |target value| (normalised RMS error, NRMSE),
// and its score is 1 / (1 + NRMSE): 1 for a perfect match, 0.5 at NRMSE = 1.
// The network fitness is the best species score. A failed simulation scores 0.

struct TargetReport {
    double fitness = 0.0;
    SimStatus status = SimStatus::Ok;
    int species = -1;            // best-matching species
    std::vector<double> nrmse;   // per floating species; empty on failure
};

TargetReport evaluate_target(const Network& net, const FitnessParams& p);

// The fitness used for selection and pruning under p.objective, with the
// tight-tolerance re-check (as in oscillation_fitness).
double network_fitness(const Network& net, const FitnessParams& p);

// --- Long-run check (SPEC M7) -----------------------------------------------
//
// Distinguishes sustained oscillation from slow damping. The network is
// re-simulated to factor * t_end with the same sample spacing. For the
// species with the best fitness score, the amplitude (mean peak-to-trough
// difference of the filtered extrema) is measured in two windows of length
// (1 - transient_fraction) * t_end: the fitness window ending at t_end, and
// the window ending at factor * t_end. The check passes if both amplitudes
// are positive and differ by at most `tolerance` relative to the first.

struct LongRunReport {
    bool passed = false;
    SimStatus status = SimStatus::Ok;  // status of the long simulation
    int species = -1;                  // species index checked
    double amp_ref = 0.0;              // amplitude in the window ending at t_end
    double amp_last = 0.0;             // amplitude in the last window
    double rel_change = 0.0;           // (amp_last - amp_ref) / amp_ref
};

LongRunReport long_run_check(const Network& net, const FitnessParams& p, double tolerance = 0.10,
                             double factor = 5.0);

// Mean |peak - trough| over consecutive filtered extrema of `w`, using the
// same threshold as signal_score. 0 if there is no complete pair.
double window_amplitude(std::span<const double> w, const OscillationParams& p);

}  // namespace ne
