// ODE solvers for a Model: adaptive Dormand-Prince RK45 and the L-stable
// Rosenbrock method RODAS3, with step-budget, wall-clock and concentration
// guards.
#pragma once

#include <span>
#include <string>
#include <vector>

#include "model.hpp"

namespace ne {

enum class Method {
    Auto,        // RK45, switching to Rosenbrock on stiffness or step budget
    RK45,        // Dormand-Prince 5(4) only
    Rosenbrock,  // RODAS3 only
};

enum class SimStatus {
    Ok,
    StepBudget,    // exceeded max_steps
    WallClock,     // exceeded max_wall_seconds
    StepTooSmall,  // step size underflow
    Stiff,         // RK45 detected stiffness (only when method == RK45)
    MaxConc,       // a concentration exceeded max_conc
    Negative,      // a concentration fell below min_conc
    NonFinite,     // NaN or infinity in the state
};

const char* to_string(SimStatus s);

struct SolverOptions {
    Method method = Method::Auto;
    double rtol = 1e-6;
    double atol = 1e-9;
    long max_steps = 50000;          // per solver (Auto may use both)
    double max_wall_seconds = 10.0;  // per simulation
    double max_conc = 1e6;
    double min_conc = -1e-8;
};

struct SimResult {
    SimStatus status = SimStatus::Ok;
    int n_state = 0;
    std::vector<double> t;  // sample times; filled up to the point of failure
    std::vector<double> y;  // samples, row-major: y[i * n_state + s]
    long steps = 0;         // accepted steps (both solvers)
    long rejected = 0;      // rejected steps (both solvers)
    bool used_stiff = false;  // Rosenbrock solver was used

    bool ok() const { return status == SimStatus::Ok; }
    int n_samples() const { return static_cast<int>(t.size()); }
    double at(int sample, int species) const { return y[static_cast<std::size_t>(sample) * n_state + species]; }
};

// Integrates from t = 0 to t_end and samples at n_points uniform times
// t_i = t_end * i / (n_points - 1), including both ends. n_points >= 2.
SimResult simulate(const Model& m, std::span<const double> y0, double t_end, int n_points,
                   const SolverOptions& opts = {});

// Time course as CSV: header "time,S0,S1,...", then one row per sample.
// Values are written in the shortest form that reads back exactly.
std::string simulation_csv(const SimResult& r, const Network& net);

// In-house LU factorisation with partial pivoting of an n x n row-major
// matrix, in place. Returns false if the matrix is singular.
bool lu_factor(std::span<double> A, int n, std::span<int> piv);
// Solves A x = b using the factors from lu_factor; b is overwritten by x.
void lu_solve(std::span<const double> LU, int n, std::span<const int> piv, std::span<double> b);

}  // namespace ne
