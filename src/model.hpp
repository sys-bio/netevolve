// ODE model built from a Network: ds/dt = N v(s) and its analytic Jacobian.
//
// The state vector holds only the floating species S0 .. S(F-1). Boundary
// species are not state variables: their fixed concentrations are folded
// into constant factors of the rates.
#pragma once

#include <span>
#include <vector>

#include "network.hpp"

namespace ne {

struct Model {
    int n_state = 0;      // F
    int n_reactions = 0;  // R

    // Stoichiometry matrix restricted to floating species, F x R, row-major:
    // N[i * R + j] = (products - reactants) of species i in reaction j.
    std::vector<double> N;

    // v_j(s) = rate_factor[j] * prod over floating_reactants[j] of s_i^coeff.
    // rate_factor[j] = k_j * prod over boundary reactants of X^coeff.
    std::vector<double> rate_factor;
    std::vector<std::vector<Term>> floating_reactants;  // species < F only
};

Model make_model(const Network& net);

// Mass-action rates v(s). `v` must have size R.
void rates(const Model& m, std::span<const double> s, std::span<double> v);

// Right-hand side ds/dt = N v(s). `dsdt` must have size F.
void rhs(const Model& m, std::span<const double> s, std::span<double> dsdt);

// Analytic Jacobian J = N dv/ds, F x F row-major: J[i * F + l] = d(ds_i/dt)/ds_l.
void jacobian(const Model& m, std::span<const double> s, std::span<double> J);

// Convenience wrappers returning new vectors.
std::vector<double> rhs(const Model& m, std::span<const double> s);
std::vector<double> jacobian(const Model& m, std::span<const double> s);

}  // namespace ne
