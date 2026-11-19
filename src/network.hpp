// Network representation, validity checking, random generation and
// Antimony I/O.
#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "rng.hpp"

namespace ne {

// Species are addressed by a single index:
//   0 .. F-1        floating species S0 .. S(F-1)
//   F .. F+B-1      boundary species X0 .. X(B-1)
struct Term {
    int species = 0;
    int coeff = 1;  // stoichiometric coefficient, >= 1
    bool operator==(const Term&) const = default;
};

// Irreversible mass-action reaction. `reactants` and `products` are kept in
// canonical form: sorted by species index, one Term per species.
struct Reaction {
    std::vector<Term> reactants;
    std::vector<Term> products;
    double k = 1.0;
    bool operator==(const Reaction&) const = default;
};

struct Network {
    int n_floating = 0;
    int n_boundary = 0;
    std::vector<Reaction> reactions;
    std::vector<double> init;      // initial concentrations, size n_floating
    std::vector<double> boundary;  // fixed concentrations, size n_boundary
    bool operator==(const Network&) const = default;
};

// Parameters for random generation.
struct GenParams {
    int n_floating = 4;
    int n_boundary = 2;
    int r_initial = 6;
    double k_min = 0.01;
    double k_max = 100.0;
    double init_conc = 1.0;
    double boundary_conc = 1.0;
};

inline bool is_boundary(const Network& net, int species) { return species >= net.n_floating; }
std::string species_name(const Network& net, int species);

// Builds a canonical Term list from a list of species indices (repeats allowed).
std::vector<Term> make_terms(std::vector<int> species);

// True if both reactions have the same reactant and product multisets.
bool same_stoichiometry(const Reaction& a, const Reaction& b);

// Returns an empty string if `net` is valid, otherwise a description of the
// first problem found. Checks structure (sizes, indices, canonical terms,
// k > 0, concentrations >= 0) and the SPEC §1.1 constraints.
std::string check_validity(const Network& net, int r_min, int r_max);

// One random reaction of type UniUni, UniBi, BiUni or BiBi, with participants
// drawn uniformly from all species and k log-uniform in [k_min, k_max).
// The result may be invalid (e.g. no floating species); callers check.
Reaction random_reaction(Rng& rng, int n_floating, int n_boundary, double k_min, double k_max);

// True if `r` alone satisfies the per-reaction validity constraints and does
// not duplicate any reaction already in `net`.
bool reaction_acceptable(const Network& net, const Reaction& r);

// Random valid network with p.r_initial reactions. Throws std::runtime_error
// if it cannot find enough distinct valid reactions.
Network random_network(Rng& rng, const GenParams& p);

// --- Antimony I/O (SPEC §3) -------------------------------------------------

// Writes `net` in the restricted Antimony subset. If `header` is non-empty it
// is written as a first line comment: "// <header>".
std::string write_antimony(const Network& net, const std::string& header = "");

struct ParseError : std::runtime_error {
    int line;
    ParseError(int line_, const std::string& msg)
        : std::runtime_error("line " + std::to_string(line_) + ": " + msg), line(line_) {}
};

// Parses the restricted Antimony subset. Throws ParseError on any error.
Network parse_antimony(const std::string& text);

// Reads a file and parses it. Throws std::runtime_error if it can't be read.
Network load_antimony_file(const std::string& path);

}  // namespace ne
