// Seeded random number source. Every stochastic component of netevolve draws
// from one of these, so a given seed reproduces a run exactly.
//
// The helpers below are written by hand rather than using <random>
// distributions, whose output is implementation-defined. Only the
// mt19937_64 engine itself (fully specified by the standard) is used.
#pragma once

#include <cmath>
#include <cstdint>
#include <random>

namespace ne {

struct Rng {
    std::mt19937_64 eng;

    explicit Rng(std::uint64_t seed) : eng(seed) {}

    // Uniform double in [0, 1), 53 bits of precision.
    double uniform01() { return static_cast<double>(eng() >> 11) * 0x1.0p-53; }

    // Uniform double in [lo, hi).
    double uniform(double lo, double hi) { return lo + (hi - lo) * uniform01(); }

    // Uniform integer in [lo, hi] (inclusive). Unbiased (rejection sampling).
    int uniform_int(int lo, int hi) {
        const std::uint64_t span = static_cast<std::uint64_t>(hi - lo) + 1;
        const std::uint64_t limit = UINT64_MAX - (UINT64_MAX % span);
        std::uint64_t x;
        do { x = eng(); } while (x >= limit);
        return lo + static_cast<int>(x % span);
    }

    // Standard normal deviate (Box-Muller; consumes two engine draws).
    double normal() {
        double u1;
        do { u1 = uniform01(); } while (u1 <= 0.0);
        const double u2 = uniform01();
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }

    // Log-uniform double in [lo, hi), lo > 0.
    double log_uniform(double lo, double hi) {
        return std::exp(uniform(std::log(lo), std::log(hi)));
    }
};

}  // namespace ne
