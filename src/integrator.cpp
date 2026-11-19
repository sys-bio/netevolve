#include "integrator.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace ne {

const char* to_string(SimStatus s) {
    switch (s) {
        case SimStatus::Ok: return "ok";
        case SimStatus::StepBudget: return "step budget exceeded";
        case SimStatus::WallClock: return "wall-clock limit exceeded";
        case SimStatus::StepTooSmall: return "step size too small";
        case SimStatus::Stiff: return "stiffness detected";
        case SimStatus::MaxConc: return "concentration exceeded max_conc";
        case SimStatus::Negative: return "negative concentration";
        case SimStatus::NonFinite: return "non-finite value";
    }
    return "?";
}

std::string simulation_csv(const SimResult& r, const Network& net) {
    std::string out = "time";
    for (int s = 0; s < r.n_state; ++s) out += "," + species_name(net, s);
    out += "\n";
    char buf[64];
    auto put = [&](double x) {
        auto res = std::to_chars(buf, buf + sizeof buf, x);
        out.append(buf, res.ptr);
    };
    for (int i = 0; i < r.n_samples(); ++i) {
        put(r.t[i]);
        for (int s = 0; s < r.n_state; ++s) {
            out += ",";
            put(r.at(i, s));
        }
        out += "\n";
    }
    return out;
}

// ---------------------------------------------------------------------------
// LU decomposition (partial pivoting)
// ---------------------------------------------------------------------------

bool lu_factor(std::span<double> A, int n, std::span<int> piv) {
    for (int k = 0; k < n; ++k) {
        int p = k;
        double best = std::abs(A[k * n + k]);
        for (int i = k + 1; i < n; ++i) {
            const double v = std::abs(A[i * n + k]);
            if (v > best) { best = v; p = i; }
        }
        piv[k] = p;
        if (best == 0.0 || !std::isfinite(best)) return false;
        if (p != k)
            for (int j = 0; j < n; ++j) std::swap(A[k * n + j], A[p * n + j]);
        const double inv = 1.0 / A[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double f = A[i * n + k] * inv;
            A[i * n + k] = f;
            if (f != 0.0)
                for (int j = k + 1; j < n; ++j) A[i * n + j] -= f * A[k * n + j];
        }
    }
    return true;
}

void lu_solve(std::span<const double> LU, int n, std::span<const int> piv, std::span<double> b) {
    // lu_factor swaps whole rows (L multipliers included), so the stored
    // factors satisfy P A = L U for the product P of all swaps. Apply every
    // swap to b first, then forward-substitute with unit-diagonal L.
    for (int k = 0; k < n; ++k)
        if (piv[k] != k) std::swap(b[k], b[piv[k]]);
    for (int k = 0; k < n; ++k)
        for (int i = k + 1; i < n; ++i) b[i] -= LU[i * n + k] * b[k];
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= LU[i * n + j] * b[j];
        b[i] = s / LU[i * n + i];
    }
}

namespace {

using Clock = std::chrono::steady_clock;
using Vec = std::vector<double>;

constexpr double kFacSafe = 0.9;
constexpr double kFacMin = 0.2;
constexpr double kFacMax = 6.0;

// Shared progress of one simulation: current point, sample bookkeeping and
// guards. Each solver advances it from `t` towards `t_end`.
struct Run {
    const Model& m;
    const SolverOptions& o;
    SimResult& res;
    int n;  // state dimension
    double t_end;
    Vec sample_t;
    int next = 0;  // index of the next sample to fill
    double t = 0.0;
    Vec y;
    Clock::time_point start;
    long attempts = 0;  // steps tried, for the periodic wall-clock check

    Vec f(const Vec& s) {
        Vec out(n);
        rhs(m, s, out);
        return out;
    }

    // Weighted RMS norm of `err` with weights atol + rtol * max(|a|, |b|).
    double err_norm(const Vec& err, const Vec& a, const Vec& b) const {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            const double sc = o.atol + o.rtol * std::max(std::abs(a[i]), std::abs(b[i]));
            const double e = err[i] / sc;
            sum += e * e;
        }
        const double r = std::sqrt(sum / std::max(n, 1));
        return std::isfinite(r) ? r : 1e10;
    }

    // Ok, or the guard that `s` violates.
    SimStatus check(std::span<const double> s) const {
        for (double v : s) {
            if (!std::isfinite(v)) return SimStatus::NonFinite;
            if (v > o.max_conc) return SimStatus::MaxConc;
            if (v < o.min_conc) return SimStatus::Negative;
        }
        return SimStatus::Ok;
    }

    bool out_of_time() {
        if (++attempts % 64 != 0) return false;
        const double secs = std::chrono::duration<double>(Clock::now() - start).count();
        return secs > o.max_wall_seconds;
    }

    // Records sample `next` with value `s`. Returns the guard status.
    SimStatus emit(std::span<const double> s) {
        const SimStatus st = check(s);
        if (st != SimStatus::Ok) return st;
        res.t.push_back(sample_t[next]);
        res.y.insert(res.y.end(), s.begin(), s.end());
        ++next;
        return SimStatus::Ok;
    }

    // Initial step size (Hairer, Norsett & Wanner, Solving ODEs I, II.4).
    double initial_step(const Vec& f0, int order) const {
        Vec scaled(n);
        double d0 = 0, d1 = 0;
        for (int i = 0; i < n; ++i) {
            const double sc = o.atol + o.rtol * std::abs(y[i]);
            d0 += (y[i] / sc) * (y[i] / sc);
            d1 += (f0[i] / sc) * (f0[i] / sc);
        }
        d0 = std::sqrt(d0 / n);
        d1 = std::sqrt(d1 / n);
        double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
        h0 = std::min(h0, t_end - t);
        Vec y1(n);
        for (int i = 0; i < n; ++i) y1[i] = y[i] + h0 * f0[i];
        Vec f1(n);
        rhs(m, y1, f1);
        double d2 = 0;
        for (int i = 0; i < n; ++i) {
            const double sc = o.atol + o.rtol * std::abs(y[i]);
            const double e = (f1[i] - f0[i]) / sc;
            d2 += e * e;
        }
        d2 = std::sqrt(d2 / n) / h0;
        const double dm = std::max(d1, d2);
        const double h1 = (dm <= 1e-15) ? std::max(1e-6, h0 * 1e-3) : std::pow(0.01 / dm, 1.0 / (order + 1));
        double h = std::min(100 * h0, h1);
        if (!std::isfinite(h) || h <= 0) h = 1e-6;
        return std::min(h, t_end - t);
    }

    // Limits h so the last step lands exactly on t_end.
    double clamp_step(double h) const {
        if (t + h >= t_end || t + 1.01 * h > t_end) return t_end - t;
        return h;
    }

    bool step_too_small(double h) const { return h <= 1e-14 * std::max(1.0, std::abs(t)); }
};

// ---------------------------------------------------------------------------
// Dormand-Prince 5(4) (Hairer, Norsett & Wanner, DOPRI5), with the
// 4th-order continuous extension and the DOPRI5 stiffness test.
// ---------------------------------------------------------------------------

SimStatus run_rk45(Run& run, bool stop_on_stiff) {
    constexpr double c2 = 1.0 / 5, c3 = 3.0 / 10, c4 = 4.0 / 5, c5 = 8.0 / 9;
    constexpr double a21 = 1.0 / 5;
    constexpr double a31 = 3.0 / 40, a32 = 9.0 / 40;
    constexpr double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
    constexpr double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
    constexpr double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176,
                     a65 = -5103.0 / 18656;
    constexpr double a71 = 35.0 / 384, a73 = 500.0 / 1113, a74 = 125.0 / 192, a75 = -2187.0 / 6784,
                     a76 = 11.0 / 84;
    constexpr double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200,
                     e6 = 22.0 / 525, e7 = -1.0 / 40;
    constexpr double d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
                     d4 = -10690763975.0 / 1880347072.0, d5 = 701980252875.0 / 199316789632.0,
                     d6 = -1453857185.0 / 822651844.0, d7 = 69997945.0 / 29380423.0;
    (void)c2; (void)c3; (void)c4; (void)c5;  // autonomous system: stage times unused

    const int n = run.n;
    Vec k1 = run.f(run.y), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n);
    Vec ytmp(n), ystage6(n), ynew(n), err(n);
    Vec r1(n), r2(n), r3(n), r4(n), r5(n);  // dense output coefficients
    double h = run.initial_step(k1, 5);
    long steps = 0;
    int stiff_hits = 0, nonstiff_run = 0;
    bool last_rejected = false;

    while (run.t < run.t_end) {
        if (steps >= run.o.max_steps) return SimStatus::StepBudget;
        if (run.out_of_time()) return SimStatus::WallClock;
        if (run.step_too_small(h)) return SimStatus::StepTooSmall;
        h = run.clamp_step(h);
        const auto& y = run.y;

        for (int i = 0; i < n; ++i) ytmp[i] = y[i] + h * a21 * k1[i];
        rhs(run.m, ytmp, k2);
        for (int i = 0; i < n; ++i) ytmp[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        rhs(run.m, ytmp, k3);
        for (int i = 0; i < n; ++i) ytmp[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        rhs(run.m, ytmp, k4);
        for (int i = 0; i < n; ++i) ytmp[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        rhs(run.m, ytmp, k5);
        for (int i = 0; i < n; ++i)
            ystage6[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        rhs(run.m, ystage6, k6);
        for (int i = 0; i < n; ++i)
            ynew[i] = y[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        rhs(run.m, ynew, k7);
        for (int i = 0; i < n; ++i)
            err[i] = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);

        const double en = run.err_norm(err, y, ynew);
        if (en > 1.0) {
            ++run.res.rejected;
            h *= std::max(kFacMin, kFacSafe * std::pow(en, -0.2));
            last_rejected = true;
            continue;
        }

        // Accepted. Stiffness test (DOPRI5): estimate h * |lambda|.
        ++steps;
        ++run.res.steps;
        {
            double num = 0, den = 0;
            for (int i = 0; i < n; ++i) {
                num += (k7[i] - k6[i]) * (k7[i] - k6[i]);
                den += (ynew[i] - ystage6[i]) * (ynew[i] - ystage6[i]);
            }
            if (den > 0 && h * std::sqrt(num / den) > 3.25) {
                nonstiff_run = 0;
                if (++stiff_hits >= 15 && stop_on_stiff) return SimStatus::Stiff;
            } else if (++nonstiff_run >= 6) {
                stiff_hits = 0;
            }
        }

        // Dense output for the samples inside (t, t + h].
        const double t_old = run.t;
        const double t_new = (run.t + h >= run.t_end) ? run.t_end : run.t + h;
        for (int i = 0; i < n; ++i) {
            const double ydiff = ynew[i] - y[i];
            const double bspl = h * k1[i] - ydiff;
            r1[i] = y[i];
            r2[i] = ydiff;
            r3[i] = bspl;
            r4[i] = ydiff - h * k7[i] - bspl;
            r5[i] = h * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * k7[i]);
        }
        while (run.next < static_cast<int>(run.sample_t.size()) && run.sample_t[run.next] <= t_new) {
            const double ts = run.sample_t[run.next];
            SimStatus st;
            if (ts == t_new) {
                st = run.emit(ynew);
            } else {
                const double th = (ts - t_old) / h, th1 = 1.0 - th;
                for (int i = 0; i < n; ++i) ytmp[i] = r1[i] + th * (r2[i] + th1 * (r3[i] + th * (r4[i] + th1 * r5[i])));
                st = run.emit(ytmp);
            }
            if (st != SimStatus::Ok) return st;
        }
        if (const SimStatus st = run.check(ynew); st != SimStatus::Ok) return st;

        run.t = t_new;
        run.y = ynew;
        k1 = k7;  // FSAL
        double fac = kFacSafe * std::pow(std::max(en, 1e-10), -0.2);
        fac = std::min(kFacMax, std::max(kFacMin, fac));
        if (last_rejected) fac = std::min(fac, 1.0);
        h *= fac;
        last_rejected = false;
    }
    return SimStatus::Ok;
}

// ---------------------------------------------------------------------------
// RODAS3 (Sandu et al., Atmos. Environ. 31:3459, 1997), in the formulation
// used by KPP. Order 3 with an embedded order-2 estimate; L-stable and
// stiffly accurate. Dense output by cubic Hermite interpolation.
//
// For each stage s:  (I/(h*gamma) - J) K_s = f(y + sum_j A_sj K_j) + sum_j (C_sj / h) K_j
// y_new = y + sum M_s K_s,   error = sum E_s K_s.
// ---------------------------------------------------------------------------

SimStatus run_rosenbrock(Run& run) {
    constexpr int S = 4;
    constexpr double gamma = 0.5;
    // Lower-triangular A and C, row-wise: stage 2 uses [0], stage 3 [1..2], stage 4 [3..5].
    constexpr double A[6] = {0.0, 2.0, 0.0, 2.0, 0.0, 1.0};
    constexpr double C[6] = {4.0, 1.0, -1.0, 1.0, -1.0, -8.0 / 3.0};
    constexpr double M[S] = {2.0, 0.0, 1.0, 1.0};
    constexpr double E[S] = {0.0, 0.0, 0.0, 1.0};
    constexpr int order = 3;

    const int n = run.n;
    Vec f0 = run.f(run.y), f1(n);
    Vec Jac(static_cast<std::size_t>(n) * n), G(static_cast<std::size_t>(n) * n);
    std::vector<int> piv(n);
    std::vector<Vec> K(S, Vec(n));
    Vec ytmp(n), ftmp(n), ynew(n), err(n);
    double h = run.initial_step(f0, order);
    long steps = 0;
    bool last_rejected = false;

    while (run.t < run.t_end) {
        if (steps >= run.o.max_steps) return SimStatus::StepBudget;
        if (run.out_of_time()) return SimStatus::WallClock;
        if (run.step_too_small(h)) return SimStatus::StepTooSmall;
        h = run.clamp_step(h);
        const auto& y = run.y;

        jacobian(run.m, y, Jac);
        for (std::size_t i = 0; i < G.size(); ++i) G[i] = -Jac[i];
        for (int i = 0; i < n; ++i) G[static_cast<std::size_t>(i) * n + i] += 1.0 / (h * gamma);
        if (!lu_factor(G, n, piv)) {
            ++run.res.rejected;
            h *= 0.5;
            last_rejected = true;
            continue;
        }

        for (int s = 0; s < S; ++s) {
            const int off = (s - 1) * s / 2;  // start of row s in A and C (s >= 1)
            if (s == 0) {
                K[0] = f0;
            } else {
                for (int i = 0; i < n; ++i) {
                    double acc = y[i];
                    for (int j = 0; j < s; ++j) acc += A[off + j] * K[j][i];
                    ytmp[i] = acc;
                }
                rhs(run.m, ytmp, ftmp);
                for (int i = 0; i < n; ++i) {
                    double acc = ftmp[i];
                    for (int j = 0; j < s; ++j) acc += (C[off + j] / h) * K[j][i];
                    K[s][i] = acc;
                }
            }
            lu_solve(G, n, piv, K[s]);
        }
        for (int i = 0; i < n; ++i) {
            double yn = y[i], e = 0.0;
            for (int s = 0; s < S; ++s) {
                yn += M[s] * K[s][i];
                e += E[s] * K[s][i];
            }
            ynew[i] = yn;
            err[i] = e;
        }

        const double en = run.err_norm(err, y, ynew);
        if (en > 1.0) {
            ++run.res.rejected;
            h *= std::max(kFacMin, kFacSafe * std::pow(en, -1.0 / order));
            last_rejected = true;
            continue;
        }

        ++steps;
        ++run.res.steps;
        rhs(run.m, ynew, f1);

        // Cubic Hermite dense output for the samples inside (t, t + h].
        const double t_old = run.t;
        const double t_new = (run.t + h >= run.t_end) ? run.t_end : run.t + h;
        while (run.next < static_cast<int>(run.sample_t.size()) && run.sample_t[run.next] <= t_new) {
            const double ts = run.sample_t[run.next];
            SimStatus st;
            if (ts == t_new) {
                st = run.emit(ynew);
            } else {
                const double th = (ts - t_old) / h, th2 = th * th, th3 = th2 * th;
                const double h00 = 2 * th3 - 3 * th2 + 1, h10 = th3 - 2 * th2 + th;
                const double h01 = -2 * th3 + 3 * th2, h11 = th3 - th2;
                for (int i = 0; i < n; ++i)
                    ytmp[i] = h00 * y[i] + h10 * h * f0[i] + h01 * ynew[i] + h11 * h * f1[i];
                st = run.emit(ytmp);
            }
            if (st != SimStatus::Ok) return st;
        }
        if (const SimStatus st = run.check(ynew); st != SimStatus::Ok) return st;

        run.t = t_new;
        run.y = ynew;
        f0 = f1;
        double fac = kFacSafe * std::pow(std::max(en, 1e-10), -1.0 / order);
        fac = std::min(kFacMax, std::max(kFacMin, fac));
        if (last_rejected) fac = std::min(fac, 1.0);
        h *= fac;
        last_rejected = false;
    }
    return SimStatus::Ok;
}

}  // namespace

SimResult simulate(const Model& m, std::span<const double> y0, double t_end, int n_points, const SolverOptions& opts) {
    if (n_points < 2) throw std::invalid_argument("simulate: n_points must be >= 2");
    if (!(t_end > 0)) throw std::invalid_argument("simulate: t_end must be > 0");
    if (static_cast<int>(y0.size()) != m.n_state) throw std::invalid_argument("simulate: y0 has the wrong size");

    SimResult res;
    res.n_state = m.n_state;
    res.t.reserve(n_points);
    res.y.reserve(static_cast<std::size_t>(n_points) * m.n_state);

    Run run{m, opts, res, m.n_state, t_end, {}, 0, 0.0, Vec(y0.begin(), y0.end()), Clock::now(), 0};
    run.sample_t.resize(n_points);
    for (int i = 0; i < n_points; ++i) run.sample_t[i] = t_end * i / (n_points - 1);
    run.sample_t.back() = t_end;

    if ((res.status = run.emit(run.y)) != SimStatus::Ok) return res;
    if (m.n_state == 0) {
        while (run.next < n_points) run.emit(run.y);
        return res;
    }

    switch (opts.method) {
        case Method::RK45:
            res.status = run_rk45(run, /*stop_on_stiff=*/false);
            break;
        case Method::Rosenbrock:
            res.used_stiff = true;
            res.status = run_rosenbrock(run);
            break;
        case Method::Auto: {
            const SimStatus st = run_rk45(run, /*stop_on_stiff=*/true);
            if (st == SimStatus::Stiff || st == SimStatus::StepBudget || st == SimStatus::StepTooSmall) {
                // Continue from the last accepted point with the stiff solver.
                res.used_stiff = true;
                res.status = run_rosenbrock(run);
            } else {
                res.status = st;
            }
            break;
        }
    }
    return res;
}

}  // namespace ne
