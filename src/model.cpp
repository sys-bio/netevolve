#include "model.hpp"

#include <algorithm>

namespace ne {

namespace {

double ipow(double x, int n) {
    double r = 1.0;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

}  // namespace

Model make_model(const Network& net) {
    Model m;
    const int F = net.n_floating;
    const int R = static_cast<int>(net.reactions.size());
    m.n_state = F;
    m.n_reactions = R;
    m.N.assign(static_cast<std::size_t>(F) * R, 0.0);
    m.rate_factor.resize(R);
    m.floating_reactants.resize(R);

    for (int j = 0; j < R; ++j) {
        const Reaction& r = net.reactions[j];
        double factor = r.k;
        for (const Term& t : r.reactants) {
            if (t.species < F) {
                m.floating_reactants[j].push_back(t);
                m.N[static_cast<std::size_t>(t.species) * R + j] -= t.coeff;
            } else {
                factor *= ipow(net.boundary[t.species - F], t.coeff);
            }
        }
        for (const Term& t : r.products)
            if (t.species < F) m.N[static_cast<std::size_t>(t.species) * R + j] += t.coeff;
        m.rate_factor[j] = factor;
    }
    return m;
}

void rates(const Model& m, std::span<const double> s, std::span<double> v) {
    for (int j = 0; j < m.n_reactions; ++j) {
        double rate = m.rate_factor[j];
        for (const Term& t : m.floating_reactants[j]) rate *= ipow(s[t.species], t.coeff);
        v[j] = rate;
    }
}

void rhs(const Model& m, std::span<const double> s, std::span<double> dsdt) {
    const int F = m.n_state, R = m.n_reactions;
    std::vector<double> v(R);
    rates(m, s, v);
    for (int i = 0; i < F; ++i) {
        double sum = 0.0;
        const double* row = m.N.data() + static_cast<std::size_t>(i) * R;
        for (int j = 0; j < R; ++j) sum += row[j] * v[j];
        dsdt[i] = sum;
    }
}

void jacobian(const Model& m, std::span<const double> s, std::span<double> J) {
    const int F = m.n_state, R = m.n_reactions;
    std::fill(J.begin(), J.end(), 0.0);
    for (int j = 0; j < R; ++j) {
        const auto& terms = m.floating_reactants[j];
        // dv_j/ds_l for each floating reactant l of reaction j.
        for (std::size_t a = 0; a < terms.size(); ++a) {
            const int l = terms[a].species;
            double d = m.rate_factor[j] * terms[a].coeff * ipow(s[l], terms[a].coeff - 1);
            for (std::size_t b = 0; b < terms.size(); ++b)
                if (b != a) d *= ipow(s[terms[b].species], terms[b].coeff);
            // J[:, l] += N[:, j] * dv_j/ds_l
            for (int i = 0; i < F; ++i) {
                const double nij = m.N[static_cast<std::size_t>(i) * R + j];
                if (nij != 0.0) J[static_cast<std::size_t>(i) * F + l] += nij * d;
            }
        }
    }
}

std::vector<double> rhs(const Model& m, std::span<const double> s) {
    std::vector<double> out(m.n_state);
    rhs(m, s, out);
    return out;
}

std::vector<double> jacobian(const Model& m, std::span<const double> s) {
    std::vector<double> out(static_cast<std::size_t>(m.n_state) * m.n_state);
    jacobian(m, s, out);
    return out;
}

}  // namespace ne
