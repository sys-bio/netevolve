// M2: stoichiometry matrix, rates, right-hand side and analytic Jacobian.
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest/doctest.h"
#include "model.hpp"
#include "network.hpp"

using namespace ne;

namespace {

Network make_net(int F, int B, std::vector<Reaction> reactions, std::vector<double> init = {},
                 std::vector<double> boundary = {}) {
    Network net;
    net.n_floating = F;
    net.n_boundary = B;
    net.reactions = std::move(reactions);
    net.init = init.empty() ? std::vector<double>(F, 1.0) : init;
    net.boundary = boundary.empty() ? std::vector<double>(B, 1.0) : boundary;
    REQUIRE(check_validity(net, 0, 1000).empty());
    return net;
}

// Central finite-difference Jacobian of the RHS.
std::vector<double> fd_jacobian(const Model& m, std::vector<double> s) {
    const int F = m.n_state;
    std::vector<double> J(static_cast<std::size_t>(F) * F);
    for (int l = 0; l < F; ++l) {
        const double s0 = s[l];
        const double h = 1e-5 * std::max(1.0, std::abs(s0));
        s[l] = s0 + h;
        const std::vector<double> fp = rhs(m, s);
        s[l] = s0 - h;
        const std::vector<double> fm = rhs(m, s);
        s[l] = s0;
        const double step = (s0 + h) - (s0 - h);  // the step actually taken
        for (int i = 0; i < F; ++i) J[static_cast<std::size_t>(i) * F + l] = (fp[i] - fm[i]) / step;
    }
    return J;
}

// max |A - B| / max |A|, the matrix-norm relative error. Entry-wise relative
// error is not meaningful for entries that are zero or cancel to near zero.
double rel_error(const std::vector<double>& A, const std::vector<double>& B) {
    double diff = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < A.size(); ++i) {
        diff = std::max(diff, std::abs(A[i] - B[i]));
        scale = std::max(scale, std::abs(A[i]));
    }
    return scale > 0 ? diff / scale : diff;
}

}  // namespace

TEST_CASE("stoichiometry matrix is built from floating species only") {
    // X0 -> S0;  S0 + S1 -> S2;  S2 -> S1 + S1;  S0 -> X1
    const Network net = make_net(3, 2,
                                 {{make_terms({3}), make_terms({0}), 1.0},
                                  {make_terms({0, 1}), make_terms({2}), 1.0},
                                  {make_terms({2}), make_terms({1, 1}), 1.0},
                                  {make_terms({0}), make_terms({4}), 1.0}});
    const Model m = make_model(net);
    REQUIRE(m.n_state == 3);
    REQUIRE(m.n_reactions == 4);
    const std::vector<double> expected = {
        1, -1, 0, -1,   // S0
        0, -1, 2, 0,    // S1
        0, 1, -1, 0,    // S2
    };
    CHECK(m.N == expected);
}

TEST_CASE("hand-computed RHS for a 3-reaction network") {
    // J0: X0 -> S0           k0 = 2,   v0 = 2 * X0          = 2 * 4         = 8
    // J1: S0 + S1 -> S2      k1 = 0.5, v1 = 0.5 * S0 * S1   = 0.5 * 1.5 * 2 = 1.5
    // J2: S2 -> S1 + S1      k2 = 3,   v2 = 3 * S2          = 3 * 0.7       = 2.1
    // dS0 = v0 - v1 = 6.5;  dS1 = -v1 + 2 v2 = 2.7;  dS2 = v1 - v2 = -0.6
    const Network net = make_net(3, 1,
                                 {{make_terms({3}), make_terms({0}), 2.0},
                                  {make_terms({0, 1}), make_terms({2}), 0.5},
                                  {make_terms({2}), make_terms({1, 1}), 3.0}},
                                 {}, {4.0});
    const Model m = make_model(net);
    const std::vector<double> s = {1.5, 2.0, 0.7};

    std::vector<double> v(3);
    rates(m, s, v);
    CHECK(std::abs(v[0] - 8.0) < 1e-12);
    CHECK(std::abs(v[1] - 1.5) < 1e-12);
    CHECK(std::abs(v[2] - 2.1) < 1e-12);

    const std::vector<double> d = rhs(m, s);
    REQUIRE(d.size() == 3);
    CHECK(std::abs(d[0] - 6.5) < 1e-12);
    CHECK(std::abs(d[1] - 2.7) < 1e-12);
    CHECK(std::abs(d[2] - (-0.6)) < 1e-12);
}

TEST_CASE("Brusselator RHS and Jacobian match the textbook equations") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    const Model m = make_model(net);
    const double A = 1.0, B = 3.0, x = 0.7, y = 2.3;
    const std::vector<double> s = {x, y};
    const std::vector<double> d = rhs(m, s);
    CHECK(std::abs(d[0] - (A + x * x * y - B * x - x)) < 1e-12);
    CHECK(std::abs(d[1] - (B * x - x * x * y)) < 1e-12);

    const std::vector<double> J = jacobian(m, s);
    CHECK(std::abs(J[0] - (2 * x * y - B - 1)) < 1e-12);
    CHECK(std::abs(J[1] - (x * x)) < 1e-12);
    CHECK(std::abs(J[2] - (B - 2 * x * y)) < 1e-12);
    CHECK(std::abs(J[3] - (-x * x)) < 1e-12);
}

TEST_CASE("analytic Jacobian matches central finite differences for 100 random networks") {
    Rng rng(31337);
    double worst = 0.0;
    for (int n = 0; n < 100; ++n) {
        GenParams p;
        p.n_floating = rng.uniform_int(2, 8);
        p.n_boundary = rng.uniform_int(0, 4);
        p.r_initial = rng.uniform_int(3, 12);
        p.boundary_conc = rng.log_uniform(0.1, 10.0);
        const Network net = random_network(rng, p);
        const Model m = make_model(net);

        std::vector<double> s(p.n_floating);
        for (double& x : s) x = rng.log_uniform(0.01, 10.0);

        const double err = rel_error(jacobian(m, s), fd_jacobian(m, s));
        worst = std::max(worst, err);
        REQUIRE_MESSAGE(err < 1e-6, "network " << n << ", relative error " << err << "\n" << write_antimony(net));
    }
    MESSAGE("worst relative Jacobian error: " << worst);
}

TEST_CASE("analytic Jacobian matches finite differences with higher-order stoichiometry") {
    // 2 S0 + S1 -> 3 S0 (cubic rate), 3 S1 -> S0 (cubic), S0 + X0 -> S1 + S1
    const Network net = make_net(2, 1,
                                 {{make_terms({0, 0, 1}), make_terms({0, 0, 0}), 0.7},
                                  {make_terms({1, 1, 1}), make_terms({0}), 1.3},
                                  {make_terms({0, 2}), make_terms({1, 1}), 2.1}},
                                 {}, {1.7});
    const Model m = make_model(net);
    const std::vector<double> s = {0.8, 1.9};
    CHECK(rel_error(jacobian(m, s), fd_jacobian(m, s)) < 1e-6);
}

TEST_CASE("boundary species are not state variables") {
    // J0: X0 -> S0 (k = 2, X0 = 3),  J1: S0 -> X1 (k = 0.5)
    Network net = make_net(1, 2, {{make_terms({1}), make_terms({0}), 2.0}, {make_terms({0}), make_terms({2}), 0.5}}, {},
                           {3.0, 1.0});
    const Model m = make_model(net);

    // The state holds only S0; N has a row only for S0.
    CHECK(m.n_state == 1);
    CHECK(m.N.size() == 2);
    CHECK(m.N == std::vector<double>{1.0, -1.0});

    // X0 enters only as a constant: ds0/dt = 2 * 3 - 0.5 * s0.
    const std::vector<double> s = {4.0};
    const std::vector<double> d = rhs(m, s);
    REQUIRE(d.size() == 1);
    CHECK(std::abs(d[0] - (6.0 - 2.0)) < 1e-12);

    // The Jacobian is F x F, with no column for boundary species.
    const std::vector<double> J = jacobian(m, s);
    REQUIRE(J.size() == 1);
    CHECK(std::abs(J[0] - (-0.5)) < 1e-12);

    // Boundary concentrations are fixed: being consumed (X0) or produced (X1)
    // never changes them. Their derivative is zero by construction: the only
    // place they can change is net.boundary, which the model only reads.
    const std::vector<double> before = net.boundary;
    (void)rhs(m, s);
    (void)jacobian(m, s);
    CHECK(net.boundary == before);

    // A different X0 changes only the constant input term.
    net.boundary[0] = 5.0;
    CHECK(std::abs(rhs(make_model(net), s)[0] - (10.0 - 2.0)) < 1e-12);
}

TEST_CASE("mass conservation: closed A -> B -> A sums to zero") {
    const Network net = make_net(2, 0, {{make_terms({0}), make_terms({1}), 1.7}, {make_terms({1}), make_terms({0}), 0.3}});
    const Model m = make_model(net);
    Rng rng(5);
    for (int i = 0; i < 1000; ++i) {
        const std::vector<double> s = {rng.log_uniform(0.01, 100.0), rng.log_uniform(0.01, 100.0)};
        const std::vector<double> d = rhs(m, s);
        REQUIRE(std::abs(d[0] + d[1]) < 1e-12);
    }
}
