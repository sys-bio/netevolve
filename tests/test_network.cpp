// M1: network representation, RNG, random generation, Antimony I/O.
#include <string>

#include "doctest/doctest.h"
#include "network.hpp"

using namespace ne;

namespace {

// Random generation parameters that vary from network to network, so the
// bulk tests cover more than the default shape.
GenParams varied_params(Rng& rng) {
    GenParams p;
    p.n_floating = rng.uniform_int(2, 8);
    p.n_boundary = rng.uniform_int(0, 4);
    p.r_initial = rng.uniform_int(3, 12);
    return p;
}

// Expects parse_antimony(text) to throw a ParseError on `line` whose message
// contains `fragment`.
void check_parse_error(const std::string& text, int line, const std::string& fragment) {
    try {
        parse_antimony(text);
        FAIL("expected a ParseError for:\n" << text);
    } catch (const ParseError& e) {
        CHECK_MESSAGE(e.line == line, "message: " << e.what());
        CHECK_MESSAGE(std::string(e.what()).find(fragment) != std::string::npos, "message: " << e.what());
        CHECK(std::string(e.what()).rfind("line " + std::to_string(line) + ":", 0) == 0);
    }
}

}  // namespace

TEST_CASE("rng: same seed gives the same stream") {
    Rng a(123), b(123);
    for (int i = 0; i < 1000; ++i) {
        REQUIRE(a.uniform01() == b.uniform01());
        REQUIRE(a.normal() == b.normal());
        REQUIRE(a.uniform_int(-3, 7) == b.uniform_int(-3, 7));
    }
}

TEST_CASE("rng: helpers stay in range") {
    Rng r(7);
    for (int i = 0; i < 100000; ++i) {
        const double u = r.uniform01();
        REQUIRE(u >= 0.0);
        REQUIRE(u < 1.0);
        const int k = r.uniform_int(2, 5);
        REQUIRE(k >= 2);
        REQUIRE(k <= 5);
        const double l = r.log_uniform(0.01, 100.0);
        REQUIRE(l >= 0.01);
        REQUIRE(l <= 100.0);
    }
}

TEST_CASE("random generation is deterministic for a fixed seed") {
    GenParams p;
    Rng a(42), b(42);
    for (int i = 0; i < 100; ++i) {
        const Network na = random_network(a, p);
        const Network nb = random_network(b, p);
        REQUIRE(na == nb);
        REQUIRE(write_antimony(na) == write_antimony(nb));
    }
    // A different seed gives a different network.
    Rng c(43), d(42);
    CHECK_FALSE(random_network(c, p) == random_network(d, p));
}

TEST_CASE("10,000 random networks pass the validity checker") {
    Rng rng(2024);
    for (int i = 0; i < 10000; ++i) {
        const GenParams p = varied_params(rng);
        const Network net = random_network(rng, p);
        const std::string err = check_validity(net, 3, 12);
        REQUIRE_MESSAGE(err.empty(), "network " << i << ": " << err << "\n" << write_antimony(net));
        REQUIRE(static_cast<int>(net.reactions.size()) == p.r_initial);
        for (const Reaction& r : net.reactions) {
            REQUIRE(r.k >= p.k_min);
            REQUIRE(r.k <= p.k_max);
        }
    }
}

TEST_CASE("validity checker rejects invalid networks") {
    Network net;
    net.n_floating = 2;
    net.n_boundary = 1;
    net.init = {1, 1};
    net.boundary = {1};
    net.reactions.push_back({make_terms({2}), make_terms({0}), 1.0});     // X0 -> S0
    net.reactions.push_back({make_terms({0}), make_terms({1}), 1.0});     // S0 -> S1
    net.reactions.push_back({make_terms({1}), make_terms({2}), 1.0});     // S1 -> X0
    REQUIRE(check_validity(net, 3, 12).empty());

    SUBCASE("too few / too many reactions") {
        CHECK_FALSE(check_validity(net, 4, 12).empty());
        CHECK_FALSE(check_validity(net, 1, 2).empty());
    }
    SUBCASE("identical reactants and products") {
        net.reactions.push_back({make_terms({0, 1}), make_terms({1, 0}), 1.0});
        CHECK(check_validity(net, 3, 12).find("identical") != std::string::npos);
    }
    SUBCASE("no floating species") {
        net.n_boundary = 2;
        net.boundary = {1, 1};
        net.reactions.push_back({make_terms({2}), make_terms({3}), 1.0});  // X0 -> X1
        CHECK(check_validity(net, 3, 12).find("no floating") != std::string::npos);
    }
    SUBCASE("duplicate reaction") {
        net.reactions.push_back({make_terms({0}), make_terms({1}), 5.0});
        CHECK(check_validity(net, 3, 12).find("duplicates") != std::string::npos);
    }
    SUBCASE("non-positive rate constant") {
        net.reactions[1].k = 0.0;
        CHECK_FALSE(check_validity(net, 3, 12).empty());
    }
    SUBCASE("species index out of range") {
        net.reactions[1].products = make_terms({5});
        CHECK_FALSE(check_validity(net, 3, 12).empty());
    }
}

TEST_CASE("round trip: write -> parse -> write is identical for 1,000 random networks") {
    Rng rng(99);
    for (int i = 0; i < 1000; ++i) {
        const Network net = random_network(rng, varied_params(rng));
        const std::string text1 = write_antimony(net, "netevolve: generation " + std::to_string(i) + ", fitness 0.5");
        const Network parsed = parse_antimony(text1);
        const std::string text2 = write_antimony(parsed, "netevolve: generation " + std::to_string(i) + ", fitness 0.5");
        REQUIRE_MESSAGE(text1 == text2, "network " << i << "\n" << text1 << "---\n" << text2);
        REQUIRE(parsed == net);  // bit-exact rate constants
    }
}

TEST_CASE("writer output matches the documented format") {
    Network net;
    net.n_floating = 3;
    net.n_boundary = 2;
    net.init = {1, 1, 1};
    net.boundary = {1, 1};
    net.reactions.push_back({make_terms({3}), make_terms({1}), 1.23});        // X0 -> S1
    net.reactions.push_back({make_terms({1, 2}), make_terms({0}), 0.45});     // S1 + S2 -> S0
    net.reactions.push_back({make_terms({0}), make_terms({2, 2}), 3.1});      // S0 -> S2 + S2
    net.reactions.push_back({make_terms({2}), make_terms({4}), 0.8});         // S2 -> X1
    const std::string expected =
        "// netevolve: generation 213, fitness 0.97\n"
        "const X0, X1;\n"
        "J0: X0 -> S1; k0*X0;\n"
        "J1: S1 + S2 -> S0; k1*S1*S2;\n"
        "J2: S0 -> S2 + S2; k2*S0;\n"
        "J3: S2 -> X1; k3*S2;\n"
        "k0 = 1.23; k1 = 0.45; k2 = 3.1; k3 = 0.8;\n"
        "S0 = 1; S1 = 1; S2 = 1;\n"
        "X0 = 1; X1 = 1;\n";
    CHECK(write_antimony(net, "netevolve: generation 213, fitness 0.97") == expected);
}

TEST_CASE("parsing examples/brusselator.ant") {
    const Network net = load_antimony_file(NETEVOLVE_EXAMPLES_DIR "/brusselator.ant");
    REQUIRE(net.n_floating == 2);
    REQUIRE(net.n_boundary == 3);
    REQUIRE(net.reactions.size() == 4);
    // Species indices: S0 = 0, S1 = 1, X0 = 2, X1 = 3, X2 = 4.
    const auto& R = net.reactions;
    CHECK(R[0].reactants == std::vector<Term>{{2, 1}});
    CHECK(R[0].products == std::vector<Term>{{0, 1}});
    CHECK(R[1].reactants == std::vector<Term>{{0, 2}, {1, 1}});
    CHECK(R[1].products == std::vector<Term>{{0, 3}});
    CHECK(R[2].reactants == std::vector<Term>{{0, 1}, {3, 1}});
    CHECK(R[2].products == std::vector<Term>{{1, 1}, {4, 1}});
    CHECK(R[3].reactants == std::vector<Term>{{0, 1}});
    CHECK(R[3].products == std::vector<Term>{{4, 1}});
    for (const Reaction& r : R) CHECK(r.k == 1.0);
    CHECK(net.init == std::vector<double>{1.0, 1.0});
    CHECK(net.boundary == std::vector<double>{1.0, 3.0, 0.0});
    CHECK(check_validity(net, 0, 100).empty());
}

TEST_CASE("parser accepts coefficients, comments and free layout") {
    const Network net = parse_antimony(
        "const X0; // boundary\n"
        "J0: X0 -> 2 S0; k*X0;  J1: S0 + S0 -> S1;\n"
        "  kk*S0*S0;\n"
        "k = 2.5e-1; kk = 3; S0 = 0.5; S1 = 0; X0 = 10;");
    REQUIRE(net.reactions.size() == 2);
    CHECK(net.reactions[0].products == std::vector<Term>{{0, 2}});
    CHECK(net.reactions[0].k == 0.25);
    CHECK(net.reactions[1].reactants == std::vector<Term>{{0, 2}});
    CHECK(net.reactions[1].k == 3.0);
    CHECK(net.init == std::vector<double>{0.5, 0.0});
}

TEST_CASE("malformed input gives an error with a line number") {
    const std::string tail = "k0 = 1; k1 = 1;\nS0 = 1; S1 = 1;\nX0 = 1;\n";

    SUBCASE("unknown species in an equation") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: S0 -> Foo; k1*S0;\n" + tail, 3, "unknown species 'Foo'");
    }
    SUBCASE("unknown species in a rate law") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: S0 -> S1; k1*S7;\n" + tail, 3, "unknown species 'S7'");
    }
    SUBCASE("missing semicolon after a rate law") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0\nJ1: S0 -> S1; k1*S0;\n" + tail, 2, "missing ';'");
    }
    SUBCASE("missing semicolon after an equation") {
        check_parse_error("const X0;\nJ0: X0 -> S0 k0*X0;\nJ1: S0 -> S1; k1*S0;\n" + tail, 2, "missing ';'");
    }
    SUBCASE("missing semicolon at end of input") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: S0 -> S1; k1*S0;\nk0 = 1; k1 = 1;\nS0 = 1; S1 = 1;\nX0 = 1", 6,
                          "missing ';'");
    }
    SUBCASE("rate law does not match the reactants") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: S0 + S1 -> S1; k1*S0;\n" + tail, 3, "does not match");
    }
    SUBCASE("rate law with a repeated species missing") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: 2 S0 -> S1; k1*S0;\n" + tail, 3, "does not match");
    }
    SUBCASE("rate law with two parameters") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: S0 -> S1; k1*k0*S0;\n" + tail, 3, "more than one parameter");
    }
    SUBCASE("boundary species not declared const") {
        check_parse_error("J0: X0 -> S0; k0*X0;\nJ1: S0 -> S1; k1*S0;\n" + tail, 1, "must be declared const");
    }
    SUBCASE("parameter with no value") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: S0 -> S1; k9*S0;\nk0 = 1;\nS0 = 1; S1 = 1;\nX0 = 1;\n", 3,
                          "has no value");
    }
    SUBCASE("unsupported syntax") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0/2;\n", 2, "unexpected character '/'");
    }
    SUBCASE("duplicate reaction") {
        check_parse_error("const X0;\nJ0: X0 -> S0; k0*X0;\nJ1: X0 -> S0; k1*X0;\n" + tail, 3, "duplicates");
    }
}
