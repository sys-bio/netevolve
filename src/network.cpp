#include "network.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>

namespace ne {

// ---------------------------------------------------------------------------
// Basics
// ---------------------------------------------------------------------------

std::string species_name(const Network& net, int species) {
    if (species < net.n_floating) return "S" + std::to_string(species);
    return "X" + std::to_string(species - net.n_floating);
}

std::vector<Term> make_terms(std::vector<int> species) {
    std::sort(species.begin(), species.end());
    std::vector<Term> out;
    for (int s : species) {
        if (!out.empty() && out.back().species == s) {
            ++out.back().coeff;
        } else {
            out.push_back({s, 1});
        }
    }
    return out;
}

bool same_stoichiometry(const Reaction& a, const Reaction& b) {
    return a.reactants == b.reactants && a.products == b.products;
}

namespace {

bool involves_floating(const Reaction& r, int n_floating) {
    for (const Term& t : r.reactants)
        if (t.species < n_floating) return true;
    for (const Term& t : r.products)
        if (t.species < n_floating) return true;
    return false;
}

// Empty string if `terms` is canonical with indices in [0, n_species).
std::string check_terms(const std::vector<Term>& terms, int n_species) {
    if (terms.empty()) return "empty side";
    for (std::size_t i = 0; i < terms.size(); ++i) {
        if (terms[i].species < 0 || terms[i].species >= n_species) return "species index out of range";
        if (terms[i].coeff < 1) return "coefficient < 1";
        if (i > 0 && terms[i - 1].species >= terms[i].species) return "terms not canonical";
    }
    return "";
}

}  // namespace

std::string check_validity(const Network& net, int r_min, int r_max) {
    if (net.n_floating < 1) return "network has no floating species";
    if (net.n_boundary < 0) return "negative boundary species count";
    if (static_cast<int>(net.init.size()) != net.n_floating) return "init size != n_floating";
    if (static_cast<int>(net.boundary.size()) != net.n_boundary) return "boundary size != n_boundary";
    for (double c : net.init)
        if (!std::isfinite(c) || c < 0) return "invalid initial concentration";
    for (double c : net.boundary)
        if (!std::isfinite(c) || c < 0) return "invalid boundary concentration";

    const int R = static_cast<int>(net.reactions.size());
    if (R < r_min || R > r_max)
        return "reaction count " + std::to_string(R) + " outside [" + std::to_string(r_min) + ", " +
               std::to_string(r_max) + "]";

    const int n_species = net.n_floating + net.n_boundary;
    for (int i = 0; i < R; ++i) {
        const Reaction& r = net.reactions[i];
        const std::string where = "reaction " + std::to_string(i) + ": ";
        if (auto e = check_terms(r.reactants, n_species); !e.empty()) return where + "reactants: " + e;
        if (auto e = check_terms(r.products, n_species); !e.empty()) return where + "products: " + e;
        if (!std::isfinite(r.k) || r.k <= 0) return where + "rate constant must be finite and > 0";
        if (r.reactants == r.products) return where + "reactants and products are identical";
        if (!involves_floating(r, net.n_floating)) return where + "involves no floating species";
        for (int j = 0; j < i; ++j)
            if (same_stoichiometry(r, net.reactions[j]))
                return where + "duplicates reaction " + std::to_string(j);
    }
    return "";
}

// ---------------------------------------------------------------------------
// Random generation
// ---------------------------------------------------------------------------

Reaction random_reaction(Rng& rng, int n_floating, int n_boundary, double k_min, double k_max) {
    const int type = rng.uniform_int(0, 3);  // 0 UniUni, 1 UniBi, 2 BiUni, 3 BiBi
    const int n_react = (type >= 2) ? 2 : 1;
    const int n_prod = (type == 1 || type == 3) ? 2 : 1;
    const int n_species = n_floating + n_boundary;

    std::vector<int> react, prod;
    for (int i = 0; i < n_react; ++i) react.push_back(rng.uniform_int(0, n_species - 1));
    for (int i = 0; i < n_prod; ++i) prod.push_back(rng.uniform_int(0, n_species - 1));

    Reaction r;
    r.reactants = make_terms(std::move(react));
    r.products = make_terms(std::move(prod));
    r.k = rng.log_uniform(k_min, k_max);
    return r;
}

bool reaction_acceptable(const Network& net, const Reaction& r) {
    if (r.reactants == r.products) return false;
    if (!involves_floating(r, net.n_floating)) return false;
    for (const Reaction& other : net.reactions)
        if (same_stoichiometry(r, other)) return false;
    return true;
}

Network random_network(Rng& rng, const GenParams& p) {
    Network net;
    net.n_floating = p.n_floating;
    net.n_boundary = p.n_boundary;
    net.init.assign(p.n_floating, p.init_conc);
    net.boundary.assign(p.n_boundary, p.boundary_conc);

    constexpr int kMaxAttempts = 1000;
    for (int i = 0; i < p.r_initial; ++i) {
        bool added = false;
        for (int attempt = 0; attempt < kMaxAttempts && !added; ++attempt) {
            Reaction r = random_reaction(rng, p.n_floating, p.n_boundary, p.k_min, p.k_max);
            if (reaction_acceptable(net, r)) {
                net.reactions.push_back(std::move(r));
                added = true;
            }
        }
        if (!added) throw std::runtime_error("random_network: could not find enough distinct valid reactions");
    }
    return net;
}

// ---------------------------------------------------------------------------
// Antimony writer
// ---------------------------------------------------------------------------

namespace {

// Shortest decimal text that parses back to exactly `x`.
std::string format_double(double x) {
    char buf[64];
    auto res = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, res.ptr);
}

std::string write_side(const Network& net, const std::vector<Term>& terms) {
    std::string out;
    for (const Term& t : terms) {
        for (int c = 0; c < t.coeff; ++c) {
            if (!out.empty()) out += " + ";
            out += species_name(net, t.species);
        }
    }
    return out;
}

}  // namespace

std::string write_antimony(const Network& net, const std::string& header) {
    std::ostringstream os;
    if (!header.empty()) {
        std::string h = header;
        std::replace(h.begin(), h.end(), '\n', ' ');
        std::replace(h.begin(), h.end(), '\r', ' ');
        os << "// " << h << "\n";
    }
    if (net.n_boundary > 0) {
        os << "const ";
        for (int j = 0; j < net.n_boundary; ++j) os << (j ? ", " : "") << "X" << j;
        os << ";\n";
    }
    const int R = static_cast<int>(net.reactions.size());
    for (int i = 0; i < R; ++i) {
        const Reaction& r = net.reactions[i];
        os << "J" << i << ": " << write_side(net, r.reactants) << " -> " << write_side(net, r.products)
           << "; k" << i;
        for (const Term& t : r.reactants)
            for (int c = 0; c < t.coeff; ++c) os << "*" << species_name(net, t.species);
        os << ";\n";
    }
    if (R > 0) {
        for (int i = 0; i < R; ++i) os << (i ? " " : "") << "k" << i << " = " << format_double(net.reactions[i].k) << ";";
        os << "\n";
    }
    if (net.n_floating > 0) {
        for (int i = 0; i < net.n_floating; ++i) os << (i ? " " : "") << "S" << i << " = " << format_double(net.init[i]) << ";";
        os << "\n";
    }
    if (net.n_boundary > 0) {
        for (int j = 0; j < net.n_boundary; ++j) os << (j ? " " : "") << "X" << j << " = " << format_double(net.boundary[j]) << ";";
        os << "\n";
    }
    return os.str();
}

// ---------------------------------------------------------------------------
// Antimony parser
// ---------------------------------------------------------------------------

namespace {

enum class Tok { Ident, Number, Colon, Arrow, Semi, Plus, Star, Equals, Comma, End };

struct Token {
    Tok type;
    std::string text;
    int line;
};

const char* tok_desc(Tok t) {
    switch (t) {
        case Tok::Ident: return "identifier";
        case Tok::Number: return "number";
        case Tok::Colon: return "':'";
        case Tok::Arrow: return "'->'";
        case Tok::Semi: return "';'";
        case Tok::Plus: return "'+'";
        case Tok::Star: return "'*'";
        case Tok::Equals: return "'='";
        case Tok::Comma: return "','";
        case Tok::End: return "end of input";
    }
    return "?";
}

bool is_ident_start(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_'; }
bool is_digit(char c) { return c >= '0' && c <= '9'; }

std::vector<Token> lex(const std::string& s) {
    std::vector<Token> toks;
    int line = 1;
    std::size_t i = 0;
    const std::size_t n = s.size();
    while (i < n) {
        const char c = s[i];
        if (c == '\n') { ++line; ++i; continue; }
        if (c == ' ' || c == '\t' || c == '\r') { ++i; continue; }
        if (c == '/' && i + 1 < n && s[i + 1] == '/') {
            while (i < n && s[i] != '\n') ++i;
            continue;
        }
        if (is_ident_start(c)) {
            std::size_t j = i;
            while (j < n && (is_ident_start(s[j]) || is_digit(s[j]))) ++j;
            toks.push_back({Tok::Ident, s.substr(i, j - i), line});
            i = j;
            continue;
        }
        if (is_digit(c) || c == '.') {
            std::size_t j = i;
            while (j < n && is_digit(s[j])) ++j;
            if (j < n && s[j] == '.') { ++j; while (j < n && is_digit(s[j])) ++j; }
            if (j < n && (s[j] == 'e' || s[j] == 'E')) {
                std::size_t k = j + 1;
                if (k < n && (s[k] == '+' || s[k] == '-')) ++k;
                if (k < n && is_digit(s[k])) {
                    while (k < n && is_digit(s[k])) ++k;
                    j = k;
                }
            }
            toks.push_back({Tok::Number, s.substr(i, j - i), line});
            i = j;
            continue;
        }
        if (c == '-' && i + 1 < n && s[i + 1] == '>') { toks.push_back({Tok::Arrow, "->", line}); i += 2; continue; }
        Tok t;
        switch (c) {
            case ':': t = Tok::Colon; break;
            case ';': t = Tok::Semi; break;
            case '+': t = Tok::Plus; break;
            case '*': t = Tok::Star; break;
            case '=': t = Tok::Equals; break;
            case ',': t = Tok::Comma; break;
            default: throw ParseError(line, std::string("unexpected character '") + c + "'");
        }
        toks.push_back({t, std::string(1, c), line});
        ++i;
    }
    toks.push_back({Tok::End, "", line});
    return toks;
}

double parse_number(const Token& t) {
    double v = 0;
    const char* first = t.text.data();
    const char* last = first + t.text.size();
    auto res = std::from_chars(first, last, v);
    if (res.ec != std::errc() || res.ptr != last) throw ParseError(t.line, "invalid number '" + t.text + "'");
    return v;
}

// Species names are S<n> (floating) and X<n> (boundary), n without leading zeros.
struct SpeciesRef {
    bool boundary;
    int index;
};

std::optional<SpeciesRef> parse_species_name(const std::string& name) {
    if (name.size() < 2 || (name[0] != 'S' && name[0] != 'X')) return std::nullopt;
    const std::string digits = name.substr(1);
    if (digits.size() > 5) return std::nullopt;
    for (char c : digits)
        if (!is_digit(c)) return std::nullopt;
    if (digits.size() > 1 && digits[0] == '0') return std::nullopt;
    return SpeciesRef{name[0] == 'X', std::stoi(digits)};
}

struct RawTerm {
    std::string name;
    int coeff;
    int line;
};

struct RawReaction {
    std::string name;
    int line;
    std::vector<RawTerm> lhs, rhs;
    std::vector<Token> rate;  // identifiers of the rate law, in order
};

struct RawAssign {
    std::string name;
    double value;
    int line;
};

class Parser {
public:
    explicit Parser(std::vector<Token> toks) : toks_(std::move(toks)) {}

    void run() {
        while (peek().type != Tok::End) statement();
    }

    std::vector<std::pair<std::string, int>> consts;
    std::vector<RawReaction> reactions;
    std::vector<RawAssign> assigns;

private:
    std::vector<Token> toks_;
    std::size_t pos_ = 0;

    const Token& peek(std::size_t ahead = 0) const { return toks_[std::min(pos_ + ahead, toks_.size() - 1)]; }
    const Token& prev() const { return toks_[pos_ == 0 ? 0 : pos_ - 1]; }
    const Token& next() { return toks_[pos_ < toks_.size() - 1 ? pos_++ : pos_]; }

    [[noreturn]] void unexpected(const char* wanted) const {
        const Token& t = peek();
        std::string got = t.type == Tok::End ? "end of input" : "'" + t.text + "'";
        throw ParseError(t.line, std::string("expected ") + wanted + ", found " + got);
    }

    const Token& expect(Tok type, const char* wanted) {
        if (peek().type != type) unexpected(wanted);
        return next();
    }

    // A missing ';' is reported on the line of the token it should follow.
    void expect_semi(const char* after) {
        if (peek().type != Tok::Semi) throw ParseError(prev().line, std::string("missing ';' after ") + after);
        next();
    }

    void statement() {
        const Token& first = peek();
        if (first.type != Tok::Ident) unexpected("a statement");
        if (first.text == "const") {
            next();
            do {
                const Token& id = expect(Tok::Ident, "a species name");
                consts.emplace_back(id.text, id.line);
            } while (peek().type == Tok::Comma && (next(), true));
            expect_semi("const declaration");
            return;
        }
        if (peek(1).type == Tok::Colon) {
            reaction();
            return;
        }
        if (peek(1).type == Tok::Equals) {
            const Token id = next();
            next();  // '='
            const Token& num = expect(Tok::Number, "a number");
            assigns.push_back({id.text, parse_number(num), id.line});
            expect_semi("assignment");
            return;
        }
        throw ParseError(first.line, "unrecognised statement starting with '" + first.text + "'");
    }

    std::vector<RawTerm> side() {
        std::vector<RawTerm> terms;
        for (;;) {
            int coeff = 1;
            if (peek().type == Tok::Number) {
                const Token& num = next();
                const double v = parse_number(num);
                if (v < 1 || v != std::floor(v) || v > 1000)
                    throw ParseError(num.line, "stoichiometric coefficient must be a positive integer, got '" + num.text + "'");
                coeff = static_cast<int>(v);
            }
            const Token& id = expect(Tok::Ident, "a species name");
            terms.push_back({id.text, coeff, id.line});
            if (peek().type != Tok::Plus) break;
            next();
        }
        return terms;
    }

    void reaction() {
        RawReaction r;
        const Token& name = next();
        r.name = name.text;
        r.line = name.line;
        next();  // ':'
        r.lhs = side();
        expect(Tok::Arrow, "'->'");
        r.rhs = side();
        expect_semi("reaction equation");
        r.rate.push_back(expect(Tok::Ident, "a rate law"));
        while (peek().type == Tok::Star) {
            next();
            r.rate.push_back(expect(Tok::Ident, "a name after '*'"));
        }
        expect_semi("rate law");
        reactions.push_back(std::move(r));
    }
};

}  // namespace

Network parse_antimony(const std::string& text) {
    Parser p(lex(text));
    p.run();

    // First line on which each species name appears, and the largest indices.
    std::map<std::string, int> first_use;
    int max_s = -1, max_x = -1;
    auto note_species = [&](const std::string& name, int line) -> SpeciesRef {
        auto ref = parse_species_name(name);
        if (!ref) throw ParseError(line, "unknown species '" + name + "'");
        auto [it, inserted] = first_use.emplace(name, line);
        if (!inserted) it->second = std::min(it->second, line);
        int& max_index = ref->boundary ? max_x : max_s;
        max_index = std::max(max_index, ref->index);
        return *ref;
    };

    // const declarations
    std::set<std::string> const_names;
    for (const auto& [name, line] : p.consts) {
        auto ref = parse_species_name(name);
        if (!ref || !ref->boundary)
            throw ParseError(line, "only boundary species (X<n>) may be declared const, got '" + name + "'");
        if (!const_names.insert(name).second) throw ParseError(line, "'" + name + "' declared const twice");
        note_species(name, line);
    }

    // Assignments: split into species values and parameter values.
    std::map<std::string, RawAssign> species_vals, param_vals;
    for (const RawAssign& a : p.assigns) {
        auto ref = parse_species_name(a.name);
        auto& dest = ref ? species_vals : param_vals;
        if (dest.count(a.name) || (ref ? param_vals : species_vals).count(a.name))
            throw ParseError(a.line, "'" + a.name + "' assigned twice");
        if (ref) note_species(a.name, a.line);
        dest.emplace(a.name, a);
    }

    // Reactions, first pass: collect species so n_floating / n_boundary are known.
    std::set<std::string> reaction_names;
    for (const RawReaction& r : p.reactions) {
        if (parse_species_name(r.name)) throw ParseError(r.line, "reaction name '" + r.name + "' clashes with a species name");
        if (!reaction_names.insert(r.name).second) throw ParseError(r.line, "duplicate reaction name '" + r.name + "'");
        for (const RawTerm& t : r.lhs) note_species(t.name, t.line);
        for (const RawTerm& t : r.rhs) note_species(t.name, t.line);
    }

    Network net;
    net.n_floating = max_s + 1;
    net.n_boundary = max_x + 1;
    if (net.n_floating < 1) throw ParseError(1, "model has no floating species");

    auto index_of = [&](const SpeciesRef& ref) { return ref.boundary ? net.n_floating + ref.index : ref.index; };

    for (const auto& [name, line] : first_use) {
        auto ref = parse_species_name(name);
        if (ref->boundary && !const_names.count(name))
            throw ParseError(line, "boundary species '" + name + "' must be declared const");
    }

    std::set<std::string> used_params;
    for (const RawReaction& raw : p.reactions) {
        Reaction r;
        auto build = [&](const std::vector<RawTerm>& terms) {
            std::vector<int> idx;
            for (const RawTerm& t : terms)
                for (int c = 0; c < t.coeff; ++c) idx.push_back(index_of(*parse_species_name(t.name)));
            return make_terms(std::move(idx));
        };
        r.reactants = build(raw.lhs);
        r.products = build(raw.rhs);

        // Rate law: exactly one parameter times the reactant species.
        std::string param;
        std::vector<int> rate_species;
        for (const Token& t : raw.rate) {
            if (auto ref = parse_species_name(t.text)) {
                if (!first_use.count(t.text))
                    throw ParseError(t.line, "unknown species '" + t.text + "' in rate law");
                rate_species.push_back(index_of(*ref));
            } else if (param.empty()) {
                param = t.text;
            } else {
                throw ParseError(t.line, "rate law of '" + raw.name + "' has more than one parameter ('" + param +
                                             "' and '" + t.text + "')");
            }
        }
        if (param.empty()) throw ParseError(raw.line, "rate law of '" + raw.name + "' has no rate constant");
        if (make_terms(rate_species) != r.reactants)
            throw ParseError(raw.line, "rate law of '" + raw.name + "' does not match mass action for its reactants");
        if (reaction_names.count(param))
            throw ParseError(raw.line, "rate constant '" + param + "' clashes with a reaction name");
        auto it = param_vals.find(param);
        if (it == param_vals.end()) throw ParseError(raw.line, "rate constant '" + param + "' has no value");
        r.k = it->second.value;
        if (!std::isfinite(r.k) || r.k <= 0)
            throw ParseError(it->second.line, "rate constant '" + param + "' must be > 0");
        used_params.insert(param);

        if (r.reactants == r.products) throw ParseError(raw.line, "reaction '" + raw.name + "' has identical reactants and products");
        if (!involves_floating(r, net.n_floating))
            throw ParseError(raw.line, "reaction '" + raw.name + "' involves no floating species");
        for (const Reaction& other : net.reactions)
            if (same_stoichiometry(r, other))
                throw ParseError(raw.line, "reaction '" + raw.name + "' duplicates an earlier reaction");
        net.reactions.push_back(std::move(r));
    }

    for (const auto& [name, a] : param_vals)
        if (!used_params.count(name)) throw ParseError(a.line, "'" + name + "' is assigned but not used");

    // Concentrations: every species must have a value.
    net.init.assign(net.n_floating, 0.0);
    net.boundary.assign(net.n_boundary, 0.0);
    auto fill = [&](const std::string& name) {
        auto it = species_vals.find(name);
        if (it == species_vals.end()) {
            auto u = first_use.find(name);
            if (u == first_use.end())
                throw ParseError(1, "species '" + name + "' is missing (species must be numbered S0, S1, ... and X0, X1, ... without gaps)");
            throw ParseError(u->second, "species '" + name + "' has no initial value");
        }
        if (!std::isfinite(it->second.value) || it->second.value < 0)
            throw ParseError(it->second.line, "concentration of '" + name + "' must be >= 0");
        return it->second.value;
    };
    for (int i = 0; i < net.n_floating; ++i) net.init[i] = fill("S" + std::to_string(i));
    for (int j = 0; j < net.n_boundary; ++j) net.boundary[j] = fill("X" + std::to_string(j));
    return net;
}

Network load_antimony_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open '" + path + "'");
    std::ostringstream ss;
    ss << in.rdbuf();
    return parse_antimony(ss.str());
}

}  // namespace ne
