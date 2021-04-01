#include "Bind.h"
#include <nupack/types/Sequence.h>
#include <nupack/types/PairList.h>
#include <nupack/types/Domain.h>

namespace nupack {

/******************************************************************************************/

std::optional<PairList> request(Type<PairList>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto z = r.request<std::string_view>()) return PairList(*z);
    return msg.error("Cannot convert to PairList");
}

/******************************************************************************************/

void render(Document &doc, Type<PairList> t) {
    doc.type(t, "core.PairList");
    render_public(doc, t);
    static_assert(has_lt<PairList>);
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
    doc.method(t, "new", rebind::construct<std::string_view>(t));
    doc.method(t, "new", rebind::construct<typename PairList::data_type>(t));
    doc.method(t, "^", [](PairList const &v, PairList const &w) {
        NUPACK_REQUIRE(len(v), ==, len(w));
        return v ^ w;
    });
    doc.method(t, "dp", &PairList::dp<pair_data_type const &>);
    doc.method(t, "pseudoknots", &PairList::pseudoknots);
}

/******************************************************************************************/

// 0 = dna weak, 1 = dna strong, 2 = rna weak, 3 = rna strong
std::atomic<unsigned> print_as_rna{0};
static constexpr const std::array<char, 16> rna_names{
    {'A', 'C', 'G', 'U', 'R', 'M', 'S', 'W', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '_'}
};

int set_sequence_type_strong(int rna) noexcept {
    return print_as_rna.exchange(rna);
}

void set_sequence_type_weak(bool rna) noexcept {
    auto last = print_as_rna.exchange(rna ? 2 : 0);
    if (last | 1) print_as_rna.store(last); // undo if it was strong
}

/******************************************************************************************/

void render(Document &doc, Type<Base> t) {
    doc.type(t, "core.Base");
    doc.method(t, "new", rebind::construct<char>(t));
    doc.method(t, "new", rebind::construct<Base>(t));
    doc.method(t, "letter", [](Base b) {
        NUPACK_REQUIRE(b, <, 16, "Invalid base", b);
        return (print_as_rna.load() < 2) ? rna_names[b] : Base::names[b];
    });

    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

std::optional<Base> request(Type<Base>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto p = r.request<char>()) return Base(*p);
    if (auto p = r.request<std::string_view>()) if (p->size() == 1) return Base((*p)[0]);
    if (auto p = r.request<std::string>()) if (p->size() == 1) return Base((*p)[0]);
    return msg.error("not convertible to Base");
}

/******************************************************************************************/

std::string sequence_string_impl(Sequence const &v) {
    std::string s;
    s.reserve(len(s));
    for (auto b : v) NUPACK_REQUIRE(b, <, 16, "Invalid base", b);

    if (print_as_rna.load() < 2) {
        for (auto b : v) s.push_back(rna_names[b]);
    } else {
        for (auto b : v) s.push_back(Base::names[b]);
    }
    return s;
}

void render(Document &doc, Type<Sequence> t) {
    doc.type(t, "core.Sequence");
    render_comparisons(doc, t);
    render_json(doc, t);
    doc.method(t, "new", rebind::construct<std::string_view>(t));
    doc.method(t, "new", rebind::construct<Sequence>(t));
    doc.method(t, "{}", sequence_string_impl);
    doc.method(t, "__hash__", [](Sequence const &s) {return std::hash<Sequence>()(s);});
    doc.method(t, "__len__", [](Sequence const &s) {return s.size();});
    doc.method(t, "^", [](Sequence const &x, Sequence const &y) {
        NUPACK_REQUIRE(len(x), ==, len(y));
        return hamming_distance(x, y);
    });
    doc.method(t, "nt", [](Sequence const &s) {return s.size();});
    doc.method(t, "__getitem__", [](Sequence const &s, std::size_t i) {return s.at(i);});
    // doc.method(t, "__setitem__", [](Sequence &s, std::size_t i, Base b) {s.at(i) = b;});
    doc.method(t, "__contains__", [] (Sequence const &s, Base b) {return contains(s, b);});
    doc.method(t, "reverse_complement", &reverse_complement);
    doc.method(t, "is_determined", [](Sequence const &s) {return all_of(s, is_determined);});
    doc.method(t, "has_wildcard", [](Sequence const &s) {return any_of(s, has_wildcard);});
    doc.method(t, "is_canonical", [](Sequence const &s) {return all_of(s, is_canonical);});

    doc.function("core.to_sequences", [](SequenceList const &v) {return v;});
    doc.function("core.to_sequences", [](vec<std::string_view> const &v) {return to_sequences(v);});
    doc.function("core.to_sequences", [](std::string_view s) {return to_sequences(s);});
    doc.function("core.to_sequences", [](Sequence s) {return SequenceList{std::move(s)};});

    doc.function("core.set_sequence_type", set_sequence_type_strong);

    doc.function("core.seq_distance", [](SequenceList const &a, SequenceList const &b) -> std::size_t {
        NUPACK_REQUIRE(len(a), ==, len(b), "Complexes are differently sized");
        return std::inner_product(a.begin(), a.end(), b.begin(), std::size_t(0), std::plus<>(), [](Sequence const &a, Sequence const &b) {
                NUPACK_REQUIRE(len(a), ==, len(b), "Sequence length does not match");
                return std::inner_product(a.begin(), a.end(), b.begin(), std::size_t(0), std::plus<>(), [](Base a, Base b) {
                    auto const &A = Base::masks.at(a);
                    auto const &B = Base::masks.at(b);
                    return !((A[0] && B[0]) || (A[1] && B[1]) || (A[2] && B[2]) || (A[3] && B[3]));
                });
            });
    });
}


std::optional<Sequence> request(Type<Sequence>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto p = r.request<std::string_view>()) return Sequence(*p);
    if (auto p = r.request<std::string>()) return Sequence(*p);
    return msg.error("not convertible to Sequence");
}

static_assert(std::is_same_v<rebind::response_method<Sequence>, rebind::ADL>);

/******************************************************************************************/

void render(Document &doc, Type<Strand> t) {
    doc.type(t, "core.RawStrand");
    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
    doc.method(t, "new", rebind::construct<Sequence>(t));
    doc.method(t, "new", rebind::construct<Strand>(t));
}

std::optional<Strand> request(Type<Strand>, rebind::Variable const &r, rebind::Dispatch &msg) {
    if (auto p = r.request<std::string_view>()) return Strand(*p);
    if (auto p = r.request<std::string>()) return Strand(*p);
    return msg.error("not convertible to Strand");
}

/******************************************************************************************/

void render(Document &doc, Type<Complex> t) {
    doc.type(t, "core.RawComplex");
    doc.method(t, "new", rebind::construct<StrandList>(t));
    doc.method(t, "new", rebind::construct<Complex>(t));
    doc.method(t, "new", [](std::string s) {return Complex(to_sequences<StrandList>(s));});

    render_json(doc, t);
    render_hash(doc, t);

    doc.method(t, "strands", [](Complex const &x) {return x.strands();});
    doc.method(t, "__getitem__", [](Complex const &x, std::size_t i) {return x.strands().at(i);});
    doc.method(t, "__contains__", [](Complex const &x, Strand const &s) {return contains(x.strands(), s);});
    doc.method(t, "__xor__", [](Complex const &a, Complex const &b) {
        NUPACK_REQUIRE(a.positions, ==, b.positions, "Complexes are differently sized");
        return hamming_distance(a.catenated, b.catenated);
    });
    doc.method(t, "__len__", [](Complex const &x) {return x.n_strands();});
    doc.method(t, "nt", [](Complex const &x) {return x.size();});
    doc.method(t, "symmetry", [](Complex const &x) {return rotational_symmetry(x.strands());});
    doc.method(t, "lowest_rotation", [](Complex x) {x.rotate_lowest(); return x;});

    doc.method(t, "__hash__", [](Complex const &a) {return hash_of(lowest_rotation(a.strands()));});
    doc.method(t, "==", [](Complex const &a, Complex const &b) {return lowest_rotation(a.views()) == lowest_rotation(b.views());});
    doc.method(t, "!=", [](Complex const &a, Complex const &b) {return lowest_rotation(a.views()) != lowest_rotation(b.views());});
    doc.method(t, "<", [](Complex const &a, Complex const &b)  {return lowest_rotation(a.views()) <  lowest_rotation(b.views());});
    doc.method(t, ">", [](Complex const &a, Complex const &b)  {return lowest_rotation(a.views()) >  lowest_rotation(b.views());});
    doc.method(t, "<=", [](Complex const &a, Complex const &b) {return lowest_rotation(a.views()) >= lowest_rotation(b.views());});
    doc.method(t, ">=", [](Complex const &a, Complex const &b) {return lowest_rotation(a.views()) <= lowest_rotation(b.views());});
}

std::optional<Complex> request(Type<Complex> t, rebind::Variable const &v, rebind::Dispatch &msg) {
    if (auto p = v.request<StrandList>()) return Complex(*p);
    return msg.error("Cannot convert to Complex", t);
}

/******************************************************************************************/

bool complex_eq(NamedComplex &t, NamedComplex &u) {
    bool eq = t == u;
    if (eq) {
        if (t.name.empty() && !u.name.empty()) t.name = u.name;
        if (u.name.empty() && !t.name.empty()) u.name = t.name;

        if (t.bonus == 0 && u.bonus != 0) t.bonus = 0;
        if (u.bonus == 0 && t.bonus != 0) u.bonus = 0;
    }
    return eq;
}

auto get_strands(NamedComplex const &c) {return c.strands();}
auto get_strands(TargetComplex const &c) {return c.strands;}

template <class C>
std::string complex_name(C const &c) {
    if (!c.name.empty()) return c.name;
    string o;
    o.push_back('(');
    bool first = true;
    for (auto const &s : get_strands(c)) {
        if (!std::exchange(first, false)) o.push_back('+');
        o += s.name;
    }
    o.push_back(')');
    return o;
}

/******************************************************************************************/

void render(Document &doc, Type<NamedComplex> t) {
    doc.type(t, "named.Complex");
    doc.method(t, "new", rebind::construct<Complex, vec<string>, vec<Strand>, string, real>(t));
    doc.method(t, "new", rebind::construct<NamedComplex>(t));

    doc.method(t, "name", complex_name<NamedComplex>);
    doc.method(t, "strands", &NamedComplex::strands);
    doc.method(t, ".bonus", &NamedComplex::bonus);
    doc.method(t, "symmetry", &NamedComplex::symmetry);

    // This is a nasty business below but is inescapable if dealing with colliding optional names...
    doc.method(t, "==", [](NamedComplex &t, NamedComplex &u) {return complex_eq(t, u);});
    doc.method(t, "!=", [](NamedComplex &t, NamedComplex &u) {return !complex_eq(t, u);});

    doc.method(t, "<", std::less<NamedComplex>());
    doc.method(t, ">", std::greater<NamedComplex>());

    doc.method(t, "<=", std::less_equal<NamedComplex>());
    doc.method(t, ">=", std::greater_equal<NamedComplex>());

    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

void render(Document &doc, Type<NamedStrand> t) {
    doc.type(t, "named.Strand");
    doc.method(t, "new", rebind::construct<Strand, Strand, string>(t));
    doc.method(t, "new", rebind::construct<NamedStrand>(t));
    doc.method(t, "~", &NamedStrand::operator~);
    doc.method(t, "reverse_complement", &NamedStrand::reverse_complement);

    doc.method(t, ".name", &NamedStrand::name);

    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

void render(Document &doc, Type<Domain> t) {
    doc.type(t, "core.Domain");
    doc.method(t, "new", rebind::construct<Sequence, Sequence, string>(t));
    doc.method(t, "new", rebind::construct<Domain>(t));
    doc.method(t, "~", &Domain::operator~);
    doc.method(t, "reverse_complement", &Domain::reverse_complement);

    doc.method(t, ".name", &Domain::name);

    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

void render(Document &doc, Type<TargetStrand> t) {
    doc.type(t, "core.TargetStrand");
    doc.method(t, "new", rebind::construct<DomainList, string>(t));
    doc.method(t, "new", rebind::construct<TargetStrand>(t));
    doc.method(t, "~", [](TargetStrand d) {
        std::reverse(d.domains.begin(), d.domains.end());
        for (auto &d : d.domains) d = ~std::move(d);
        std::reverse(static_cast<Sequence&>(d).begin(), static_cast<Sequence&>(d).end());
        return d;
    });

    doc.method(t, ".domains", &TargetStrand::domains);
    doc.method(t, ".name", &TargetStrand::name);

    render_comparisons(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

// This is a nasty business below but is inescapable if dealing with colliding optional names and complexes D:
bool target_complex_eq(TargetComplex &t, TargetComplex &u) {
    bool eq = t == u;
    if (eq) {
        if (t.name.empty() && !u.name.empty()) t.name = u.name;
        if (u.name.empty() && !t.name.empty()) u.name = t.name;

        if (t.structure.empty() && !u.structure.empty()) t.structure = u.structure;
        if (u.structure.empty() && !t.structure.empty()) u.structure = t.structure;

        if (t.bonus == 0 && u.bonus != 0) t.bonus = 0;
        if (u.bonus == 0 && t.bonus != 0) u.bonus = 0;
    }
    return eq;
}

void render(Document &doc, Type<TargetComplex> t) {
    doc.type(t, "core.TargetComplex");
    doc.method(t, "new", rebind::construct<vec<TargetStrand>, string, Structure, real>(t));
    doc.method(t, "new", rebind::construct<TargetComplex>(t));
    doc.method(t, "nt", &TargetComplex::nt);
    doc.method(t, "__len__", [](TargetComplex const &d) {return len(d.strands);});
    doc.method(t, "__getitem__", [](TargetComplex const &d, std::size_t i) {return d.strands.at(i);});

    doc.method(t, "==", target_complex_eq);

    doc.method(t, "!=", [](TargetComplex &t, TargetComplex &u) {return !target_complex_eq(t, u);});

    doc.method(t, "<", std::less<TargetComplex>());
    doc.method(t, ">", std::greater<TargetComplex>());

    doc.method(t, "<=", std::less_equal<TargetComplex>());
    doc.method(t, ">=", std::greater_equal<TargetComplex>());

    doc.method(t, "name", complex_name<TargetComplex>);

    NUPACK_PUBLIC(TargetComplex, strands, structure, bonus);
    render_json(doc, t);
    render_hash(doc, t);
}


/******************************************************************************************/

void render(Document &doc, Type<Structure> t) {
    doc.type(t, "core.Structure");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "new", rebind::construct<string>(t));
    doc.method(t, "new", rebind::construct<Structure>(t));
    doc.method(t, "new", rebind::construct<PairList, Nicks>(t));
    doc.method(t, ".values", &PairList::values);

    doc.method(t, "dp",     &Structure::dp);
    doc.method(t, "dp_rle", &Structure::dp_rle);

    doc.method(t, "dotparensplus",     &Structure::dp);
    doc.method(t, "rle_dotparensplus", &Structure::dp_rle);
    doc.method(t, "nicks",             &Structure::nicks);

    render_comparisons(doc, t);
    // render_public(doc, t);
    render_json(doc, t);
    render_hash(doc, t);
}

/******************************************************************************************/

}
