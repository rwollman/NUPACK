#pragma once
#include "Sequence.h"
#include "Structure.h"
#include "Complex.h"

namespace nupack {

/******************************************************************************************/

struct NamedStrand : Strand, MemberOrdered {
    using base_type = Strand;
    using save_repr = void; // need to hide these from the serializer...
    using load_repr = void;
    Strand complement;
    string name;

    NUPACK_REFLECT_BASE(NamedStrand, Strand, name);

    NamedStrand() = default;
    NamedStrand(Strand s, Strand c, string n) : Strand(std::move(s)), complement(std::move(c)), name(std::move(n)) {}

    NamedStrand reverse_complement(bool wobble=false) const;
    NamedStrand operator~() const {return reverse_complement();}
};

/******************************************************************************************/

struct NamedComplex : Complex {
    using base_type = Complex;
    using save_repr = void; // need to hide these from the serializer...
    using load_repr = void;
    string name;
    vec<string> strand_names;
    vec<Strand> complements;
    real bonus = 0;

    NUPACK_REFLECT_BASE(NamedComplex, Complex, name, strand_names, bonus);

    NamedComplex() = default;

    NamedComplex(Complex c, vec<string> s, vec<Strand> m, string n, real b)
        : Complex(std::move(c)), name(std::move(n)), strand_names(std::move(s)), complements(std::move(m)), bonus(b) {
            NUPACK_REQUIRE(strand_names.size(), ==, Complex::n_strands());
            NUPACK_REQUIRE(complements.size(), ==, Complex::n_strands());
        }

    vec<NamedStrand> strands() const;

    std::size_t symmetry() const;
    vec<std::string_view> lowest_names() const;

    bool operator<(NamedComplex const &) const;
    bool operator==(NamedComplex const &) const;
};

/******************************************************************************************/

// A sequence with an associated name
struct Domain : Sequence, MemberOrdered {
    using base_type = Sequence;
    using save_repr = void; // need to hide these from the serializer...
    using load_repr = void;
    string name;
    Sequence complement;

    Domain() = default;
    Domain(Sequence s, Sequence c, string n)
        : Sequence(std::move(s)), complement(std::move(c)), name(std::move(n)) {}

    Domain reverse_complement(bool wobble=false) const;
    Domain operator~() const {return reverse_complement();}

    NUPACK_REFLECT_BASE(Domain, Sequence, name);
};

using DomainList = small_vec<Domain, 4>;

/******************************************************************************************/

// A strand with a list of component domains
struct TargetStrand : Sequence, MemberOrdered {
    using base_type = Sequence;
    using save_repr = void; // need to hide these from the serializer...
    using load_repr = void;

    DomainList domains;
    string name;

    NUPACK_REFLECT_BASE(TargetStrand, Sequence, domains, name);

    TargetStrand() = default;

    TargetStrand(DomainList d, string name)
        : Sequence(join<Sequence>(d)), domains(std::move(d)), name(std::move(name)) {}
};

/******************************************************************************************/

struct TargetComplex : TotallyOrdered {
    vec<TargetStrand> strands;
    Structure structure;
    string name;
    real bonus = 0;

    TargetComplex() = default;
    TargetComplex(vec<TargetStrand> v, string name, Structure s={}, real bonus=0);

    std::size_t nt() const {return sum(strands, len);}

    explicit operator Complex() const {return Complex(strands);}

    bool operator<(TargetComplex const &) const;
    bool operator==(TargetComplex const &) const;

    NUPACK_REFLECT(TargetComplex, strands, structure, name, bonus);
};

/******************************************************************************************/

}


namespace std {
    template <> struct hash<nupack::NamedComplex> {
        size_t operator()(nupack::NamedComplex const &) const;
    };

    template <> struct hash<nupack::TargetComplex> {
        size_t operator()(nupack::TargetComplex const &) const;
    };

    template <> struct hash<nupack::NamedStrand> : nupack::MemberHash {};

    template <> struct hash<nupack::Domain> : nupack::MemberHash {};

    template <> struct hash<nupack::TargetStrand> : nupack::MemberHash {};
}

/******************************************************************************************/
