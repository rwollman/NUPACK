#pragma once
#include "../standard/Optional.h"
#include "../reflect/Reflection.h"
#include "../standard/Vec.h"
#include "../iteration/Range.h"

namespace nupack { namespace newdesign {

struct Weight {
    Optional<string> tube;
    Optional<string> complex;
    Optional<string> strand;
    Optional<string> domain;
    real weight;

    Weight() = default;
    Weight(Optional<string>, Optional<string>, Optional<string>, Optional<string>, real);

    NUPACK_REFLECT(Weight, tube, complex, strand, domain, weight);
};

struct Design; // forward declaration to avoid Design.h include

/* for mapping complexes back to strand and domain names */
/* NOTE: assumes that all */
struct ReversedComplex {
    /* range of nucleotide indices -> name of domain */
    std::map<std::pair<uint, uint>, string> _domains;
    /* range of nucleotide indices -> name of strand */
    std::map<std::pair<uint, uint>, string> _strands;
    /* length of the complex */
    uint _N;

    ReversedComplex() = default;
    ReversedComplex(Design const &, uint);

    void reverse_map(Design const &, uint);

    vec<string> domains() const;
    vec<string> strands() const;

    NUPACK_REFLECT(ReversedComplex, _domains, _strands, _N);
};


using ComplexWeights = std::map<uint, vec<real>>;

struct Weights {
    vec<Weight> specifications;

    ComplexWeights per_complex;
    std::map<uint, ComplexWeights> per_tube;
    std::map<uint, ReversedComplex> reversed_complexes;
    /* objectives are only weighted at the top-level */
    /* only the multitubeobjective uses the rest of the weights */
    vec<real> objective_weights;

    void add(Weight w) {specifications.emplace_back(std::move(w));}
    void add_objective_weight(real w) {objective_weights.emplace_back(w);}

    void resolve_weights(Design const &);
    void resolve_single_complex(ComplexWeights &cws, uint index, Weight const &w);
    void make_reversed_complexes(Design const &, vec<uint> const &);

    explicit operator bool() const {return len(specifications) > 0;}

    NUPACK_REFLECT(Weights, specifications, per_complex, per_tube, reversed_complexes, objective_weights);
};

}}
