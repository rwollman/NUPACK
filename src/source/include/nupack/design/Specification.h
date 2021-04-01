#pragma once
#include "Designer.h"
#include "Weights.h"

namespace nupack { namespace newdesign {

/** @brief specification convertible to CompConstraint or IdentConstraint */
struct DualListSpec {
    vec<string> left;
    vec<string> right;

    std::pair<vec<int>, vec<int>> get_variables(DesignSequence const &) const;

    NUPACK_REFLECT(DualListSpec, left, right);
};


/** @brief specification convertible to PatternConstraint */
struct PatternSpec {
    vec<string> domains;
    string pattern;

    void add_constraint(DesignSequence &) const;

    NUPACK_REFLECT(PatternSpec, domains, pattern);
};


struct DiversitySpec {
    vec<string> domains;
    int word_length;
    int min_nucleotide_types;

    void add_constraint(DesignSequence &) const;

    NUPACK_REFLECT(DiversitySpec, domains, word_length, min_nucleotide_types);
};


/** @brief specification convertible to WordConstraint */
struct WordSpec {
    vec<string> domains;
    vec<vec<string>> comparisons;

    void add_constraint(DesignSequence &) const;

    NUPACK_REFLECT(WordSpec, domains, comparisons);
};


/** @brief specification convertible to MatchConstraint */
struct SimilaritySpec {
    vec<string> domains;
    string reference;
    std::pair<real, real> range;

    void add_constraint(DesignSequence &) const;

    NUPACK_REFLECT(SimilaritySpec, domains, reference, range);
};


/** @brief specification of a Complex */
struct ComplexSpec {
    string name;
    vec<string> strands;
    Structure structure;

    NUPACK_REFLECT(ComplexSpec, name, strands, structure);
};


/** @brief specification of a Tube */
struct TubeSpec {
    string name;
    vec<std::pair<vec<string>, real>> targets;

    NUPACK_REFLECT(TubeSpec, name, targets);
};


/** @brief collection of specifications for individual constraint types */
struct ConstraintSpec {
    vec<DualListSpec> complementarity;
    vec<DualListSpec> match;
    vec<PatternSpec> pattern;
    vec<DiversitySpec> diversity;
    vec<WordSpec> word;
    vec<SimilaritySpec> similarity;

    NUPACK_REFLECT(ConstraintSpec, complementarity, match, pattern, diversity, word, similarity);
};


/** @brief specification of all components of a Design and its encapsulating Designer */
struct Specification {
    vec<DomainSpec> domains;
    vec<StrandSpec> strands;
    vec<ComplexSpec> complexes;
    vec<TubeSpec> tubes;
    Model<real> model;

    Weights weights;
    ConstraintSpec constraints;
    vec<Objective> objectives;

    DesignParameters parameters;
    bool wobble_mutations = false;

    Specification() = default;
    Specification(Model<real> m, bool w) : model(std::move(m)), wobble_mutations(w) {}

    NUPACK_REFLECT(Specification, domains, strands, complexes, tubes, model, weights, constraints, objectives, parameters, wobble_mutations)

    /** @brief look up complex by name or strands and return the index in the list
     * of complexes in the Design.
     */
    auto complex_index(vec<string> const &x) const {
        if (len(x) == 1) {
            auto it = find_if(complexes, [&](auto const &c) {return c.name == front(x);});
            if (it != end_of(complexes)) return it - begin_of(complexes);
        }

        auto low = lowest_rotation(x);
        auto it = find_if(complexes, [&](auto const &c) {return low == lowest_rotation(c.strands);});
        if (it == end_of(complexes)) NUPACK_ERROR("unknown complex", x);

        return it - begin_of(complexes);
    }

    static vec<uint> ensure_compatibility(Specification const &spec, SingleResult const &res);

    explicit operator Designer() const;
};

vec<int> extract_variables(vec<string> const &names, DesignSequence const &seqs);
vec<int> extract_element(string name, DesignSequence const &seqs);

}}
