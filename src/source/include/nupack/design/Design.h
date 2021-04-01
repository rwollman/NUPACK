#pragma once
#include "SequenceAdapter.h"
#include "Complex.h"
#include "Tube.h"
#include "Models.h"
#include "Weights.h"
#include "Logging.h"

namespace nupack { namespace newdesign {

/**
 * @brief Contains the tangible elements of the design
 * @details [long description]
 *
 */
struct Design : MemberOrdered {
    DesignSequence sequences;
    vec<Complex> complexes; /// on- and off-targets
    vec<Tube> tubes;
    ModelMap models;
    NUPACK_REFLECT(Design, sequences, complexes, tubes, models);
    using is_member_comparable = True;

    /* Constructors */
    Design() = default;
    Design(DesignSequence seq) : sequences(std::move(seq)) {}

    /* add components to design */
    void add_complex(vec<string> const &strands, Model<real> model,
            string const &name={}, Structure struc={}, DecompositionParameters={}, real bonus=0);
    void add_tube(vec<uint> const &indices, vec<real> const &concs, string const &name);

    /* forward to sequences */
    void initialize_sequence() {sequences.initialize_sequence();}
    void set_sequence(Sequence const &s) {sequences.set_sequence(s);}
    bool mutate_sequence(vec<uint> const &vars) {return sequences.mutate_sequence(vars);}
    void add_structure_complementarity();

    Sequence const & sequence() const {return sequences.nucleotides;}

    vec<real> log_pfuncs(Local const &env, uint depth=0, EnsemblePartition const &part={},
            EnsembleLevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver) const;
    vec<Defect> complex_defects(Local const &env, uint depth=0, EnsemblePartition const &part={},
            EnsembleLevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver) const;

    Defect normalized_defect(Local const &env, uint depth=0, EnsemblePartition const &part={},
            EnsembleLevelSpecification const &indiv={}, Weights const &weights={}, EngineObserver &obs=NullEngineObserver) const;

    uint max_depth() const {return maximum(complexes, [](auto const &c) {return c.depth();}).depth();}
    void initialize_decomposition(EnsemblePartition const &part={});
    void redecompose_active(Local const &env, uint depth, EnsemblePartition const &part={});

    /* JSON serialization without the model since the whole point is deferred/cached creation */
    static constexpr auto repr_names() {return make_names("sequences", "complexes", "tubes");}
    auto save_repr() const {return make_members(sequences, complexes, tubes);}
    void load_repr(DesignSequence s, vec<Complex> c, vec<Tube> t) {
        sequences = std::move(s); complexes = std::move(c); tubes = std::move(t);
    }
};

uint find_tube(string name, Design const &design);
uint find_complex(string name, Design const &design);

Variant<DomainView, StrandView> find_sequence_element(Design const &, string const &);


}}
