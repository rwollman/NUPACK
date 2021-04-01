#pragma once
#include "SequenceAdapter.h"
#include "../types/Structure.h"
#include "Models.h"
#include "Defect.h"
#include "Decomposition.h"
#include "TypeImports.h"
#include "../types/Sequence.h"
#include "Logging.h"

namespace nupack { namespace newdesign {

struct Target : MemberOrdered {
    Model<real> model;
    Structure structure;
    NUPACK_REFLECT(Target, structure, model);

    bool has_structure() const {return structure.valid();}

    /* Constructors */
    Target() = default;
    Target(Model<real> m, Structure s={}) : model(std::move(m)), structure(std::move(s)) {}

    auto const & cached_models(ModelMap const &map) const {return map.get(model).models;}
    auto & environment(ModelMap const &map) const {return map.get(model);}
};

struct Complex : MemberOrdered {
    // when linking to "global" sequence state, use defining string labels for
    // domains and strands to pull out appropriate views, which need not
    // change after that point. Any changes to this state will be
    // automatically reflected in the views without the Complex ever needing
    // to interface with the global sequence object again
    vec<StrandView> strands;
    Target target; // start with just one
    string name;
    DecompositionParameters params;
    ComplexNode decomposition;
    real bonus = 0;

    /* Most recently computed sequence-dependent values */
    NUPACK_REFLECT(Complex, name, strands, target, params, decomposition, bonus);

    /* Constructors */
    Complex() = default;
    Complex(vec<StrandView> s, Target t, string name={}, DecompositionParameters params={}, real bonus=0) :
            strands(std::move(s)), target(std::move(t)), name(std::move(name)), params(std::move(params)),
            decomposition(strands, target.structure, {}), bonus(bonus) {}

    /* Status */
    bool is_on_target() const {return target.has_structure();}

    /* Sequence dependent properties */
    real log_pfunc(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs=NullEngineObserver) const;
    Tensor<real, 2> pair_probabilities(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs=NullEngineObserver) const;
    Defect defect(Local env, ModelMap const &map, Sequence const &, EngineObserver &obs=NullEngineObserver) const;

    real join_penalty(ModelMap const &map) const {
        return -target.model.beta * (len(strands) - 1) * target.model.join_penalty();
    }

    /* Representational */
    vec<vec<uint>> strands_as_indices() const {
        return vmap(strands, [](auto const &s) {return vec<uint>(s.to_indices());});
    }

    vec<uint> to_indices() const { return join(strands_as_indices()); }

    vec<uint> nucleotide_counts() const {
        auto indices = to_indices();
        vec<uint> counts(maximum(indices)+1);
        for (auto i : indices) ++at(counts, i);
        return counts;
    }

    /* Estimate */
    real log_pf_single_strands(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs=NullEngineObserver) const;

    /* Decomposition */
    real log_pfunc(Local env, ModelMap const &map, Sequence const &s,
            uint depth, LevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver) const;
    ProbabilityMatrix pair_probabilities(Local env, ModelMap const &map, Sequence const &s,
            uint depth, LevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver) const;
    Defect defect(Local env, ModelMap const &map, Sequence const &s,
            uint depth, LevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver) const;

    void structure_decompose();
    bool probability_decompose(Sequence const &s, ModelMap const &map, uint depth=0,
            LevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver);
    auto depth() const {return decomposition.depth();}
    void index_nodes();
    vec<int> get_node_indices(uint depth, bool include_leaves=true) const;
    string hierarchical_pfunc(ModelMap const &map, Sequence const &s, uint depth, EngineObserver &obs=NullEngineObserver);
    string decomposition_connectivity() const;

    /**
     * @brief sum up strand nucleotides to get the number of nucleotides in the complex
     *
     * @return number of physical nucleotides in the complex
     */
    auto size() const {return sum(strands, len);}

    auto symmetry_correction() const {return std::log(rotational_symmetry(strands));}

    string json_decomposition() const;

};


}}
