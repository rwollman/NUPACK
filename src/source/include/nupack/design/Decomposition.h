#pragma once
#include "Split.h"
#include "Models.h"
#include "ThermoWrapper.h"
#include "Defect.h"
#include "Granularity.h"
#include "Logging.h"

namespace nupack { namespace newdesign {

struct DecompositionParameters {
    int H_split {2};
    int N_split {12};
    real f_split {0.99};
    real f_sparse {0.00001};
    real dG_clamp {-20};

    NUPACK_REFLECT(DecompositionParameters, H_split, N_split, f_split, f_sparse, dG_clamp);
};

/**
 * @brief replace a shared environment with a serial one if the size of the
 * complex (or subcomplex) is smaller than the threshold.
 *
 * @tparam V environment type
 * @tparam S a complex or decomposition node type
 * @param v an environment
 * @param s either a complex or a decomposition node
 */
template <class V, class S>
auto threshold(V const &v, S const &s) {
    static constexpr uint thresh = 500;
    if (v.n_workers() == 1) return v;
    if (len(s) <= thresh) return Local();
    return v;
}

struct ComplexNode;
using PairedChildren = std::pair<SplitPoint, std::pair<ComplexNode, ComplexNode>>;

using ThermoData = std::pair<ProbabilityMatrix, real>;
ThermoData join_children(SplitPoint, ThermoData const &, ThermoData const &);
ThermoData merge_alternatives(vec<ThermoData> const &, real f_sparse);


struct ComplexNode {
    struct Cache {
        mutable vec<std::pair<::nupack::Complex, ThermoData>> map;

        void add(::nupack::Complex seq, ThermoData data, uint depth) const {
            if (depth + 1 > len(map)) map.resize(depth + 1);
            at(map, depth) = {std::move(seq), std::move(data)};
        }

        bool match(::nupack::Complex const &s, uint depth) const {
            return len(map) > depth && at(map, depth).first == s;
        }

        void revoke_non_root() const { map.resize(1); }

        ThermoData get(uint depth) const {return at(map, depth).second;}

        NUPACK_REFLECT(Cache, map);
    };

    vec<SplitPoint> enforced_pairs;
    vec<StrandView> sequence;
    Structure structure;
    vec<PairedChildren> children;
    int index = -1; // unset value
    Cache cache;

    NUPACK_REFLECT(ComplexNode, enforced_pairs, structure, sequence, children, index, cache);

    ComplexNode() = default;
    ComplexNode(vec<StrandView> sequence, Structure structure, vec<SplitPoint> enforced) :
            enforced_pairs(std::move(enforced)),
            sequence(std::move(sequence)),
            structure(std::move(structure)) {}

    // pfunc -> merge or compute if end of tree or at that depth
    ThermoData dynamic_program(Local env, ThermoEnviron &t_env, Sequence const &s, uint depth,
            DecompositionParameters const &params, LevelSpecification const &indiv={},
            EngineObserver &obs=NullEngineObserver) const;

    void add_child(SplitPoint sp);
    void structure_decompose(uint min_size, uint min_helix);
    bool probability_decompose(DecompositionParameters const &params,
            Sequence const &s, ThermoEnviron &t_env, int depth=0,
            LevelSpecification const &indiv={}, EngineObserver &obs=NullEngineObserver);

    auto size() const {return sum(sequence, len);}
    uint depth() const;

    void register_indices(vec<int> &registered, uint depth, bool include_leaves=true) const;

    /**
     * @brief allows an arbitrary function to applied to all left and right
     *     children of a node
     *
     * @tparam F the functor type being passed in
     * @param f the functor passed in to call on all child nodes
     */
    template <class F>
    void child_op(F &&f) {
        for (auto &c : item_view(children)) {
            f(c.first); f(c.second);
        }
    }

    /**
     * @brief same as other child_op, but for non-modifying behavior
     */
    template <class F>
    void child_op(F const &f) const {
        for (auto const &c : item_view(children)) {
            f(c.first); f(c.second);
        }
    }

    /**
     * @brief apply an arbitrary function to this node and all descendent nodes.
     *
     * @tparam F the functor type
     * @param f the function to apply to the current nodes and all its descendents
     */
    template <class F>
    void apply_recursive(F &&f) {
        f(*this);
        child_op([&](auto &c) {c.apply_recursive(f);});
    }

    template <class F>
    void apply_recursive(F const &f) const {
        f(*this);
        child_op([&](auto const &c) {c.apply_recursive(f);});
    }

    /* code for JSON import/export without including the cache */
    static constexpr auto repr_names() {return make_names("structure", "sequence", "enforced_pairs", "index", "children", "length");}
    auto save_repr() const {return make_members(structure, sequence, enforced_pairs, index, children, size());}
    void load_repr(Structure str, vec<StrandView> s, vec<SplitPoint> spl, int ind, vec<PairedChildren> chi, uint) {
        structure = std::move(str);
        sequence = std::move(s);
        enforced_pairs = std::move(spl);
        index = std::move(ind);
        children = std::move(chi);
        cache = Cache();
    }
};


}}
