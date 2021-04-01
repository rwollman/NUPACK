/** \file StaticState.h
 * @brief Contains StaticState, containing all State functionality excluding kinetics
 */
#pragma once
#include "StateBase.h"
#include "../model/Model.h"
#include "../loop/StaticLoop.h"
#include "../standard/Set.h"
#include "../standard/Ptr.h"

namespace nupack {

template <class Base_, class Loop_>
struct StaticState : public Base_, public Indexable<StaticState<Base_, Loop_>> {

    using Loop = Loop_;
    using base_type = Base_;
    using base_type::complexes;
    using base_type::pairs;
    using base_type::sys;

    StaticState() = default;

    template <class Y=vec<string>, class ...Ts>
    explicit StaticState(Y &&sys, PairList p, Ts const &...ts)
        : base_type(fw<Y>(sys), std::move(p)) {if (!pairs.empty()) build(ts...);}

    /**************************************************************************************/

    /// Unordered vector of loops
    vec<Loop> loops;
    NUPACK_EXTEND_REFLECT(StaticState, base_type, loops);

    /**************************************************************************************/

    /// Iterate through the loops
    auto & iter() {return loops;};

    bool check_structure(Pairable, bool throw_error) const;

    /// Check that the structure is valid for the user-supplied energy model
    void check_structure(Pairable p) const {check_structure(p, true);}
    bool is_valid(Pairable p) const {return check_structure(p, false);}

    /// Calculate the energy w.r.t. to a user-supplied energy model
    template <class Model> auto calculate_energy(Model const &) const;
    /// Build all the loops
    template <class ...Ts> void build(Ts const &...);
    /// Add a loop
    template <class ...Ts> void emplace(Ts&& ...ts) {loops.emplace_back(fw<Ts>(ts)...);}

    auto with_structure(PairList p={}) const {return StaticState(this->sys, std::move(p));}
    auto with_structure(StateBase const &w) const {return with_structure(this->aligned_pairs(w));}
};

NUPACK_DEFINE_TEMPLATE(isState, StaticState, class, class);

/******************************************************************************************/

template <class W, class O> template <class ...Ts>
void StaticState<W, O>::build(Ts const &...ts) {
    NUPACK_ASSERT(sys, "Strand system pointer should not be null");
    loops.reserve(pairs.n_pairs() + len(*sys));
    complexes = build_complex_set(loops, *sys, pairs, ts...);
    for (auto &o : loops) {
        o.finalize();
        if (o.exterior()) complexes.set_loop_index(o.strand_index(*this), o.index());
    }
}

/******************************************************************************************/

template <class W, class O>
bool StaticState<W, O>::check_structure(Pairable p, bool throw_error) const {
    bool ok = true;
    pairs.for_each_pair([&](iseq i, iseq j) {
        auto const strand1 = sys->strand_map[i], strand2 = sys->strand_map[j];
        if (!p(strand1 != strand2, next(sys->total_sequence, i), next(sys->total_sequence, j))) {
            ok = false;

            if (throw_error) {
                auto const dotparens = this->dp();
                auto const sequences = this->sequence();
                auto const index1 = i - sys->nicks[strand1] - 1, index2 = j - sys->nicks[strand2] - 1;
                Base const base1 = sys->total_sequence[i], base2 = sys->total_sequence[j];
                NUPACK_ERROR("Invalid secondary structure for the given energy model",
                            sequences, dotparens, strand1, strand2, index1, index2, base1, base2);
            }
        }
    });
    return ok;
}

/// Calculate the energy for any Model and return it, not modifying *this
template <class W, class O> template <class Model>
auto StaticState<W, O>::calculate_energy(Model const &model) const {
    return (len(sys->strands) - len(complexes)) * model.join_penalty()
        + sum(loops, [&](auto const &o){return model.loop_energy(o.sequences(), o.nick());});
}

/******************************************************************************************/

pair_data_type make_pairs_map(StateBase const &w);

/******************************************************************************************/

template <class State, class F>
void for_exterior_loops_in_complex(State &&w, ComplexSet::Indices const &v, F &&f) {
    for (auto &&o : w) if (o.exterior() && contains(v, o.strand_index(w))) f(o);
}

/******************************************************************************************/

}
