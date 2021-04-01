/** \file StateBase.h
 * @brief Contains non-templated base classes for StaticState and JumpState objects
 */
#pragma once
#include "ComplexSet.h"
#include "System.h"
#include "../types/PairList.h"

namespace nupack {

/******************************************************************************************/

struct StateBase : TotallyOrdered {
    StateBase() {};

    explicit StateBase(std::shared_ptr<System const>, PairList p={});

    explicit StateBase(System s, PairList p={}) : StateBase(std::make_shared<System const>(std::move(s)), std::move(p)) {}

    NUPACK_REFLECT(StateBase, sys, pairs, complexes);

    /// Pointer to system for access to strands
    std::shared_ptr<System const> sys;
    /// Vector of which bases are paired
    PairList pairs;
    /// Complex data for the state
    ComplexSet complexes;
    /// Unpseudoknotted dot-parens
    string dp() const;
    /// Base sequences ordered to match dot-parens
    string sequence() const;
    /// Number of bases
    auto n_bases() const {return sys->n_bases();};
    /// Align another State's pairs to agree with the strand ordering of this one
    PairList aligned_pairs(StateBase const &) const;

    std::size_t symmetry() const;

    /// Print the state
    friend std::ostream & operator<<(std::ostream &os, StateBase const &w) {
        dump_os(os, "State('", w.sequence(), "', '", w.dp(), "')");
        return os;
    }

    /**************************************************************************************/

    bool operator<(StateBase const &w) const {
        if (less_ptr(sys, w.sys)) return true;
        if (less_ptr(w.sys, sys)) return false;
        return pairs < w.pairs;
    }

    bool operator==(StateBase const &w) const {
        return (pairs == w.pairs) && equal_ptr(sys, w.sys);
    }

    auto hash() const {return sys ? hash_of(*sys, pairs) : hash_of(pairs);}

    auto operator^(StateBase const &w) const {return pairs ^ w.pairs;}

    static constexpr auto repr_names() {return make_names("pairs", "sys");}

    auto save_repr() const {return make_members(pairs, sys);}

    void load_repr(decltype(pairs) p, decltype(sys) s) {
        if (s) *this = StateBase(std::move(s), std::move(p));
    }
};


/******************************************************************************************/

struct JumpStateBase : public StateBase {
    using base_type = StateBase;
    using StateBase::StateBase;

    NUPACK_EXTEND_REFLECT(JumpStateBase, base_type, last_move, add_rates, del_rates, energy);

    /// Last move the state took
    StateMove last_move;
    /// Base pair addition rates summed by loop
    Fenwick<long double> add_rates=Fenwick<long double>(0);
    /// Base pair deletion rates summed by loop
    Fenwick<long double> del_rates=Fenwick<long double>(0);
    /// Total free energy of state
    real energy=0;

    /// Register change in energy, base pair changes (called from Loop)
    void register_move(BaseIter, BaseIter, real dE, real rate);
};


/******************************************************************************************/

NUPACK_DEFINE_TYPE(is_state_base, StateBase);
NUPACK_EXTEND_TYPE(is_state_base, JumpStateBase);

}
