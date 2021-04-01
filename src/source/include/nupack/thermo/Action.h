/**
 * @brief Recursion adapter for modifying QB recursion, used by design mostly
 *
 * @file Action.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../algorithms/TypeSupport.h"
#include <functional>

namespace nupack {namespace thermo {

/******************************************************************************************/

/**
 * Actions are called once for each base i, j during the calculation.
 * The boolean can_pair whether the thermodynamic model permits a pair between i and j
 * The f is the recursion yielding the value of QB
 * The third argument is the slice of sequences being calculated
 * The fourth and fifth arguments are the indices of i and j
 * The last argument is the partition function matrices block Q
 */
struct PairingAction {
    std::function<bool(int, int)> predicate; // treated as if it returns true if empty

    template <class Block, class Algebra, class F, class Model, class S>
    auto operator()(int i, int j, bool can_pair, Algebra A, Block const &, S const &s, Model const &, F &&recursion) const {
        bool bad = !can_pair || (predicate && !predicate(i + s.offset, j + s.offset));
        return bad ? A.zero() : A.maybe() & recursion();
    }
};

using DefaultAction = PairingAction;

/******************************************************************************************/

}}
