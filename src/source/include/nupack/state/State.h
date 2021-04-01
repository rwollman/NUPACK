#pragma once
#include "StaticState.h"
#include "../model/Model.h"

namespace nupack {

/******************************************************************************************/

template <class S, class Model>
auto structure_energy(S &&sequences, PairList p, Model const &em, bool distinguishable=true) {
    StaticState<> w(std::forward<S>(sequences), std::move(p));
    w.check_structure(em.pairable);
    auto out = w.calculate_energy(em);

    if (!distinguishable) out += std::log(real(w.symmetry())) / em.beta;

    return out;
}

/******************************************************************************************/

}
