#pragma once
#include "../standard/Vec.h"

namespace nupack {

// Cost of calculating the lengths[:-1] x lengths[1:] subblock
// Reduces to n**3 if len(lengths) is 1
template <class V>
std::size_t subblock_cost(V const &v) {
    if (len(v) == 0) return 0;
    if (len(v) == 1) return cube(std::size_t(front(v)));
    return 3 * front(v) * back(v) * (2 * sum(v, caster<std::size_t>) - front(v) - back(v));
}

inline constexpr std::size_t unit_subblock_cost(std::size_t n) {return n > 1 ? 6 * (n - 1) : n;}

std::array<std::size_t, 3> unit_evaluation_costs(uint n, std::size_t lmax);

using EvaluationCostTable = vec<vec<std::array<std::size_t, 3>>>;
EvaluationCostTable unit_evaluation_cost_table(uint n, real timeout);

}