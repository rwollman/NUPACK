#pragma once
#include "SequenceAdapter.h"
#include "TypeImports.h"
#include "../types/Structure.h"

namespace nupack { namespace newdesign {

using SplitPoint = std::pair<uint, uint>;

ProbabilityMatrix sparsify(Tensor<real, 2> const &in, real f_sparse);

struct ProbabilitySplit {
    uint first;
    uint second;
    real prob;
    real cost;
    NUPACK_REFLECT(ProbabilitySplit, first, second, prob, cost);

    ProbabilitySplit() = default;
    ProbabilitySplit(uint f, uint s, real p, real c) : first(f), second(s), prob(p), cost(c) {}
};

/**
 * @brief check if two split points lead to disjoint ensembles.
 * @details the function will return the same answer regardless of argument order.
 *
 * @param a one split
 * @param b the other split
 *
 * @return true if the two split points cross
 */
template <class A=SplitPoint, class B=SplitPoint>
bool crosses(A const &a, B const &b) {
    int i, j, d, e;
    std::tie(i, j) = std::minmax(a.first, a.second);
    std::tie(d, e) = std::minmax(b.first, b.second);
    if (i == d && j == e) return false;
    // return (i <= d && d < j && j <= e) || (d <= i && i < e && e <= j);
    return d == j || e == i || (i <= d && d < j && j <= e) || (d <= i && i < e && e <= j);
}


vec<SplitPoint> valid_split_points(Structure const &s, uint min_size, uint min_helix);
bool is_padded(int i, Nicks const &bounds, int min_helix);
bool is_large_enough(SplitPoint sp, int n, int min_size);
bool is_valid(SplitPoint sp, Structure const &s, uint min_size, uint min_helix);
vec<SplitPoint> ascending_cost_splits(vec<SplitPoint> splits, uint n);
real children_cost(SplitPoint, uint n);
std::tuple<vec<ProbabilitySplit>, vec<ProbabilitySplit>>
        possible_splits(ProbabilityMatrix const &probs, uint min_size, uint min_helix, Structure const &s);
vec<SplitPoint> minimal_splits(ProbabilityMatrix const &probs, real f_split, uint min_size, uint min_helix, Structure const &s);
std::pair<Structure, Structure> split(SplitPoint const &sp, Structure const &s);
std::pair<vec<StrandView>, vec<StrandView>> split(SplitPoint const &sp, vec<StrandView> const &seq);
std::pair<vec<SplitPoint>, vec<SplitPoint>> split(SplitPoint const &sp, vec<SplitPoint> const &pairs);
}}
