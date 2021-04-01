/**
 * @brief Non-backtracking functions to calculate pair probability and MFE gap matrix
 *
 * @file PairProbability.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Tensor.h"

namespace nupack { namespace thermo {

/// Return the pair probability matrix with the unpaired probability on the diagonal
/// Assumed that the lower half QB is nonzero
template <class Out, class Mat, class T>
auto pairs_from_QB(PF, T const q, Mat const &QB) {
    auto n = len(QB) / 2;
    Tensor<Out, 2> PP(n, n, *zero);
    if (mantissa(q) == 0) return PP;
    NUPACK_ASSERT(std::isfinite(mantissa(q)), q);
    auto const iq = PF::invert()(q);
    for (auto i : range(n)) for (auto j : range(i+1, n)) {
        bool err = false;
        *PP(i, j) = PF::element_value(err, fold(PF::times(), *QB(j, i), iq, *QB(i + n, j)), Zero());
        NUPACK_ASSERT(!err, "Overflow during pair probability calculation", *QB(j, i), iq, *QB(i + n, j));
        *PP(j, i) = *PP(i, j);
    }
    for (auto i : range(n)) *PP(i, i) = 1 - sum(PP(i, span(0, n)));
    return PP;
}

/// Return the base pair MFE cost matrix with the minimum cost on the diagonal
/// Cost(i, j) = mfe(given i and j paired) - mfe(no constraints)
template <class Out, class Mat, class T>
auto pairs_from_QB(MFE, T const q, Mat const &QB) {
    auto n = len(QB) / 2;
    Tensor<Out, 2> PP(n, n, *inf);
    for (auto i : range(n)) for (auto j : range(i+1, n))
        *PP(j, i) = *PP(i, j) = *QB(j, i) + *QB(i + n, j) - q;
    for (auto i : range(n)) *PP(i, i) = minimum(PP(i, span(0, n)));
    return PP;
}

/******************************************************************************************/

}}
