/**
 * @brief Spreadsort algorithm from Boost, plus a binary radix sort
 *
 * @file Spreadsort.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../algorithms/Traits.h"
#include "../algorithms/Functor.h"
#include "../algorithms/Numeric.h"

#include <boost/sort/spreadsort/spreadsort.hpp>

namespace nupack {

/******************************************************************************************/

template <class V, class R, class ...Ts>
void spreadsort_float(V &&v, R &&r, Ts &&...ts) {
    auto f = [r=fw<R>(r)](auto &&x, auto offset) {
        using T = no_qual<decltype(r(fw<decltype(x)>(x)))>;
        return boost::sort::spreadsort::float_mem_cast<T, uint_of_size<sizeof(T)>>(r(fw<decltype(x)>(x))) >> offset;
    };
    boost::sort::spreadsort::float_sort(begin_of(v), end_of(v), std::move(f), fw<Ts>(ts)...);
}

template <class V, NUPACK_IF(is_floating_point<value_type_of<no_qual<V>>>)>
void spreadsort(V &&v) {boost::sort::spreadsort::float_sort(begin_of(v), end_of(v));}

template <class V, class ...Ts, NUPACK_IF(is_integral<value_type_of<no_qual<V>>>)>
void spreadsort(V &&v, Ts &&...ts) {boost::sort::spreadsort::integer_sort(begin_of(v), end_of(v), fw<Ts>(ts)...);}

/******************************************************************************************/

/// Simple binary radix sort - most significant digit first, recursive
template <class Iter>
void binary_sort_ms(Iter b, Iter e, std::size_t const i=0) {
    auto it = std::partition(b, e, [n=bitsof(*b)-1-i](double const &x) {return !bit_at(x, n);});
    if (i+1 != bitsof(*b) && b != it) binary_sort_ms(b, it, i+1);
    if (i+1 != bitsof(*b) && e != it) binary_sort_ms(it, e, i+1);
}

/// Simple binary radix sort - most significant digit first, not recursive
template <class Iter>
void binary_sort_ms2(Iter b, Iter e) {
    vec<std::tuple<Iter, Iter, std::size_t>> v;
    v.reserve(e-b);
    v.emplace_back(b, e, bitsof(*b)-1);
    while (!v.empty()) {
        auto const i = third_of(v.back());
        auto const p = first_of(v.back()), q = second_of(v.back());
        v.pop_back();
        auto const m = std::partition(p, q, [i](double const &x) {return !bit_at(x, i);});
        if (i && p != m) v.emplace_back(p, m, i-1);
        if (i && q != m) v.emplace_back(m, q, i-1);
    }
}

/// Simple binary radix sort - least significant digit first
template <class Iter>
void binary_sort_ls(Iter b, Iter e) {
    for (auto i : range(bitsof(*b)))
        std::stable_partition(b, e, [i](double const &x) {return !bit_at(x, i);});
}

/******************************************************************************************/

}
