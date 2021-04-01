/**
 * @brief Common searching, counting, and finding algorithms
 *
 * @file Search.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Iterator.h"

#include "../algorithms/Utility.h"
#include "../algorithms/Operators.h"
#include "../algorithms/Tuple.h"
#include "../common/Constants.h"
#include <stdexcept>

namespace nupack {

/******************************************************************************************/

/// Count number of matches to element "t" in container "v"
template <class V, class T>
auto count(V const &v, T &&t) {return std::count(begin_of(v), end_of(v), fw<T>(t));}

/// Count number of elements for which predicate returns true
template <class V, class F=Identity>
auto count_if(V const &v, F &&f={}) {return std::count_if(begin_of(v), end_of(v), fw<F>(f));}

/******************************************************************************************/

template <class V, class F>
void collect(V &v, F &&f) {for (auto &i : v) i = f();}

/******************************************************************************************/

template <class V, class F=Identity> bool any_of(V const &v, F &&f={}) {return std::any_of(begin_of(v), end_of(v), fw<F>(f));}
template <class V, class F=Identity> bool all_of(V const &v, F &&f={}) {return std::all_of(begin_of(v), end_of(v), fw<F>(f));}
template <class V, class F=Identity> bool none_of(V const &v, F &&f={}) {return std::none_of(begin_of(v), end_of(v), fw<F>(f));}

/// Find index of an element in a container
template <class V, class T> usize find_index(V const &vec, T &&t) {
    return std::find(begin_of(vec), end_of(vec), fw<T>(t)) - begin_of(vec);
}

/******************************************************************************************/

/// Test if an iterator is in a random-access container
template <class V, NUPACK_IF(is_random_access<V>)>
bool contains_iter(V const &v, const_iterator_of<V> t) {
    return (t <= end_of(v)) && (t >= begin_of(v));
}

template <class V, NUPACK_IF(!is_random_access<V>)>
bool contains_iter(V const &v, const_iterator_of<V> t) {
    auto const e = end_of(v);
    if (t == e) return true;
    for (auto b = begin_of(v); b != e; ++b) if (t == b) return true;
    return false;
}

/******************************************************************************************/

/// Test if a value is in a container
template <class V>
bool contains(V const &v, value_type_of<V> const &t) {
    return std::find(begin_of(v), end_of(v), t) != end_of(v);
}

/******************************************************************************************/

/// Binary search between two iterators
template <class It1, class It2, class T, class F=Identity>
It1 binary_it_search(It1 b, It2 e, T const &t, F const &f=Identity()) {
   auto mask = e - b - 1; auto ret = b;
   while (mask > 0) if ((mask >>= 1) + ret < e && f(ret[mask]) < t) ret += mask + 1;
   return ret;
}

/******************************************************************************************/

template <class V, class F=not_equals_t>
bool all_equal(V const &v, F &&f={}) {
    return std::adjacent_find(begin_of(v), end_of(v), fw<F>(f)) == end_of(v);
}

/******************************************************************************************/

template <class V, class Iter>
auto nonconst_iter(V &v, Iter it) -> decltype(v.erase(it, it)) {return v.erase(it, it);}

/******************************************************************************************/

template <class V, class T> auto upper_bound(V &&v, T const &t) {
    return std::upper_bound(begin_of(v), end_of(v), t);
}
template <class V, class T, class F> auto upper_bound(V &&v, T const &t, F &&f) {
    return std::upper_bound(begin_of(v), end_of(v), t, [&](auto const &a, auto const &b){
        return a < f(b);
    });
}

template <class V, class T> auto lower_bound(V &&v, T const &t) {
    return std::lower_bound(begin_of(v), end_of(v), t);
}
template <class V, class T, class F> auto lower_bound(V &&v, T const &t, F &&f) {
    return std::lower_bound(begin_of(v), end_of(v), t, [&](auto const &a, auto const &b){
        return f(a) < b;
    });
}

/// Binary search in a container
template <class V, class T, class F=Identity>
auto binary_search(V &&v, T const &t, F const &f=Identity()) -> decltype(begin_of(v)) {
   return binary_it_search(begin_of(v), end_of(v), t, f);
}

/// Binary search in a container
template <class V, class ...Ts>
auto binary_search_index(V &&v, Ts &&...ts) -> decltype(end_of(v) - begin_of(v)) {
   return binary_search(v, fw<Ts>(ts)...) - begin_of(v);
}

/******************************************************************************************/

/// Find the first index where the partial sum of it is greater than a value
template <class It1, class It2, class T, class F=Identity>
std::pair<It2, T> find_it_cumulative(It1 b, It2 e, T t, F const &f=Identity()) {
    for (auto i = b; i != e; ++i) if (minus_if(t, f(*i))) return {i, t};
    throw std::out_of_range("find_it_cumulative");
}

/******************************************************************************************/

/// Find the first index where the partial sum of it is greater than a value
template <class V, class T, class F=Identity>
auto find_cumulative(V const &v, T t, F const &f=Identity()) {
    return find_it_cumulative(begin_of(v), end_of(v), t, f);
}

/******************************************************************************************/

/// Find a value in a container where only indices where mask is true are considered
template <class V, class T, class M>
auto find_with_mask(V const &v, T const &t, M const &mask) {
    auto it = begin_of(v);
    auto m = begin_of(mask);
    while (!(it == end_of(v) || (*m && *it == t))) {++m; ++it;}
    return std::make_pair(it, m);
}

/******************************************************************************************/

/// Find cumulative for a 2D container, with some forgiveness for precision issues
template <class M, class T>
std::tuple<usize, usize, T> find_cumulative_mat(M const &m, T t) {
    usize b = 0, c = 0;
    for (b = 0; b != 4; ++b) for (c = 0; c != 4; ++c)
        if (minus_if(t, m(b, c))) return std::make_tuple(b, c, t);
    throw std::out_of_range("find_cumulative_mat");
}

/******************************************************************************************/

template <class It1, class It2, class F=Identity>
auto find_first_mismatch(It1 b, It2 e, F &&f=Identity()) {
    auto ret = std::adjacent_find(b, e, [=](auto const &i, auto const &j){return f(i) != f(j);});
    return ret == e ? ret : std::next(ret);
}

/******************************************************************************************/

template <class V, class F, class Comp=less_t>
auto extremum(V &&v, F &&f, Comp &&comp={}) {
    auto it = begin_of(v);
    auto const e = end_of(v);
    auto ret = f(*it);
    while (++it != e) {auto tmp = f(*it); if (comp(tmp, ret)) ret = tmp;}
    return ret;
}

template <class V, class F=Identity> auto max_element(V &&v, F &&f=Identity()) {
    return std::max_element(begin_of(v), end_of(v), reduce_op<less_t>(fw<F>(f)));
}

template <class V, class F=Identity> decltype(auto) maximum(V &&v, F &&f=Identity()) {
    return *max_element(fw<V>(v), fw<F>(f));
}

/******************************************************************************************/

template <class V, class F=Identity> auto min_element(V &&v, F &&f=Identity()) {
    return std::min_element(begin_of(v), end_of(v), reduce_op<less_t>(fw<F>(f)));
}

template <class V, class F=Identity> decltype(auto) minimum(V &&v, F &&f=Identity()) {
    return *min_element(fw<V>(v), fw<F>(f));
}

template <class V, class F> auto max_value(V &&v, F &&f) {
    auto b = begin_of(v);
    auto const e = end_of(v);
    auto ret = f(*b);
    while (++b != e) {auto tmp = f(*b); if (tmp > ret) ret = tmp;}
    return ret;
}

/******************************************************************************************/

template <class T> auto minmax(T const &t) {
    return unpack_as<T>(std::minmax(std::get<0>(t), std::get<1>(t)));
}

/******************************************************************************************/

template <class V, class F=Identity> auto min_index(V const &v, F &&f=Identity()) {
    return min_element(v, fw<F>(f)) - begin_of(v);
}

template <class V, class F=Identity> auto max_index(V const &v, F &&f=Identity()) {
    return max_element(v, fw<F>(f)) - begin_of(v);
}

/******************************************************************************************/

template <class V, class F> auto find_if(V &&v, F &&f) {
    return std::find_if(begin_of(v), end_of(v), fw<F>(f));
}

/******************************************************************************************/

template <class V, class F>
auto log_scan(V &&v, F &&f) {
    auto n = len(v);
    while (n >>= 1) if (f(view(end_of(v) - 2 * n, end_of(v)))) return n;
    return n;
}

/******************************************************************************************/

template <class V, class F>
auto last2_scan(V &&v, F &&f) {return f(view(end_of(v) - 2, end_of(v)));}

/******************************************************************************************/

template <class V, class ...Ts, NUPACK_IF(is_findable<V, Ts...>)>
auto find(V &&v, Ts &&...ts) {return v.find(fw<Ts>(ts)...);}

template <class V, class ...Ts, NUPACK_IF(!is_findable<V, Ts...>)>
auto find(V &&v, Ts &&...ts) {
    return std::find(begin_of(v), end_of(v), fw<Ts>(ts)...);
}

/******************************************************************************************/

template <class V, class ...Ts, NUPACK_IF(is_findable<V, Ts...>)>
auto ordered_find(V &&v, Ts &&...ts) {return v.find(fw<Ts>(ts)...);}

template <class V, class ...Ts, NUPACK_IF(!is_findable<V, Ts...>)>
auto ordered_find(V &&v, Ts &&...ts) {
    return std::lower_bound(begin_of(v), end_of(v), fw<Ts>(ts)...);
}

template <class V, class ...Ts>
auto ordered_index(V &&v, Ts &&...ts) {
    auto const b = begin_of(v);
    return ordered_find(fw<V>(v), fw<Ts>(ts)...) - b;
}

template <class U, class V, class F=std::equal_to<>>
constexpr bool equal_ranges(U const &u, V const &v, F &&f={}) {
    return std::size(u) == std::size(v) && std::equal(begin_of(u), end_of(u), begin_of(v), fw<F>(f));
}

/******************************************************************************************/

}
