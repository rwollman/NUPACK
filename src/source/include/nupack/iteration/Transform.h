/**
 * @brief Common algorithms for transforming, reducing ranges
 *
 * @file Transform.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../algorithms/Traits.h"
#include "../algorithms/Functor.h"
#include "../algorithms/Operators.h"
#include "../algorithms/Numeric.h"

#include <vector>
#include <algorithm>
#include <numeric>

namespace nupack {

/******************************************************************************************/

template <class V, NUPACK_IF(has_reserve<V>)> void reserve_space(V &v, std::size_t t) {v.reserve(t);}
template <class V, NUPACK_IF(!has_reserve<V>)> void reserve_space(V const &, std::size_t) {}

template <class V> V reserved(std::size_t n) {V v; reserve_space(v, n); return v;}

/******************************************************************************************/

/// Accumulate container with optional unary predicate
template <class V, class B, class F=Identity>
auto accumulate(V const &v, B const &update, F const &f={}) -> no_qual<decltype(f(*begin_of(v)))> {
    if (begin_of(v) == end_of(v)) return no_qual<decltype(f(*begin_of(v)))>{zero};
    auto ret = f(*begin_of(v));
    for (auto it = std::next(begin_of(v)); it != end_of(v); ++it) update(ret, f(*it));
    return ret;
}

/// Sum of container with optional unary predicate
template <class V, class F=Identity>
auto sum(V const &v, F const &f={}) -> no_qual<decltype(f(*begin_of(v)))> {
    return accumulate(v, plus_eq, f);
}

/// Product of container with optional unary predicate
template <class V, class F=Identity>
auto product(V const &v, F const &f={}) -> no_qual<decltype(f(*begin_of(v)))> {
    return accumulate(v, times_eq, f);
}

/// Prefix sums of a container: if keep_first is true, then the returned vector starts with a 0
template <class Out=void, class V, class F=plus_t>
auto prefixes(bool keep_first, V const &v, F const &f={}) {
    nonvoid<Out, std::vector<value_type_of<V>>> ret(static_cast<std::size_t>(keep_first), value_type_of<V>{});
    std::partial_sum(begin_of(v), end_of(v), std::back_inserter(ret), f);
    return ret;
}

/******************************************************************************************/

template <class U, class V> bool equal_range(U const &u, V const &v) {
    return std::equal(begin_of(u), end_of(u), begin_of(v), end_of(v));
}

/******************************************************************************************/

template <class V, class F=less_t> void sort(V &&v, F &&f={}) {std::sort(begin_of(v), end_of(v), fw<F>(f));}

template <class V, class F=less_t> V sorted(V v, F &&f={}) {sort(v, fw<F>(f)); return v;}

template <class V, class F=less_t> bool is_sorted(V const &v, F &&f={}) {return std::is_sorted(begin_of(v), end_of(v), fw<F>(f));}

template <class R=void, class V, class F=less_t> auto arg_sort(V const &v, F &&f={}) {
    nonvoid<R, std::vector<size_type_of<V>>> out{indices(v)};
    sort(out, [&](auto i, auto j) {return f(v[i], v[j]);});
    return out;
}

template <class R=void, class V, class F=less_t>
auto iter_sort(V &&v, F &&f={}) {
    nonvoid<R, std::vector<decltype(begin_of(v))>> out{iterators(v)};
    sort(out, [&](auto i, auto j) {return f(*i, *j);});
    return out;
}
/******************************************************************************************/

/// Swap iterator in container with the back one, and erase the new back one
template <class V> void swap_erase(V &v, iterator_of<V> it) {
    if (std::next(it) != end_of(v)) *it = std::move(v.back());
    v.pop_back();
}

/// Swap element in container with the back one, and erase the new back one
template <class V> void swap_erase(V &v, size_type_of<V> i) {swap_erase(v, begin_of(v) + i);}

/******************************************************************************************/

template <class V, class F, class T>
void replace(V &&v, F &&f, T &&t) {std::replace(begin_of(v), end_of(v), fw<F>(f), fw<T>(t));}

template <class V, class ...Ts>
void replace_if(V &&v, Ts &&...ts) {std::replace_if(begin_of(v), end_of(v), fw<Ts>(ts)...);}

template <class V, class F>
void transform(V &&v, F &&f) {std::transform(begin_of(v), end_of(v), begin_of(v), fw<F>(f));}

template <class V, class V2, class F>
void transform(V &&v, V2 &&v2, F &&f) {std::transform(begin_of(v), end_of(v), begin_of(v2), fw<F>(f));}

/******************************************************************************************/

/// Take first from container and erase that element
template <class V> auto take_first(V &v) {
    auto it = begin_of(v); auto ret = *it;
    v.erase(it); return ret;
}

/******************************************************************************************/

/// Insert at beginning of container
template <class V>
void insert_front(V &v, value_type_of<V> const &t) {v.insert(begin_of(v), t);}

/******************************************************************************************/

/// Insert at end_of of container
template <class V, class... Ts>
auto extend(V &v, Ts &&...ts) -> decltype(v.insert(end_of(v), fw<Ts>(ts)...)) {
    return v.insert(end_of(v), fw<Ts>(ts)...);
}

template <class V, class V2>
auto extend(V &v, V2 &&v2) -> decltype(v.insert(end_of(v), begin_of(fw<V2>(v2)), end_of(fw<V2>(v2)))) {
    return v.insert(end_of(v), begin_of(fw<V2>(v2)), end_of(fw<V2>(v2)));
}

template <class Out=void, class V>
auto join(V &&v) {
    nonvoid<Out, value_type_of<V>> ret;
    reserve_space(ret, sum(v, len));
    for (auto &&i : v) extend(ret, fw<decltype(i)>(i));
    return ret;
}

template <template <class...> class Out, class V>
auto join(V &&v) {
    return join<Out<value_type_of<value_type_of<V>>>>(fw<V>(v));
}

/// Concatenate consecutive values, containers, or iterator pairs into a container
template <class V> void cat(V const &) {}

/// Concatenate consecutive values, containers, or iterator pairs into a container
template <class V, class... Args>
void cat(V &v, const_iterator_of<V> b, const_iterator_of<V> e, Args &&...args) {
    extend(v, b, e); cat(v, fw<Args>(args)...);
}

/// Concatenate consecutive values, containers, or iterator pairs into a container
template <class V, class... Args>
void cat(V &v, V const &other, Args &&...args) {
    extend(v, begin_of(other), end_of(other)); cat(v, fw<Args>(args)...);
}

/// Concatenate consecutive values, containers, or iterator pairs into a container
template <class V, class... Args>
void cat(V &v, value_type_of<V> const &t, Args &&...args) {
    extend(v, t); cat(v, fw<Args>(args)...);
}

/// Concatenate consecutive values, containers, or iterator pairs into a container
template <class V, class... Args>
V catted(Args &&...args) {V v; cat(v, fw<Args>(args)...); return v;}

/******************************************************************************************/

/// Concatenate between two iterators skipping over the end_of() if necessary
template <class V1, class V2>
void circular_cat(V1 &to, V2 const &from, const_iterator_of<V2> i, const_iterator_of<V2> j) {
    if (i <= j) cat(to, i, j);
    else cat(to, i, end_of(from), begin_of(from), j);
}

/******************************************************************************************/

/// Rotate so the element with the minimum begin_of() in a range is at the beginning
template <class V> auto rotate_min_begin(V &v) {
    auto it = std::min_element(begin_of(v), end_of(v), [](auto &&s1, auto &&s2) {return begin_of(s1) < begin_of(s2);});
    std::rotate(begin_of(v), it, end_of(v));
    return it - begin_of(v);
}

/******************************************************************************************/

/// Rotate so the minimum element is at the beginning
template <class V> auto rotate_min(V &v) {
    auto it = std::min_element(begin_of(v), end_of(v));
    std::rotate(begin_of(v), it, end_of(v));
    return std::distance(begin_of(v), it);
}

/******************************************************************************************/

template <class V, class F>
auto erase_if(V &v, F const &f) {
    return v.erase(std::remove_if(begin_of(v), end_of(v), predicate_op<equals_t>(f)), end_of(v));
}

/******************************************************************************************/

/// Fill a range with a value
template <class V, class T>
auto fill(V &v, T const &t) -> decltype(v = t) {return v = t;}

/// Fill a range with a value
template <class V, class T, NUPACK_IF(!can_assign_from<V, T>)>
V & fill(V &v, T const &t) {for (auto &i : v) fill(i, t); return v;}

/******************************************************************************************/

template <class V, class F>
auto partition(V &v, F &&f) {return std::partition(begin_of(v), end_of(v), fw<F>(f));}

/******************************************************************************************/

template <class V, class F=less_t> V merge(V const &v1, V const &v2, F f={}) {
    V v;
    reserve_space(v, len(v1) + len(v2));
    std::merge(begin_of(v1), end_of(v1), begin_of(v2), end_of(v2), std::back_inserter(v), f);
    return v;
}

/******************************************************************************************/

template <class V, class F=less_t>
V take_lowest(V v, std::size_t n, F &&f={}) {
    if (n < len(v)) {
        std::partial_sort(begin_of(v), begin_of(v) + n, end_of(v), reduce_op<less_t>(fw<F>(f)));
        v.erase(std::next(begin_of(v), n), end_of(v));
    }
    return v;
}

/******************************************************************************************/

template <class V>
V unique_sorted(V v) {
    std::sort(begin_of(v), end_of(v));
    v.erase(std::unique(begin_of(v), end_of(v)), end_of(v));
    return v;
}

template <class V, class F>
V unique_sorted(V v, F &&f) {
    std::sort(begin_of(v), end_of(v), [&](auto const &m, auto const &n) {return f(m) < f(n);});
    v.erase(std::unique(begin_of(v), end_of(v), [&](auto const &m, auto const &n) {return f(m) == f(n);}), end_of(v));
    return v;
}

/******************************************************************************************/

template <class V> V duplicate(V const &v, std::size_t n=2) {
    V ret;
    reserve_space(ret, n * len(v));
    for (std::size_t i = 0; i != n; ++i) cat(ret, v);
    return ret;
}

template <class V> bool is_duplicate(V const &v, std::size_t n=2) {
    auto m_mod = div(len(v), n);
    if (m_mod.second) return false;
    std::size_t m{m_mod.first};
    auto const b = begin_of(v), e = end_of(v), c = b + m;
    auto p = e, q = p + m;
    for (; p != e; p += m, q += m)
        if (!std::equal(b, c, p, q)) return false;
    return true;
}

/******************************************************************************************/

template <class ...V> void erase_all(V &...v) {(void) std::initializer_list<int>{(swap(V(), v), 0)...};}

/******************************************************************************************/

NUPACK_DETECT_2(can_emplace_back, decltype(declval<T>().emplace_back(declval<U>())));

template <class R, class V, class F, class P=AlwaysTrue, NUPACK_IF(!can_emplace_back<R, decltype(declval<F>()(*begin_of(declval<V>())))>)>
auto vmap(V &&v, F &&map, P &&predicate={}) {
    R ret(len(v));
    auto i = begin_of(ret);
    for (auto j = begin_of(v); j != end_of(v); ++i, ++j) if (predicate(*j)) *i = map(*j);
    return ret;
}

/// Apply a functor f across a range v and collect the results in a container
template <class R, class V, class F, class P=AlwaysTrue, NUPACK_IF(can_emplace_back<R, decltype(declval<F>()(*begin_of(declval<V>())))>)>
auto vmap(V &&v, F &&map, P &&predicate={}) {
    R ret;
    reserve_space(ret, len(v));
    for (auto i = begin_of(v); i != end_of(v); ++i) if (predicate(*i)) ret.emplace_back(map(*i));
    return ret;
}

template <template <class...> class R=std::vector, class V, class F, class P=AlwaysTrue>
auto vmap(V &&v, F &&map, P &&predicate={}) {
    return vmap<R< no_qual<decltype(map(*begin_of(v)))> >>(fw<V>(v), fw<F>(map), fw<P>(predicate));
}

/// For a range of indices is... and a container v, return [v[i]...] as a new container
template <class R, class I, class V> auto imap(I const &is, V const &v) {
    return vmap<R>(is, [&](auto i) -> decltype(v[i]) {return v[i];});
}

template <class I, class V> auto imap(I const &is, V const &v) {return imap<V>(is, v);}

template <template <class ...> class R=std::vector, class T>
auto vfull(std::size_t n, T const &t) {return R<T>(n, t);}

/******************************************************************************************/

/// Copy first range into second range
template <class T, class U>
auto copy_range(T const &t, U &&u) {return std::copy(begin_of(t), end_of(t), begin_of(u));}

/******************************************************************************************/

}
