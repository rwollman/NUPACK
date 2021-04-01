/**
 * @brief Common iteration patterns including compile-time traversal
 *
 * @file Patterns.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Transform.h"
#include "../algorithms/Tuple.h"
#include "../algorithms/Operators.h"

#include <vector>

namespace nupack {

/******************************************************************************************/

/// Call a function until it returns true
template <class F> void while_false(F &&f) {while (!f()) {};}
/// Call a function until it returns false
template <class F> void while_true(F &&f) {while (f()) {};}

/******************************************************************************************/

static_assert(is_tuple<std::tuple<int, int &, int &&, int *>>, "ugh");

/// Tuple for_each helper
template <class Tuple, class F, std::size_t ...Is>
void for_each_tuple(Tuple &&tuple, F &&f, indices_t<Is...>) {
    (void) std::initializer_list<int>{1, (f(std::get<Is>(fw<Tuple>(tuple))), 0)...};
}

/// Tuple for_each helper
template <class Tuple, class F, std::size_t ...Is>
bool while_each_tuple(Tuple &&tuple, F &&f, indices_t<Is...>) {
    bool ok = true;
    (void) std::initializer_list<int>{1, (ok ? (ok = f(std::get<Is>(fw<Tuple>(tuple))), 0) : 0)...};
    return ok;
}

/// Tuple (compile time) for_each
template <class T, class F, NUPACK_IF(traits::is_like_tuple<decay<T>>)>
void for_each(T &&t, F &&f) {for_each_tuple(fw<T>(t), fw<F>(f), indices_in<T>());}

template <class T, class F, NUPACK_IF(traits::is_like_tuple<decay<T>>)>
bool while_each(T &&t, F &&f) {return while_each_tuple(fw<T>(t), fw<F>(f), indices_in<T>());}

/// Dynamic for_each
template <class T, class F, NUPACK_IF(!traits::is_like_tuple<decay<T>>)>
void for_each(T &&t, F &&f) {for (auto &&i : t) f(i);}

template <class T, class F, NUPACK_IF(!traits::is_like_tuple<decay<T>>)>
bool while_each(T &&t, F &&f) {for (auto &&i : t) if(!f(i)) return false; return true;}

/******************************************************************************************/

/// Call function with zipped pairs of arguments from 2 tuples
template <class T1, class T2, class F, std::size_t ...Is>
void for_each_zip_impl(T1 &&t1, T2 &&t2, F &&f, indices_t<Is...>) {
    (void) std::initializer_list<int>{1, (f(std::get<Is>(fw<T1>(t1)), std::get<Is>(fw<T2>(t2))), void(), int{})...};
}

/// Call function with zipped pairs of arguments from 2 tuples
template <class T1, class T2, class F>
void for_each_zip(T1 &&t1, T2 &&t2, F &&f) {
    static_assert(tuple_size<T1> == tuple_size<T2>, "tuples should be same size");
    for_each_zip_impl(fw<T1>(t1), fw<T2>(t2), fw<F>(f), indices_in<T1>());
}

template <std::size_t ...Is, class F>
void for_each_index(indices_t<Is...>, F &&f) {(void) std::initializer_list<int>{(f(std::integral_constant<std::size_t, Is>{}), 0)...};}

template <std::size_t N, class F>
void for_each_index(F &&f) {for_each_index(indices_up_to<N>(), fw<F>(f));}

template <std::size_t ...Is, class F>
bool while_each_index(indices_t<Is...>, F &&f) {
    bool ok = true;
    (void) std::initializer_list<bool>{true, ok && (ok = f(std::integral_constant<std::size_t, Is>{}))...};
    return ok;
}

template <std::size_t N, class F>
bool while_each_index(F &&f) {return while_each_index(indices_up_to<N>(), fw<F>(f));}

template <std::size_t ...Is>
constexpr auto indices_c(indices_t<Is...> = {}) {return std::make_tuple(std::integral_constant<std::size_t, Is>{}...);}

/******************************************************************************************/

/// Call a function for each argument
//https://isocpp.org/blog/2015/01/for-each-argument-sean-parent
template <class F, class ...Args>
void for_each_arg(F &&f, Args &&...args) {
    (void) std::initializer_list<int>{(f(args), 0)...};
}

/******************************************************************************************/

/// Call a function with all pairs (i, j) in a range where j >= i+o
template <class B, class E, class F> void for_ordered_pairs(B b, E e, int o, F &&f) {
    for (auto i = b; i < e; ++i) for (auto j = i + o; j < e; ++j) f(i, j);}

/// Call a function with all pairs (i, j) in a range where j > i
template <class B, class E, class F> void for_ordered_pairs(B b, E e, F &&f) {
    for (auto i = b; i < e; ++i) for (auto j = i + 1; j < e; ++j) f(i, j);}

/******************************************************************************************/

/// Call a function with all pairs (i, j) in a range
template <class B, class E, class F> void for_pairs(B b, E e, F &&f) {
    for (auto i = b; i != e; ++i) for (auto j = b; j != e; ++j) f(i, j);}

/******************************************************************************************/

template <class T1, class T2, class Integer, class F>
void for_chunks(T1 i, T2 j, Integer n, F const &f) {
    Integer space = Integer(j - i) / n;
    for (auto k = 0; k < n - 1; ++k)
        f(k, std::make_pair(i + k * space, i + k * space + space));
    f(n - 1, std::make_pair(i + (n - 1) * space, j));
}

/******************************************************************************************/

template <class T1, class T2, class T3, class T4, class Integer, class F>
void for_blocks(T1 i1, T2 i2, Integer n1, T3 j1, T4 j2, Integer n2, F const &f) {
    for_chunks(i1, i2, n1, [&](auto i, auto pi){
        for_chunks(j1, j2, n2, [&](auto j, auto pj){f(i, pi, j, pj);});
    });
}

/******************************************************************************************/

/// Call a function with each permutation of a sequence (from iterators)
template <class It, class F>
void for_permutations(It b, It e, F &&f) {
    std::sort(b, e); f(); if (b == e) return;
    while (std::next_permutation(b, e)) f();
}

/// Call a function with each permutation of a sequence (from view)
template <class V, class F>
void for_permutations(V v, F &&f) {for_permutations(begin_of(v), end_of(v), [&]{f(v);});}

/// Compute all necklaces of a given size from n unique elements
/// The input vector v should be of the desired size and filled with 0s in general
template <class V, class F=NoOp>
std::size_t compute_necklaces(V &&v, uint const n_elements, F &&f={}) {
    uint const size = std::size(v);
    if (!size || !n_elements) return 0;

    std::size_t n_necklaces = 1;
    f(static_cast<V const &>(v));

    auto find = [&] { // find the last element which is not n - 1
        uint i = 0;
        for (; i != size; ++i) if (v[size - 1 - i] + 1 != n_elements) break;
        return i;
    };

    for (uint i = find(); i != size; i = find()) {
        ++v[size - 1 - i];
        std::copy(std::begin(v), std::begin(v) + i, std::end(v) - i);

        if (size % (size - i) == 0) {++n_necklaces; f(static_cast<V const &>(v));}
    }
    return n_necklaces;
}

/******************************************************************************************/

/// Return the equivalent vector but rotated so that it is lexicographically lowest
/// Could be done without this concatenation strategy in the future
template <class V>
std::size_t lowest_rotational_order(V const &v) {
    if (len(v) <= 1) return 0;
    V v2{v};
    auto n = len(v2);
    cat(v2, v);
    auto const b = begin_of(v2);
    std::size_t best = 0;
    for (std::size_t i = 1u; i != n; ++i)
        if (std::lexicographical_compare(b + i, b + i + n, b + best, b + best + n))
            best = i;
    return best;
}

/******************************************************************************************/

template <class V>
V lowest_rotation(V v) {
    auto const i = lowest_rotational_order(v);
    std::rotate(v.begin(), std::next(v.begin(), i), v.end());
    return v;
}

/******************************************************************************************/

/// Return the rotational symmetry number of a sequence (1 if there is no symmetry)
template <class V>
std::size_t rotational_symmetry(V const &v) {
    auto const b = begin_of(v), e = end_of(v);
    std::size_t const n = std::distance(b, e);
    for (std::size_t i = 1; i <= n / 2; ++i) if (n % i == 0)
        if (std::equal(b+i, e, b, e-i) && std::equal(b, b+i, e-i, e)) return n / i;
    return 1;
}

/******************************************************************************************/

template <class F>
void prime_factorization(std::size_t n, F &&f) {
    std::size_t z = 2;
    while (z * z <= n) {
        if (n % z == 0) {f(std::size_t(n)); n /= z;}
        else ++z;
    }
    if (n > 1) f(n);
}

/******************************************************************************************/

/// iterator product of for_permutations
template <class B, class E, class F>
void for_permutations_across(B begin, E end, F &&f) {
    if (begin == end) {f(); return;}
    for_permutations(begin_of(*begin), end_of(*begin), [&] {
        for_permutations_across(begin + 1, end, f);
    });
}

/// for each place an element can go into a container
template <class V, class T, class F>
void for_splices(bool first, V &v, T &&t, F &&f) {
    v.insert(begin_of(v), fw<T>(t));
    auto a = begin_of(v), b = std::next(a);
    if (first || b == end_of(v)) f();
    for (; b != end_of(v); a = b, ++b) {
        swap(*a, *b);
        f();
    }
    v.pop_back();
}

/// for each partitioning of a group of elements -- inner loops
template <class It, class V, class F>
void for_partitions(bool first, std::size_t n, It const m, It const e, V &v, F &&f) {
    if (m == e) {f(v); return;}
    for (auto &i : v) for_splices(first, i, *m, [&] {
        for_partitions(first, n, std::next(m), e, v, f);
    });
    v.emplace_back();
    v.back().reserve(n); // prevent any reallocation so iterators don't go bad
    v.back().emplace_back(*m);
    for_partitions(first, n, std::next(m), e, v, f);
    v.pop_back(); // or *b as separate complex
}

/// for each partitioning of a group of elements
template <class T, class F>
void for_partitions(bool first, T const &t, F &&f) {
    std::vector<std::vector<value_type_of<T>>> v;
    v.reserve(len(t));
    for_partitions(first, 2*len(t), begin_of(t), end_of(t), v, f);
}

template <class B, class E, class V, class F>
void for_choose_any(bool first, B const b, E const e, V &v, F &&f) {
    auto c = std::next(b);
    if (c == e) {
        for_splices(first, v, *b, [&] {f(v);});
        if (first || len(v)) f(v); // not sure if worth adding argument for empty cases
    } else {
        for_splices(first, v, *b, [&] {for_choose_any(first, c, e, v, f);});
        for_choose_any(first, c, e, v, f);
    }
}

template <class T, class F>
void for_choose_any(bool first, T const &t, F &&f) {
    std::vector<value_type_of<T>> v;
    if (!len(t)) {f(v); return;}
    v.reserve(len(t));
    for_choose_any(first, begin_of(t), end_of(t), v, f);
}

/******************************************************************************************/

/// Next iterator, except .end_of()-1 returns begin_of()
template <class V, class It> It cyclic_next(V const &v, It const &it) {
    auto const up = std::next(it); return (up == end_of(v)) ? begin_of(v) : up;
}

/// Previous iterator, except .begin_of() returns .end_of()[-1]
template <class V, class It> It cyclic_prev(V const &v, It const &it) {
    return std::prev(it == begin_of(v) ? end_of(v) : it);
}

template <class V, class F> void for_circularly_adjacent(V &&v, F &&f) {
    auto x = begin_of(v);
    if (x == end_of(v)) return;
    auto y = std::next(x);
    while (y != end_of(v)) {
        f(*x, *y);
        x = std::move(y);
        y = std::next(x);
    }
    f(*x, *begin_of(v));
}

/******************************************************************************************/

NUPACK_DETECT(has_back, decltype(std::declval<T>().back()));
NUPACK_DETECT(has_front, decltype(std::declval<T>().front()));

/// First element in a container
template <class V, NUPACK_IF(!has_front<V &&>)> decltype(auto) front(V &&v) {return *begin_of(fw<V>(v));}
/// Last element in a container
template <class V, NUPACK_IF(!has_back<V &&>)> decltype(auto) back(V &&v) {return *std::prev(end_of(fw<V>(v)));}
/// First element in a container
template <class V> auto front(V &&v) -> decltype(fw<V>(v).front()) {return fw<V>(v).front();}
/// Last element in a container
template <class V> auto back(V &&v) -> decltype(fw<V>(v).back()) {return fw<V>(v).back();}

/// Element from indexing a container forwards
template <class V, class N> decltype(auto) front(V &&v, N &&n) {return *std::next(begin_of(fw<V>(v)), n);}
/// Iterator from indexing a container forwards
template <class V, class N> decltype(auto) next(V &&v, N &&n) {return std::next(begin_of(fw<V>(v)), n);}
/// Element from indexing a container backwards. 0 gives the last element.
template <class V, class N> decltype(auto) back_index(V &&v, N &&n) {return *std::prev(end_of(fw<V>(v)), n+1);}

/******************************************************************************************/

template <class ...Cs> void zip_iterators(Cs &&...cs) {
    constexpr auto n = sizeof...(Cs);
    auto ts = slice<0, n-1>(std::forward_as_tuple(cs...));
    for (auto its = map_each(ts, begin_of);
         first_of(its) != end_of(first_of(ts));
         for_each(its, increment))
        unpack(its, [&](auto ...is){
            std::get<n-1>(std::forward_as_tuple(cs...))(is...);
        });
}

/******************************************************************************************/

struct ignore {
    template <class T> constexpr ignore(T const &) {}
    template <class T> constexpr ignore const & operator=(T const &) const {return *this;}
};

template <std::size_t ...Is>
struct arg_t {
    template <class T, class ...Ts>
    T constexpr operator()(sink_type<ignore, decltype(Is)>..., T &&t, Ts const &...) const {return fw<T>(t);}
};

template <std::size_t ...Is> arg_t<Is...> make_arg_t(indices_t<Is...>);
template <std::size_t N> static constexpr auto arg_c = decltype(make_arg_t(indices_up_to<N>())){};

/******************************************************************************************/

/// Iterate through a tuple of iterators, the visitor function is the last argument; returns iterators when done
template <class T, std::size_t ...Is>
auto zip_tuple(indices_t<Is...>, T t) {
    auto &&f = std::get<sizeof...(Is)>(t);
    auto iters = std::make_tuple(begin_of(std::get<Is>(t))...);
    for (; first_of(iters) != end_of(first_of(t)); for_each(iters, increment))
        f((*std::get<Is>(iters))...);
    return iters;
}

/// Iterate through a zip of iterators, the visitor function is the last argument; returns iterators when done
template <class ...Cs>
auto zip(Cs &&...cs) {return zip_tuple(indices_up_to<sizeof...(Cs)-1>(), std::forward_as_tuple(cs...));}

template <class T, std::size_t ...Is>
std::size_t izip_tuple(indices_t<Is...>, T t) {
    std::size_t n = 0;
    auto &&f = std::get<sizeof...(Is)>(t);
    for (auto is = std::make_tuple(begin_of(std::get<Is>(t))...); first_of(is) != end_of(first_of(t)); for_each(is, increment), ++n)
        f(n, (*std::get<Is>(is))...);
    return n;
}

/// Iterate through a zip of iterators plus a prepended counter, the visitor function is the last argument
template <class ...Cs, std::enable_if_t<sizeof...(Cs) != 2, int> = 0>
auto izip(Cs &&...cs) {
    return izip_tuple(indices_up_to<sizeof...(Cs)-1>(), std::forward_as_tuple(cs...));
}

/// Iterate through a zip of iterators plus a prepended counter, the visitor function is the last argument
template <class C, class F>
auto izip(C &&c, F &&f) {
    std::size_t n = 0;
    for (auto it = begin_of(c); it != end_of(c); ++it) f(n++, *it);
    return n;
}

/******************************************************************************************/

/// Call a functor with unpacked elements, e.g. of a vector of tuples
template <class V, class F>
auto unzip(V &&v, F &&f) {for (auto &&i : v) unpack(fw<decltype(i)>(i), f);}

/******************************************************************************************/

/// Move iterator to beginning of container
template <class V>
auto move_begin(V &&v) {return std::make_move_iterator(begin_of(v));}
/// Move iterator to end of container
template <class V>
auto move_end(V &&v) {return std::make_move_iterator(end_of(v));}

/******************************************************************************************/

// given indices ts, make reordering such that order[ts...] = 0,1,2,3...
template <class ...Ts>
std::vector<std::size_t> shifted_to_front(std::size_t n, Ts const &...ts) {
    std::vector<std::size_t> ret(n);
    for (uint i = 0u; i != n; ++i) ret[i] = i;
    std::size_t i = 0;
    (void) std::initializer_list<int> {(swap(ret[ts], ret[i++]), 0)...};
    return ret;
}

/// permute a container according to a sequence of indices
template <class V, class Order>
void reorder(V &v, Order const &o)  {
    std::size_t done = 0;
    for (value_type_of<Order> s = 0, d; done != len(o); ++s) {
        for (d = o[s]; s < d; d = o[d]);
        if (d == s) {
            ++done;
            auto temp = std::move(v[s]);
            while (d = o[d], d != s) {swap(temp, v[d]); ++done;}
            v[s] = std::move(temp);
        }
    }
}

/******************************************************************************************/

/// construct a container from another directly if allowed, else construct with STL iterator pair
template <class O, class T, NUPACK_IF(can_construct<O, T &&>)>
O cast_container(T &&t) {return O{fw<T>(t)};}

template <class O, class T, NUPACK_IF(!can_construct<O, T &&>)>
O cast_container(T &&t) {return O{move_begin(t), move_end(t)};}

template <class O, class T, NUPACK_IF(!can_construct<O, T const &>)>
O cast_container(T const &t) {return O{begin_of(t), end_of(t)};}

/******************************************************************************************/

/// Fold a binary functor across the elements
template <class Op, class T> constexpr T fold(Op &&, T &&t) {return static_cast<T &&>(t);}
/// Fold a binary functor across the elements
template <class Op, class T, class U, class ...Ts>
constexpr decltype(auto) fold(Op &&op, T &&t, U &&u, Ts &&...ts) {
    return fold(op, op(fw<T>(t), fw<U>(u)), fw<Ts>(ts)...);
}

/******************************************************************************************/

}
