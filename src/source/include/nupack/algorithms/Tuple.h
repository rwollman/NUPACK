/**
 * @brief Algorithms to do with containers that are like std::tuple, indcluding compile-time indices
 *
 * @file Tuple.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Traits.h"
#include "Extents.h"

/******************************************************************************************/

namespace nupack {

/******************************************************************************************/

template <std::size_t N>
constexpr auto indices_up_to(size_constant<N> = {}) {return std::make_index_sequence<N>();}

template <int N, std::size_t ...Is>
constexpr indices_t<(Is+N)...> add_to_indices(indices_t<Is...>) {return {};}

template <std::size_t M, std::size_t N>
constexpr auto indices_between(size_constant<N> = {}) {return add_to_indices<M>(indices_up_to<N-M>());}

template <class N>
constexpr auto indices_up_to() {return std::make_index_sequence<N::value>();}

template <class T>
constexpr auto indices_in() {return indices_up_to<std::tuple_size<decay<T>>::value>();}

template <class T> constexpr auto indices_in(T const &) {return indices_in<T>();}

/******************************************************************************************/

namespace detail {
    template <std::size_t B, class T, std::size_t ...Is>
    constexpr auto slice(indices_t<Is...>, T &&t) {return std::forward_as_tuple(at<B + Is>(fw<T>(t))...);}

    template <class T, std::size_t ...Is>
    auto as_tie(indices_t<Is...>, T &&t) {return std::tie(std::get<Is>(fw<T>(t))...);}

    /// Unpack tuple to call a function
    template <class F, class T, std::size_t... I>
    decltype(auto) unpack(indices_t<I...>, F &&f, T &&t) {return f(at<I>(t)...);}

    template <class T, class F, std::size_t ...I>
    auto map_each(indices_t<I...>, T &&t, F &&f) {return std::make_tuple(f(at<I>(t))...);}

    template <class T, class M, class R, std::size_t ...I>
    auto map_reduce_each(indices_t<I...>, T &&t, M &&m, R &&r) {return fw<R>(r)(m(at<I>(t))...);}

    template <class T, std::size_t ...I, class ...Ts, NUPACK_IF(!any_of_c<(is_same<decltype(at<I>(std::declval<T &&>())(std::declval<Ts &&>()...)), void>) ... >)>
    auto call_each(indices_t<I...>, T &&t, Ts &&...ts) {return std::make_tuple(at<I>(fw<T>(t))(ts...)...);}

    template <class T, std::size_t ...I, class ...Ts, NUPACK_IF(any_of_c<(is_same<decltype(at<I>(std::declval<T &&>())(std::declval<Ts &&>()...)), void>) ... >)>
    void call_each(indices_t<I...>, T &&t, Ts &&...ts) {NUPACK_UNPACK(at<I>(fw<T>(t))(fw<Ts>(ts)...));}

    template <class Out, class T, std::size_t ...Is>
    Out unpack_as(indices_t<Is...>, T &&t) {return Out{static_cast<value_type_of<Out>>(at<Is>(fw<T>(t)))...};}

    template <class T, std::size_t ...Is>
    auto filled_tuple(indices_t<Is...>, T &&t) {return std::tuple<sink_type<T, decltype(Is)>...>{sink<decltype(Is)>(t)...};}

    template <class T, std::size_t ...Is>
    auto move_each(indices_t<Is...>, T &&t) {return std::forward_as_tuple(std::move(at<Is>(t))...);}

    template <class T, std::size_t ...Is>
    auto decay_each(indices_t<Is...>, T &&t) {return std::make_tuple(fw<decltype(at<Is>(t))>(at<Is>(t))...);}
}

/******************************************************************************************/

template <class ...Ts>
auto move_as_tuple(Ts &&...ts) {return std::make_tuple(std::move(ts)...);}

template <class T>
auto move_each(T &&t) {return detail::move_each(indices_in<T>(), fw<T>(t));}

template <class T>
auto decay_each(T &&t) {return detail::decay_each(indices_in<T>(), fw<T>(t));}

/******************************************************************************************/

template <class T, class F>
auto map_each(T &&t, F &&f) {return detail::map_each(indices_in<T>(), fw<T>(t), fw<F>(f));}

template <class T, class M, class R>
auto map_reduce_each(T &&t, M &&m, R &&r) {return detail::map_reduce_each(indices_in<T>(), fw<T>(t), fw<M>(m), fw<R>(r));}

template <class T, class ...Ts>
auto call_each(T &&t, Ts &&...ts) {return detail::call_each(indices_in<T>(), fw<T>(t), fw<Ts>(ts)...);}

/******************************************************************************************/

/// Return tuple of references in the index range [B, E)
template <std::size_t B, std::size_t E, class T>
constexpr auto slice(T &&t) {return detail::slice<B>(indices_up_to<E-B>(), fw<T>(t));}

/// Return tuple of references to a set of indices
template <std::size_t ...Is, class T>
constexpr auto tie_each(T &&t) {return std::forward_as_tuple(at<Is>(fw<T>(t))...);}

/******************************************************************************************/

template <class F, class T>
auto unpack(T &&t, F &&f) {return detail::unpack(indices_in<T>(), fw<F>(f), fw<T>(t));}

/******************************************************************************************/

// Make a tuple of lvalue references to an original tuple
template <class T, NUPACK_IF(is_like_tuple<no_qual<T>>)>
auto as_tie(T &&t) {return detail::as_tie(indices_in<T>(), fw<T>(t));}

template <class T, NUPACK_IF(!is_like_tuple<no_qual<T>>)>
auto as_tie(T &&t) {return std::tie(t);}

/******************************************************************************************/

// Extend a std::tie
template <class T, class ...Ts, NUPACK_IF(is_tuple<decay<T>>)>
auto extend_tie(T &&t, Ts &&...ts) {return std::tuple_cat(t, std::tie(ts...));}

// Extend a std::make_tuple
template <class T, class ...Ts, NUPACK_IF(is_tuple<decay<T>>)>
auto extend_tuple(T &&t, Ts &&...ts) {return std::tuple_cat(t, std::make_tuple(fw<Ts>(ts)...));}

/******************************************************************************************/

/// variable template for tuple size
template <class T>
static constexpr auto tuple_size = std::tuple_size<decay<T>>::value;
/// Type of element at a given position in a tuple
template <class T, std::size_t I> using tuple_type = std::tuple_element_t<I, T>;
/// length overload for compile-time size things
template <class T> struct extents::length<T, void_if<(is_tuple<T>) || (is_pair<T>)>> {
    constexpr auto operator()(T const &) const {return tuple_size<T>;}
};

/******************************************************************************************/

template <class Out, class T>
Out unpack_as(T &&t) {return detail::unpack_as<Out>(indices_in<T>(), fw<T>(t));}

template <std::size_t N, class T>
auto filled_tuple(T &&t) {return detail::filled_tuple(indices_in<T>(), fw<T>(t));}

/******************************************************************************************/

}
