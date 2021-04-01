/** \file Array.h
 * @brief Defines collective array operations for Array (aka std::array)
 */
#pragma once
#include "../algorithms/Utility.h"
#include "../reflect/Print.h"
#include "../reflect/Memory.h"
#include "../reflect/Hash.h"
#include <array>

/******************************************************************************************/

namespace nupack {

NUPACK_DEFINE_TEMPLATE(is_array, std::array, class, std::size_t);

template <class T>
struct hash<T, void_if<(has_hash<value_type_of<T>>) && (is_array<T>)>> : RangeHash<T> {};

/******************************************************************************************/

namespace detail {
    template <class T, std::size_t d1, std::size_t d2, std::size_t ...ds>
    struct multi_array {using type = std::array<typename multi_array<T, d2, ds...>::type, d1>;};

    template <class T, std::size_t d1, std::size_t d2>
    struct multi_array<T, d1, d2> {using type = std::array<std::array<T, d2>, d1>;};

    template <class T, class=void> struct root_value_type {using type = T;};

    template <class T> struct root_value_type<T, void_t<value_type_of<T>>> {
        using type = type_in<root_value_type<value_type_of<T>>>;
    };
}

/// Alias for multidimensional std::array
template <class T, std::size_t d1, std::size_t d2, std::size_t... ds>
using multi_array = type_in<detail::multi_array<T, d1, d2, ds...>>;

/// Alias for multidimensional std::array raw element type
template <class A> using root_value_type = type_in<detail::root_value_type<A>>;

/******************************************************************************************/

namespace detail {
    template <class A, NUPACK_IF(is_array<A> && !is_array<value_type_of<A>>)>
    constexpr auto dimensions_as_tuple() {return std::make_tuple(std::tuple_size<A>::value);}

    template <class A, NUPACK_IF(is_array<A> && is_array<value_type_of<A>>)>
    constexpr auto dimensions_as_tuple() {
        return std::tuple_cat(std::make_tuple(std::tuple_size<A>::value),
                              dimensions_as_tuple<value_type_of<A>>());
    }

    template <class T, std::size_t ...Is>
    constexpr auto dimensions_as_array(T const &t, indices_t<Is...>) {
        return std::array<std::size_t, sizeof...(Is)>{std::get<Is>(t)...};
    }
}

/// Get dimensions of a multidimensional std::array
template <class A, NUPACK_IF(is_array<A>)> constexpr auto shape() {
    auto t = detail::dimensions_as_tuple<A>();
    return detail::dimensions_as_array(t, indices_in(t));
}

/******************************************************************************************/

/// Print array like a list
template <class T> struct io::PrintAsContainer<T, void_if<is_array<T>>> : PrintAsList {};

/******************************************************************************************/

/// Sum memory of each element in the array
template <class T> struct memory::impl<T, void_if<is_array<T>>> {
    constexpr std::size_t operator()(T const &t) const {return sum(t, measure);}
    void erase(T &t) const {for (auto &i : t) release(i);}
};

/******************************************************************************************/

template <class T, std::size_t N> struct detail::is_like_tuple_t<std::array<T, N>> : True {};

/******************************************************************************************/

/// Get a flattened view of a multidimensional array
template <class A, NUPACK_IF(is_array<A>)>
auto flat_view(A &a) {
    using T = root_value_type<decay<A>>;
    return ptr_view(reinterpret_cast<T *>(&a), sizeof(A) / sizeof(T));
}

/******************************************************************************************/

}
