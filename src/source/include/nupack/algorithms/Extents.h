/**
 * @brief begin_of, end_of, and len operations - offers a bit more generality than std:: versions
 *
 * @file Extents.h
 * @author Mark Fornace
 * @date 2018-05-18
 */
#pragma once
#include "Constants.h"
#include <iterator>

/******************************************************************************************/

namespace nupack_adl {
    template <class T>
    auto begin_of(T &&t) -> decltype(begin(std::forward<T>(t))) {return begin(std::forward<T>(t));}
    template <class T>
    auto end_of(T &&t) -> decltype(end(std::forward<T>(t))) {return end(std::forward<T>(t));}
}

namespace nupack { namespace extents {
    template <class T, class=void> struct impl : False {};

    template <class T, class=void> struct length;

    template <class T> struct length<T, void_if<has_size<T>>> {
        constexpr auto operator()(T const &t) const {return t.size();}
    };

    template <class T> struct length<T, void_if<std::is_array<T>::value>> {
        constexpr auto operator()(T const &t) const {return std::extent<T>::value;}
    };
}

/******************************************************************************************/

template <class T, NUPACK_IF(extents::impl<no_qual<T>>::value)>
constexpr auto begin_of_f(T &&t) -> decltype(extents::impl<no_qual<T>>::begin(fw<T>(t))) {return extents::impl<no_qual<T>>::begin(fw<T>(t));}

template <class T, NUPACK_IF(extents::impl<no_qual<T>>::value)>
constexpr auto end_of_f(T &&t) -> decltype(extents::impl<no_qual<T>>::end(fw<T>(t))) {return extents::impl<no_qual<T>>::end(fw<T>(t));}

/******************************************************************************************/

NUPACK_DETECT(has_std_begin, decltype(std::begin(std::declval<T>())));
NUPACK_DETECT(has_std_end, decltype(std::end(std::declval<T>())));

template <class T, NUPACK_IF(has_std_begin<no_qual<T>>)>
decltype(auto) begin_of_f(T &&t) {return std::begin(fw<T>(t));}
template <class T, NUPACK_IF(has_std_end<no_qual<T>>)>
decltype(auto) end_of_f(T &&t) {return std::end(fw<T>(t));}

/******************************************************************************************/

template <class T, NUPACK_IF(!has_std_begin<no_qual<T>> && !extents::impl<no_qual<T>>::value)>
auto begin_of_f(T &&t) -> decltype(nupack_adl::begin_of(fw<T>(t))) {return nupack_adl::begin_of(fw<T>(t));}

template <class T, NUPACK_IF(!has_std_end<no_qual<T>> && !extents::impl<no_qual<T>>::value)>
auto end_of_f(T &&t) -> decltype(nupack_adl::end_of(fw<T>(t))) {return nupack_adl::end_of(fw<T>(t));}

/******************************************************************************************/

NUPACK_UNARY_FUNCTOR(begin_of, begin_of_f(fw<T>(t)));
NUPACK_UNARY_FUNCTOR(end_of, end_of_f(fw<T>(t)));
NUPACK_UNARY_FUNCTOR(len, extents::length<no_qual<T>>()(fw<T>(t)));

/******************************************************************************************/

NUPACK_DETECT(has_len, decltype(extents::length<no_qual<T>>()(declval<T>())));

NUPACK_DETECT(is_random_access, decltype(begin_of(declref<T>()) + 2));

NUPACK_DETECT(is_iterable, decltype(
    begin_of(declref<T>()) != end_of(declref<T>()),   // begin/end and operator !=
    ++declref<decltype(begin_of(declref<T>()))>(), // operator ++
    *begin_of(declref<T>())                        // dereference operator
));

}
