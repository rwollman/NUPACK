/**
 * @brief Useful macros, mostly for defining traits and simple functors
 *
 * @file Macro.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

/******************************************************************************************/

/// Define a type alias which is T1 if T1 is well-formed, else T2 if T2 is well-formed, else SFINAE error
#define NUPACK_TEMPLATE_FALLBACK(NAME, T1, T2) \
namespace detail { \
    template <class T, class=void> struct NAME##_t_2 {}; \
    template <class T> struct NAME##_t_2<T, void_t<T2>> {using type = T2;}; \
    template <class T, class=void> struct NAME##_t : NAME##_t_2<T> {}; \
    template <class T> struct NAME##_t<T, void_t<T1>> {using type = T1;}; \
} \
template <class T> using NAME = typename detail::NAME##_t<T>::type;

/******************************************************************************************/

/// This part is in the signature template: template <{class ...Ts}>
#define NUPACK_TEMPLATE_SIGNATURE_ELEMENT(r, data, i, elem) BOOST_PP_COMMA_IF(i) elem BOOST_PP_CAT(T, i)
/// This part is in the expansion: struct<{Ts...}>
#define NUPACK_TEMPLATE_SIGNATURE_NAME(r, data, i, elem) BOOST_PP_COMMA_IF(i) BOOST_PP_CAT(T, i)
/// Extend a trait such that a template of the given arguments yields true. The last argument is made variadic
#define NUPACK_EXTEND_TEMPLATES(ELLIPSIS, NAME, TYPE, ...) namespace detail { \
    template <BOOST_PP_SEQ_FOR_EACH_I(NUPACK_TEMPLATE_SIGNATURE_ELEMENT, 0, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__  ELLIPSIS))> \
    struct NAME##_t<TYPE<BOOST_PP_SEQ_FOR_EACH_I(NUPACK_TEMPLATE_SIGNATURE_NAME, 0, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) ELLIPSIS>> : std::true_type {}; \
}

#define NUPACK_DEFINE_TEMPLATES(ELLIPSIS, NAME, TYPE, ...) \
namespace detail {template <class T> struct NAME##_t : std::false_type {};} \
NUPACK_EXTEND_TEMPLATES(ELLIPSIS, NAME, TYPE, __VA_ARGS__) \
namespace traits {template <class T> static constexpr bool NAME = detail::NAME##_t<T>::value;}

/// Define a trait such that a template of the given arguments yields true. The last argument is not made variadic
#define NUPACK_DEFINE_TEMPLATE(NAME, TYPE, ...) NUPACK_DEFINE_TEMPLATES(, NAME, TYPE, __VA_ARGS__)
#define NUPACK_EXTEND_TEMPLATE(NAME, TYPE, ...) NUPACK_EXTEND_TEMPLATES(, NAME, TYPE, __VA_ARGS__)

/// Define a trait such that a template of the given arguments yields true. The last argument is made variadic
#define NUPACK_EXTEND_VARIADIC(NAME, TYPE, ...) NUPACK_EXTEND_TEMPLATES(..., NAME, TYPE, __VA_ARGS__)
#define NUPACK_DEFINE_VARIADIC(NAME, TYPE, ...) NUPACK_DEFINE_TEMPLATES(..., NAME, TYPE, __VA_ARGS__)

/******************************************************************************************/

/// Define a trait such that only the given type yields true
#define NUPACK_DEFINE_TYPE(NAME, TYPE) \
namespace detail { \
    template <class T> struct NAME##_t : std::false_type {}; \
    template <> struct NAME##_t<TYPE> : std::true_type {}; \
} \
namespace traits {template <class T> static constexpr bool NAME = detail::NAME##_t<T>::value;}

/// Add a specific type to a given trait
#define NUPACK_EXTEND_TYPE(NAME, TYPE) \
namespace detail {template <> struct NAME##_t<TYPE> : std::true_type {};}

/// Detect if a given type is well-formed
#define NUPACK_DETECT(NAME, expr) \
namespace detail { \
    template <class T, class=void> struct NAME##_t : std::false_type {}; \
    template <class T> struct NAME##_t<T, ::nupack::void_t<expr>> : std::true_type {}; \
} \
namespace traits {template <class T> static constexpr bool NAME = detail::NAME##_t<T>::value;}

/// Detect if a given type is well-formed
#define NUPACK_DETECT_2(NAME, expr) \
namespace detail { \
    template <class T, class U, class=void> struct NAME##_t : std::false_type {}; \
    template <class T, class U> struct NAME##_t<T, U, ::nupack::void_t<expr>> : std::true_type {}; \
} \
namespace traits {template <class T, class U> static constexpr bool NAME = detail::NAME##_t<T, U>::value;}

/// Detect if a given type is well-formed; no namespace entered
#define NUPACK_DETECT_INLINE(NAME, expr) \
template <class T, class=void> struct NAME##_t : std::false_type {}; \
template <class T> struct NAME##_t<T, ::nupack::void_t<expr>> : std::true_type {}; \
template <class T> static constexpr bool NAME = NAME##_t<T>::value;

/// Perform a detection of a given infix operator which returns a boolean
#define NUPACK_DETECT_BOOL_OP(NAME, NAME2, op) \
namespace detail { \
    template <class T, class U, class=void> struct NAME##_detect : std::false_type {}; \
    template <class T, class U> struct NAME##_detect <T, U, ::nupack::void_if<(std::is_convertible<decltype(std::declval<T>() op std::declval<U>()), bool>::value)>> : std::true_type {}; \
    template <class T, class U, class=void> struct NAME##_t : std::false_type {}; \
    template <class T, class U> struct NAME##_t<T, U, ::nupack::void_t<decltype(NAME##_detect<T, U>::value)>> : std::integral_constant<bool, NAME##_detect<T, U>::value> {}; \
} \
namespace traits { \
    template <class T> static constexpr bool NAME = detail::NAME##_t<T, T>::value; \
    template <class T, class U> static constexpr bool NAME2 = detail::NAME##_t<T, U>::value; \
}

/// Make a functor forwarding two arguments T &&t and U &&u to a given expression
#define NUPACK_BINARY_FUNCTOR(NAME, Op) \
struct NAME##_t {template <class T, class U> constexpr auto operator()(T &&t, U &&u) const -> decltype(Op) {return Op;}}; \
static constexpr auto NAME = NAME##_t()

/// Make a functor forwarding one argument T &&t to a given expression
#define NUPACK_UNARY_FUNCTOR(NAME, Op) \
struct NAME##_t {template <class T> constexpr auto operator()(T &&t) const -> decltype(Op) {return Op;}}; \
static constexpr auto NAME = NAME##_t()

/// Make a functor (lambda) forwarding all arguments to a given expression
#define NUPACK_FUNCTOR(NAME) [&](auto &&...ts) -> decltype(NAME(std::forward<decltype(ts)>(ts)...)) {return NAME(std::forward<decltype(ts)>(ts)...);}
/// Evaluate an expression in order on a pack expression
#define NUPACK_UNPACK(...) {(void) std::initializer_list<bool>{(__VA_ARGS__, false)...};}
/// Evaluate an expression in order on a pack expression, shortcircuiting if it returns false; nupack_while_ok holds if shortcircuiting did not occur
#define NUPACK_WHILE(...) bool nupack_while_ok = true; (void) std::initializer_list<bool>{(nupack_while_ok = nupack_while_ok && (__VA_ARGS__))...};

#define NUPACK_TAG(NAME) struct NAME##_t {}; static constexpr auto const NAME = NAME##_t{};

#define NUPACK_STRINGIFY(X) #X

/******************************************************************************************/

#include "TypeSupport.h"
