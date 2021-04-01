/**
 * @brief Defines common NUPACK exceptions and exception-throwing macros
 *
 * @file Error.h
 * @author Mark Fornace
 * @date 2018-05-31
 * @todo Redo macros to use fmt library or similar
 */
#pragma once
#include "Runtime.h"
#include "Config.h"
#include "../algorithms/Macro.h"

#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>

/******************************************************************************************/

/// Runtime error exception, containing message and string/value pairs
#define NUPACK_ERROR(msg, ...) ({ \
    std::stringstream NUPACK_BUFFER; \
    NUPACK_BUFFER << std::setprecision(16) << std::boolalpha << ::nupack::message_string(NUPACK_FILE, __LINE__, msg); \
    BOOST_PP_SEQ_FOR_EACH(NUPACK_APPEND, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    throw ::nupack::Error(NUPACK_BUFFER.str()); \
})

/******************************************************************************************/

/// Programming error exception, containing message and string/value pairs
#define NUPACK_BUG(msg, ...) ({ \
    std::stringstream NUPACK_BUFFER; NUPACK_BUFFER << std::boolalpha; \
    NUPACK_BUFFER << ::nupack::message_string(NUPACK_FILE, __LINE__, msg); \
    BOOST_PP_SEQ_FOR_EACH(NUPACK_APPEND, , BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    throw ::nupack::Bug(NUPACK_BUFFER.str()); \
})
/// Assert macro printing key-value pairs in ...
#define NUPACK_ASSERT(x, ...) if (BOOST_UNLIKELY(!(x))) {NUPACK_ERROR("Assertion failure", x, __VA_ARGS__);}
/// Assert with augmented printing for comparisons
#define NUPACK_REQUIRE(lhs, op, rhs, ...) NUPACK_ASSERT(lhs op rhs, lhs, rhs, __VA_ARGS__)

/// Debug versions
#define NUPACK_DASSERT(x, ...) if (Debug && BOOST_UNLIKELY(!(x))) {NUPACK_ERROR("Assertion failure", x, __VA_ARGS__);}
#define NUPACK_DREQUIRE(lhs, op, rhs, ...) NUPACK_DASSERT(lhs op rhs, lhs, rhs, __VA_ARGS__)

/******************************************************************************************/

#define NUPACK_ALL_EQUAL(msg, ...) NUPACK_ASSERT(::nupack::equal_arguments(__VA_ARGS__), msg, __VA_ARGS__)

namespace nupack {

template <class T, class ...Ts>
constexpr bool equal_arguments(T const &t, Ts const &...ts) {return ((t == ts) && ...);}

string message_string(string fn, int line, string msg);
/******************************************************************************************/

template <class ...Ts>
static string join_message(Ts const &...ts) {
    std::stringstream ss;
    (void) std::initializer_list<int>{(ss << ts, 0)...};
    string s = ss.str();
    std::replace(s.begin(), s.end(), '\0', '0');
    return s;
}

/// User caused error
struct Error: public std::runtime_error {
    Error(std::string);

    template <class ...Ts>
    Error(Ts const &...ts): Error(join_message("NUPACK error: ", ts...)) {};
};

/// Developer caused error
struct Bug: public std::runtime_error {
    Bug(string const &message): std::runtime_error(join_message("NUPACK bug", message, '\n', get_backtrace())){};
};

/// Compile-time deduced error
struct CompileError: public std::exception {
    CompileError(string const &message) {};
};

/******************************************************************************************/

NUPACK_DETECT_2(has_at, decltype(declval<T>().at(declval<U>())));

template <class V, class N, bool B=true, NUPACK_IF(B && !DebugBounds)>
auto at(V &&v, N &&n) -> decltype(fw<V>(v)[fw<N>(n)]) {return fw<V>(v)[fw<N>(n)];}

template <class V, class N, NUPACK_IF(traits::has_at<V &&, N &&> && DebugBounds)>
decltype(auto) at(V &&v, N &&n) {
    try {return fw<V>(v).at(fw<N>(n));}
    catch(...) {throw Error("Element access error with container of type ", type_name<V>, " and index of type ", type_name<N>);}
}

template <class V, class N, NUPACK_IF(!traits::has_at<V &&, N &&> && DebugBounds)>
decltype(auto) at(V &&v, N &&n) {
    if (n >= len(v)) throw Error("Element access error with container of type ", type_name<V>, " and index of type ", type_name<N>);
    return fw<V>(v)[fw<N>(n)];
}

/******************************************************************************************/

}

