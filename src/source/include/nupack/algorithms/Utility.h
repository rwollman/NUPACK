/**
 * @brief Contains common algorithms used in the program, many relying on STL (std::algorithm)
 *
 * @file Utility.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include <algorithm>
#include <functional>
#include <limits>
#include <chrono>
#include <mutex>

#include "Operators.h"
#include "Functor.h"
#include "Numeric.h"

#include <boost/config.hpp>

namespace nupack {

/// tell the compiler something is likely; override with template bool=false to get unlikely
template <bool B=true, class T, NUPACK_IF(B)>
constexpr bool likely(T &&t) {bool b(fw<T>(t)); return BOOST_LIKELY(b);}
/// tell the compiler something is likely; override with template bool=false to get unlikely
template <bool B=true, class T, NUPACK_IF(!B)>
constexpr bool likely(T &&t) {bool b(fw<T>(t)); return BOOST_UNLIKELY(b);}
template <class T> constexpr bool unlikely(T &&t) {bool b(fw<T>(t)); return BOOST_UNLIKELY(b);}

/// if the template bool is true, tell compiler that it's likely, else tell the compiler nothing
template <bool B, class T, NUPACK_IF(B)> constexpr bool likely_if(T &&t) {return likely(fw<T>(t));}
template <bool B, class T, NUPACK_IF(!B)> constexpr bool likely_if(T &&t) {return bool(fw<T>(t));}

/// if the template bool is true, tell compiler that it's unlikely, else tell the compiler nothing
template <bool B, class T, NUPACK_IF(B)> constexpr bool unlikely_if(T &&t) {return unlikely(fw<T>(t));}
template <bool B, class T, NUPACK_IF(!B)> constexpr bool unlikely_if(T &&t) {return bool(fw<T>(t));}

template <class F> decltype(auto) switch_c(bool c, F &&f) {return c ? fw<F>(f)(True()) : fw<F>(f)(False());}

/******************************************************************************************/

inline void clobber_memory() {asm volatile ("" ::: "memory");}

/******************************************************************************************/

template <class T, class U, class Comp=less_t> constexpr auto ordered_pair(T &&t, U &&u, Comp c={}) {
    using R = std::common_type_t<T, U>;
    return c(t, u) ? std::pair<R, R>(fw<T>(t), fw<U>(u)) : std::pair<R, R>(fw<U>(u), fw<T>(t));
}

/******************************************************************************************/

/// return the time it takes to run a function n times
template <class F, class ...Ts>
double time_it(std::size_t n, F &&f, Ts &&...ts) {
    auto const t0 = std::chrono::high_resolution_clock::now();
    if (n) f(ts...);
    for (std::size_t i = 1; i < n; ++i) {clobber_memory(); f(ts...);}
    return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count() / n;
}

template <class F> auto time_it(F &&f) {return time_it(1, fw<F>(f));}

/******************************************************************************************/

}
