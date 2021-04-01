/**
 * @brief Common mathematics functions
 *
 * @file Numeric.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Operators.h"
#include "Functor.h"
#include "../iteration/Patterns.h"
#include <numeric>
#include <bitset>
#include <limits>
#include <array>
#include <cstring>

namespace nupack {

NUPACK_UNARY_FUNCTOR(round,          std::round(t));
NUPACK_UNARY_FUNCTOR(abs,            std::abs(t));
NUPACK_UNARY_FUNCTOR(exp,            std::exp(t));
NUPACK_UNARY_FUNCTOR(sqroot,         std::sqrt(t));
NUPACK_UNARY_FUNCTOR(log,            std::log(t));
NUPACK_UNARY_FUNCTOR(ceil,           std::ceil(t));
NUPACK_UNARY_FUNCTOR(floor,          std::floor(t));
NUPACK_UNARY_FUNCTOR(sign,           !std::signbit(t));
NUPACK_UNARY_FUNCTOR(is_nonnegative, !std::signbit(t));
NUPACK_UNARY_FUNCTOR(is_negative,    std::signbit(t));
NUPACK_UNARY_FUNCTOR(is_positive,    (t > static_cast<std::decay_t<decltype(t)>>(0)));
NUPACK_UNARY_FUNCTOR(is_nan,         std::isnan(t));
NUPACK_UNARY_FUNCTOR(is_finite,      std::isfinite(t));

template <class T, class U>
constexpr auto pow(T t, U u) -> decltype(std::pow(t, u)) {return std::pow(t, u);}

/******************************************************************************************/

template <class T, class=void>
struct numeric_limits : std::numeric_limits<T> {static_assert(std::is_arithmetic_v<T>);};

template <class T> constexpr auto max_log() {return (numeric_limits<T>::max_exponent * 445u) / 642u;} // is about log2
template <class T> constexpr auto max_log2() {return numeric_limits<T>::max_exponent;}
template <class T> constexpr auto max_log10() {return numeric_limits<T>::max_exponent10;}

template <class T> constexpr auto min_log() {return (numeric_limits<T>::min_exponent * 445u) / 642u;} // is about log2
template <class T> constexpr auto min_log2() {return numeric_limits<T>::min_exponent;}
template <class T> constexpr auto min_log10() {return numeric_limits<T>::min_exponent10;}

template <class T>
static constexpr auto max_limit = numeric_limits<T>::max;

template <class T>
static constexpr auto min_limit = numeric_limits<T>::min;

/******************************************************************************************/

/// square a number
struct square_t {
    template <class T, NUPACK_IF(!can_multiply<T>)>
    constexpr decltype(auto) operator()(T &&t) const {return pow(fw<T>(t), 2);}
    template <class T, NUPACK_IF(can_multiply<T>)>
    constexpr decltype(auto) operator()(T &&t) const {return t * t;}
};

static constexpr auto sq = square_t();

/// cube a number
struct cube_t {
    template <class T, NUPACK_IF(!can_multiply<T>)>
    constexpr decltype(auto) operator()(T &&t) const {return pow(fw<T>(t), 3);}
    template <class T, NUPACK_IF(can_multiply<T>)>
    constexpr decltype(auto) operator()(T &&t) const {return t * t * t;}
};

static constexpr auto cube = cube_t();

/******************************************************************************************/

/// solve quadratic formula with first coefficient=1, return both roots
template <class T> auto quadratic_solve(T const b, T const c) {
    auto r1 = (-b - sqroot(sq(b) - 4 * c)) / 2; // should be optimized away
    auto r2 = (-b + sqroot(sq(b) - 4 * c)) / 2;
    return std::make_pair(r1, r2);
}

/// solve quadratic formula, return both roots
template <class T> auto quadratic_solve(T a, T b, T c) {return quadratic_solve(b / a, c / a);}

/******************************************************************************************/

template <std::size_t N, class T, NUPACK_IF(N == 0)>
constexpr T legendre(T const &) {return 1;}

template <std::size_t N, class T, NUPACK_IF(N == 1)>
constexpr T const & legendre(T const &t) {return t;}

template <std::size_t N, class T, NUPACK_IF(N >= 2)>
constexpr auto legendre(T const &t) {return ((2 * N + 1) * t * legendre<N-1>(t) - N * legendre<N-2>(t)) / (N + 1);}

template <class T, class F>
void legendres(std::size_t min, std::size_t max, T const t, F &&f) {
    static_assert(!is_integral<T>, "Legendre function expects floating point numbers");
    T a = 1;
    T b = t;
    for (std::size_t n = 2; n <= max; ++n) {
        a = ((2 * n + 1) * t * b - n * a) / (n + 1);
        if (min <= n) f(n, a);
        swap(a, b);
    }
}

template <class T, class F>
void legendres(std::size_t max, T const t, F &&f) {
    f(0, T{1});
    f(1, t);
    legendres(2, max, t, fw<F>(f));
}

template <class T, class F>
void powers(T const t, int p, F &&f) {for (T e = static_cast<T>(1u); p--; e *= t) f(e);}

template <class T, class F>
void powers(T const t, int p, int q, F &&f) {for (T e = pow(t, p); p != q; ++p, e *= t) f(e);}

template <class T, NUPACK_IF(is_integral<T>)>
constexpr T pow2(T t) {return static_cast<T>(1u) << t;}

/******************************************************************************************/

/// apply unary functor n and fill an iterator range with the successive results
template <class It, class T, class Unary=Identity>
void iota(It b, It e, T &&t, Unary &&u) {
    if (e == b) return;
    *b = fw<T>(t);
    for (auto i = std::next(b); i != e; ++i) *i = u(*std::prev(i));
}

/// apply unary functor n and fill a container with the successive results
template <class V, class T, class Unary=Identity>
void iota(V &v, T &&t, Unary &&u={}) {iota(begin_of(v), end_of(v), fw<T>(t), fw<Unary>(u));}

/******************************************************************************************/

/// return hi if input > hi, lo if input < lo, else input
template <class T, class Compare=less_t>
constexpr T const & clamp(T const &v, T const &lo, T const &hi, Compare comp=less) {
    return comp(v, hi) ? std::max(v, lo, comp) : std::min(v, hi, comp);
}

/******************************************************************************************/

/// R(i,j) = M1(i,j) * M2(j,i)
template <class M1, class M2> auto times_transpose(M1 const &m1, M2 const &m2) {return (m1 % m2.t()).eval();}

/// R(i,j) = M(i,j) * M(j,i)
template <class M> auto times_transpose(M const &m) {return (m % m.t()).eval();}

/******************************************************************************************/

/// Given integers x and y, return quotient (x / y) and remainder (x % y) in a std::pair
template <bool Signed=true, class X, class Y> auto div(X x, Y y) {
    using Z = signed_type_of<std::common_type_t<X, Y>>; // std::div wants signed types
    using R = if_t<Signed, Z, unsigned_type_of<Z>>; // user can force return of an unsigned type though
    auto z = std::div(static_cast<Z>(x), static_cast<Z>(y));
    return std::make_pair(static_cast<R>(z.quot), static_cast<R>(z.rem));
}

/// Return x/y if x is divisible by y, or else x/y+1
template <bool Signed=true, class X, class Y> auto ceil_quotient(X x, Y y) {
    auto p = div(x, y);
    return p.first + bool(p.second);
}

/// Subtract "b" from "a" unless "a" < "b". Return whether subtraction did not take place
template <class A, class B>
bool minus_if(A &a, B const &b) {if (a < b) return true; else a -= b; return false;}

/// Same as minus_if, but divide a by b when returning a true value
template <class A, class B>
bool minus_divide_if(A &a, B const &b) {
    if (a < b) {a /= b; return true;} else {a -= b; return false;}
}

/******************************************************************************************/

/// Returns the next power of 2, in a zero-index sense (e.g. 3, 7, 15, 31, 63...)
template <class T>
T next_power_of_two(T s) {T ret = 1; do {ret <<= 1;} while (s >>= 1); return ret - 1;}

/******************************************************************************************/

/// If a value is <= 0, return its smallest possible positive value
template <class T> T min_floor(T t) {return t > T() ? t : std::numeric_limits<T>::min();}
/// If a value is <= 0, return 0
template <class T> T zero_floor(T t) {return t > T() ? t : T();}

/******************************************************************************************/

/// weighted average for two elements
template <class T> T weight_avg(T const &t1, T const &t2, double f1, double f2) {
    return (t1 * f1 + t2 * f2) / (f1 + f2);
}
/// average of any number of numbers
template <class ...Ts> auto avg(Ts const &... ts) {return fold(plus, ts...) / sizeof...(Ts);}

/******************************************************************************************/

template <class V> auto hamming_distance(V const &v1, V const &v2) {
    return std::inner_product(begin_of(v1), end_of(v1), begin_of(v2), size_type_of<V>(), plus, not_equals);
}

/******************************************************************************************/

// Round a floating point number to "n" significant binary digits
template <class T> T sig_round(T t, int n=sizeof(T) * 4) {
    int exponent, junk;
    auto mantissa = std::frexp(t, &exponent);
    return std::ldexp(std::frexp(std::round(std::ldexp(mantissa, n)), &junk), exponent);
}

/******************************************************************************************/

template <class T, NUPACK_IF(is_floating_point<T>)>
auto binary_exponent(T t) {int r; std::frexp(t, &r); return r;}

/******************************************************************************************/

template <class V, class T> V linspace(T n) {V ret(n); iota(ret, T(*zero), plus_one); return ret;}
template <class V, class T> V linspace(T b, T e) {V ret(e - b); iota(ret, b, plus_one); return ret;}

/******************************************************************************************/

/// Return previous element or zero if there is none
template <class V, class It, class F=Identity>
auto prev_or_zero(V const &v, It it, F const &f=Identity()) -> decltype(f(it[0])) {
    return it == begin_of(v) ? decltype(f(it[0]))(*zero) : f(it[-1]);
}

/******************************************************************************************/

/// Unpack a partial sum array to get a marginal value at a position
template <class V, class It, class F=Identity>
auto element_from_sum(V const &v, It it, F const &f=Identity()) -> decltype(f(it[0])) {
    return it == begin_of(v) ? f(it[0]) : f(it[0]) - f(it[-1]);
}

/******************************************************************************************/

template <class U=void, class T, class ...Ts>
constexpr std::array<nonvoid<U, decay<T>>, 1+sizeof...(Ts)> make_array(T &&t, Ts &&...ts) {return {fw<T>(t), fw<Ts>(ts)...};}

/******************************************************************************************/

template <class T>
auto & as_bitset(T &t) {return reinterpret_cast<std::bitset<sizeof(T) * CHAR_BIT> &>(t);}

template <class T>
auto const & as_bitset(T const &t) {return reinterpret_cast<std::bitset<sizeof(T) * CHAR_BIT> const &>(t);}

template <class T>
auto as_bitset(T &&t) {return reinterpret_cast<std::bitset<sizeof(T) * CHAR_BIT> const &>(t);}

/******************************************************************************************/

template <class T>
constexpr bool bit_at(T const &t, std::size_t i) {return reinterpret_cast<std::bitset<sizeof(T) * CHAR_BIT> const &>(t)[i];}

/******************************************************************************************/

template <class ...Ts> auto log_sum_exp(Ts const &... ts) {
    auto t_max = fold(max, ts...);
    if (std::isinf(t_max)) return t_max;
    return t_max + log(fold(plus, exp(ts - t_max)...));
}

/******************************************************************************************/

template <class T>
void zero_memory(T *b, T const *e) {
    std::memset(static_cast<void *>(b), 0, reinterpret_cast<char const *>(e) - reinterpret_cast<char const *>(b));
}

template <class T>
void contiguous_fill(T *b, T *e, T const t) {
    if (std::is_pod_v<T> && reinterpret_cast<uint_of_size<sizeof(T)> const &>(t) == 0u) zero_memory(b, e);
    else std::fill(b, e, t);
}

/******************************************************************************************/

template <class T, class U>
constexpr auto unsigned_minus(T const &t, U const &u) {return std::max<T>(t, u) - u;}

}
