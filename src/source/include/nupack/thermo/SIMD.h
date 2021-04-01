/**
 * @brief SIMD functions and Boost.SIMD wrapper for dynamic programs
 *
 * @file SIMD.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../iteration/Range.h"
#include "../common/Error.h"
#include "../reflect/Print.h"
#include <boost/simd/pack.hpp>
#include <boost/simd/function/saturated.hpp> // needed to fix compiling on EC2 for some reason
#include <boost/simd/function/ldexp.hpp>
#include <boost/simd/function/ilogb.hpp>
#include <boost/simd/function/ifrexp.hpp>
#include <boost/simd/function/sum.hpp>
#include <boost/simd/function/multiplies.hpp>
#include <boost/simd/function/plus.hpp>
#include <boost/simd/function/is_eqz.hpp>
#include <boost/simd/function/if_zero_else.hpp>
#include <boost/simd/function/minimum.hpp>
#include <boost/simd/function/maximum.hpp>
#include <boost/simd/function/min.hpp>
#include <boost/simd/function/rec.hpp>
#include <boost/simd/function/unary_minus.hpp>
#include <boost/simd/function/store.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/reverse.hpp>
#include <boost/simd/function/is_nan.hpp>
#include <boost/simd/function/any.hpp>
#include <boost/align/aligned_allocator.hpp>

#include <boost/simd/function/inc.hpp>
#include <boost/simd/function/log2.hpp>
#include <boost/simd/function/exp2.hpp>


namespace nupack {

/******************************************************************************************/

struct Zero {
    template <class T>
    friend T operator+(T const &t, Zero) {return t;}

    template <class T>
    friend T operator+(Zero, T const &t) {return t;}

    friend Zero operator+(Zero, Zero) {return {};}

    friend False operator<(Zero, Zero) {return {};}

    template <class T>
    friend Zero operator*(T const &t, Zero) {return {};}

    template <class T>
    friend Zero operator*(Zero, T const &t) {return {};}

    friend Zero operator*(Zero, Zero) {return {};}

    template <class T>
    friend T operator-(T const &t, Zero) {return t;}

    template <class T>
    friend T operator-(Zero, T const &t) {return -t;}

    friend Zero operator-(Zero, Zero) {return {};}

    Zero operator-() const {return {};}

    constexpr operator std::int32_t() const {return 0;}
    constexpr operator std::int64_t() const {return 0;}
};

NUPACK_UNARY_FUNCTOR(always_zero_exp, Zero());

/******************************************************************************************/

namespace simd {

using namespace boost::simd;
namespace bs = boost::simd;

static auto ifrexp = bs::pedantic_(bs::ifrexp);

template <class T> static constexpr bool can_simd = (sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8);

/// pack_size is the default number of elements in a SIMD pack suggested by Boost.SIMD
template <class T>
static constexpr auto pack_size = 2 * bs::pack<decay<T>>::static_size;

/******************************************************************************************/

/* Chunk is used as an element accessor, i.e. array[Chunk()]. It tells the container
   to return N elements starting at a given position  */
template <int N>
struct Chunk {
    static constexpr auto length = N;
    int value;
    explicit constexpr Chunk(int v) : value(v) {}

    friend auto operator-(Chunk c, int i) noexcept {return Chunk(c.value-i);}
    friend auto operator-(int i, Chunk c) noexcept {return Chunk<-N>(i-c.value);}
    friend auto operator+(int i, Chunk c) noexcept {return Chunk(c.value+i);}
    friend auto operator+(Chunk c, int i) noexcept {return Chunk(c.value+i);}
};

NUPACK_DEFINE_TEMPLATE(is_chunk, Chunk, int);

template <class T>
struct pack_element_type {using type = T;};

template <class T, std::size_t N>
struct pack_element_type<bs::pack<T, N>> {using type = T;};

template <class T> using pack_element_t = typename pack_element_type<decay<T>>::type;

/******************************************************************************************/

template <int I, class T, NUPACK_IF(I >= +1)>
auto load_pack(T const *t) noexcept {return bs::load<bs::pack<T, I>>(t);}

template <int I, class T, NUPACK_IF(I <= -1)>
auto load_pack(T const *t) noexcept {return bs::reverse(load_pack<-I>(t + (1 + I)));}

/******************************************************************************************/

/// Perform a map-reduce operation with SIMD, where the output is modified in place
template <class R, class D, class M, class Op>
auto map_reduce(R reduce, D domain, M map, Op op) {
    auto it = begin_of(domain);
    auto ret = map(*it++);
    constexpr auto Z = pack_size<decltype(ret)>;
    if (it + Z <= end_of(domain)) {
        auto sum = map(Chunk<Z>(*it));
        for (it += Z; it + Z <= end_of(domain); it += Z) reduce(sum, map(Chunk<Z>(*it)));
        reduce(ret, op(sum));
    }
    for (; it != end_of(domain); ++it) reduce(ret, map(*it));
    return ret;
}

/// Shorthand for the default allocator that should be used
template <class T>
using allocator = boost::alignment::aligned_allocator<T, bs::pack<T>::alignment>;

/******************************************************************************************/

struct times_t {
    template <class T, class U, NUPACK_IF(!is_pair<T> && !is_pair<U>)>
    auto operator()(T const &t, U const &u) const noexcept {return bs::multiplies(t, u);}

    template <class T, class U, NUPACK_IF(is_pair<T> && !is_pair<U>)>
    constexpr auto operator()(T const &t, U const &u) const noexcept {return std::make_pair(t.first * u, t.second);}

    template <class T, class U, NUPACK_IF(!is_pair<T> && is_pair<U>)>
    constexpr auto operator()(T const &t, U const &u) const noexcept {return std::make_pair(t * u.first, u.second);}

    template <class T, class U, NUPACK_IF(is_pair<T> && is_pair<U>)>
    constexpr auto operator()(T const &t, U const &u) const noexcept {return std::make_pair(t.first * u.first, t.second + u.second);}
};

static constexpr auto times = times_t();

/******************************************************************************************/

/// reciprocal functor: f(x) = 1/x
struct invert_t {
    template <class T, NUPACK_IF(is_pair<T>)>
    constexpr auto operator()(T const &t) const noexcept {return std::make_pair(bs::rec(t.first), bs::unary_minus(t.second));}

    template <class T, NUPACK_IF(!is_pair<T>)>
    constexpr auto operator()(T const &t) const noexcept {return bs::rec(t);}
};

static constexpr auto invert = invert_t();

/******************************************************************************************/

struct max_t {
    constexpr auto operator()(Zero, Zero) const noexcept {return Zero();}
    constexpr auto operator()(Zero, Zero, Zero) const noexcept {return Zero();}

    template <class T, class U>
    constexpr decltype(auto) operator()(T const &t, U const &u) const noexcept {return bs::max(t, u);}

    template <class T, class U, class V>
    constexpr decltype(auto) operator()(T const &t, U const &u, V const &v) const noexcept {return bs::max(bs::max(t, u), v);}
};

static constexpr auto max = max_t();

/******************************************************************************************/

/// plus equals functor: f(x) = x + y
struct plus_t {
    template <class T, class U, NUPACK_IF(!is_same<Zero, T, U> && !is_pair<T>)>
    decltype(auto) operator()(T &&t, U &&u) const noexcept {
        static_assert(is_floating_point<pack_element_t<T>> || is_signed<pack_element_t<T>>, "first");
        static_assert(is_floating_point<pack_element_t<U>> || is_signed<pack_element_t<U>>, "second");
        return bs::plus(t, u);
        // return bs::plus(fw<T>(t), fw<U>(u));
    }

    template <class T, NUPACK_IF((!is_same<T, Zero>))>
    constexpr T operator()(T &&t, Zero) const noexcept {return t;}

    template <class T, NUPACK_IF((!is_same<T, Zero>))>
    constexpr T operator()(Zero, T &&t) const noexcept {return t;}

    constexpr Zero operator()(Zero, Zero) const noexcept {return {};}

    template <class T, class U, NUPACK_IF(is_pair<T>)>
    constexpr auto operator()(T &&t, U &&u) const noexcept {
        return (t.second < u.second) ? std::make_pair(u.first + simd::ldexp(t.first, t.second-u.second), u.second)
                                     : std::make_pair(t.first + simd::ldexp(u.first, u.second-t.second), t.second);
    }
};
static constexpr auto plus = plus_t();

/******************************************************************************************/

struct lse2_t {
    template <class T, class U>
    decltype(auto) operator()(T const &t, U const &u) const noexcept {
        return bs::if_else(t < u, u + bs::log2(bs::inc(bs::exp2(t - u))),
                                  t + bs::log2(bs::inc(bs::exp2(u - t))));
    }
    template <class T>
    decltype(auto) operator()(T const &t) const noexcept {
        auto x = bs::maximum(t);
        return x + bs::log2(bs::sum(bs::exp2(t - x)));
    }
};

static constexpr auto lse2 = lse2_t();

/******************************************************************************************/

struct min_t {
    template <class T, class U, NUPACK_IF(!is_same<Zero, decay<T>, decay<U>>)>
    decltype(auto) operator()(T &&t, U &&u) const noexcept {return bs::min(fw<T>(t), fw<U>(u));}

    Zero operator()(Zero, Zero) const noexcept {return Zero();}
};
static constexpr auto min = min_t();

/******************************************************************************************/

/// ldexp: x = min(x, y)
struct min_eq_t {
    template <class T, class U>
    decltype(auto) operator()(T &t, U &&u) const noexcept {return t = min(t, fw<U>(u));}
};
static constexpr auto min_eq = min_eq_t();

/******************************************************************************************/

/// ldexp: return f(x, y) = x * 2^y
struct ldexp_t {
    template <class T> constexpr T operator()(T &&t, Zero) const noexcept {return static_cast<T &&>(t);}

    template <class T, class U, NUPACK_IF(!is_same<decay<U>, Zero>)>
    auto operator()(T &&t, U &&u) const noexcept {
        constexpr pack_element_t<U> u0 = std::numeric_limits<pack_element_t<T>>::min_exponent * 3 / 4;
        // underflow breaks this so floor it; simd::ldexp gives non compliant
        // answer if x = 0 and 2^y -> inf (should be 0, but returns nan)
        return bs::if_zero_else(bs::is_eqz(t), bs::ldexp(fw<T>(t), bs::max(fw<U>(u), u0)));
    }
};

static constexpr auto ldexp = ldexp_t();

/******************************************************************************************/

}}
