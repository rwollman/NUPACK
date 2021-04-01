/**
 * @brief Compile-time constants; pretty rudimentary and should be swapped out (maybe for boost hana)
 *
 * @file Constants.h
 * @author Mark Fornace
 * @date 2018-05-18
 */
#pragma once
#include "Traits.h"
#include <limits>

namespace nupack {

/******************************************************************************************/

namespace detail {
    template <class T, class N, class=void> struct DefaultConstantConvert;
}

template <class T, class N, class=void> struct ConvertConstant : std::false_type {};

template <class T, class N> using ConstantConverter = if_t<derives_from<ConvertConstant<T, N>, std::false_type>,
    detail::DefaultConstantConvert<T, N>, ConvertConstant<T, N>>;

/// compile-time constant that can be explicitly converted to numeric types, generally
template <class N>
struct Constant : N {
    template <class T, NUPACK_IF(can_convert<decltype(ConstantConverter<T, N>()()), T>)>
    explicit constexpr operator T() const {return static_cast<T>(ConstantConverter<T, N>()());}

    /// provide another type that is implicitly convertible in the same way
    struct implicit {
        template <class T, NUPACK_IF(can_convert<decltype(ConstantConverter<T, N>()()), T>)>
        constexpr operator T() const {return static_cast<T>(ConstantConverter<T, N>()());}
        constexpr operator Constant() const {return Constant();}

//#define NUPACK_TEMP(op) template <class T> \
//friend constexpr auto operator op(T const t, implicit i) -> decltype(t op static_cast<T>(i)) {return t op static_cast<T>(i);} \
//template <class T> \
//friend constexpr auto operator op(implicit i, T const t) -> decltype(static_cast<T>(i) op t) {return static_cast<T>(i) op t;}
//    NUPACK_TEMP(+); NUPACK_TEMP(-); NUPACK_TEMP(*); NUPACK_TEMP(/); NUPACK_TEMP(|);
//    NUPACK_TEMP(<<); NUPACK_TEMP(>>); NUPACK_TEMP(<); NUPACK_TEMP(>);
//    NUPACK_TEMP(<=); NUPACK_TEMP(>=); NUPACK_TEMP(==); NUPACK_TEMP(!=);
//#undef NUPACK_TEMP
    };

    constexpr auto operator*() const {return implicit();}
};

/******************************************************************************************/

template <int N> using size_constant = std::integral_constant<std::size_t, N>;
template <int N> using int_constant = std::integral_constant<int, N>;

struct True : Constant<std::true_type> {};
static constexpr auto const true_c = True{};
struct False : Constant<std::false_type> {};
static constexpr auto const false_c = False{};
template <bool B> using bool_t  = std::conditional_t<B, True, False>;

struct zero_t : Constant<size_constant<0>> {
    constexpr zero_t operator-() const {return {};}
};
static constexpr auto const zero = zero_t{};

struct AlwaysZero {template <class ...Ts> constexpr auto operator()(Ts const &...) const {return *zero;}};

struct one_t : Constant<size_constant<1>> {};
static constexpr auto const one = one_t{};
struct two_t : Constant<size_constant<2>> {};
static constexpr auto const two = two_t{};
struct three_t : Constant<size_constant<3>> {};
static constexpr auto const three = three_t{};

struct inf_tag {static constexpr auto value = std::numeric_limits<std::ptrdiff_t>::max();};
using inf_t = Constant<inf_tag>;
static constexpr auto const inf = inf_t{};

struct minf_tag {static constexpr auto value = std::numeric_limits<std::ptrdiff_t>::lowest();};
using minf_t = Constant<minf_tag>;
static constexpr auto const minf = minf_t{};

//#define NUPACK_COP(op) template <class M, class N> \
//constexpr auto operator op(Constant<M> c, Constant<N> d) {return bool_t<decltype(c)::value op decltype(d)::value>{};}
//NUPACK_COP(<); NUPACK_COP(>); NUPACK_COP(<=); NUPACK_COP(>=); NUPACK_COP(==); NUPACK_COP(!=);
//#undef NUPACK_COP

namespace detail {

/******************************************************************************************/

template <class T> struct DefaultConstantConvert<T, std::false_type, void_if<can_convert<bool, T>>> {
    constexpr T operator()() const {return static_cast<T>(false);}
};

template <class T> struct DefaultConstantConvert<T, std::true_type, void_if<can_convert<bool, T>>> {
    constexpr T operator()() const {return static_cast<T>(true);}
};

template <class T> struct DefaultConstantConvert<T, size_constant<0>, void_if<is_arithmetic<T> && can_convert<std::size_t, T>>> {
    constexpr T operator()() const {return static_cast<T>(0);}
};

template <class T> struct DefaultConstantConvert<T, size_constant<0>, void_if<!(is_arithmetic<T> && can_convert<std::size_t, T>) && can_construct<T>>> {
    constexpr T operator()() const {return T{};}
};

template <std::size_t I, class T> struct DefaultConstantConvert<T, size_constant<I>, void_if<I && is_arithmetic<T>>> {
    constexpr T operator()() const {return static_cast<T>(I);}
};

template <class T> struct DefaultConstantConvert<T, inf_tag, void_if<is_floating_point<T>>> {
    constexpr T operator()() const {return std::numeric_limits<T>::infinity();}
};

template <class T> struct DefaultConstantConvert<T, inf_tag, void_if<is_integral<T>>> {
    constexpr T operator()() const {return std::numeric_limits<T>::max();}
};

template <class T> struct DefaultConstantConvert<T, minf_tag, void_if<is_floating_point<T>>> {
    constexpr T operator()() const {return -std::numeric_limits<T>::infinity();}
};

template <class T> struct DefaultConstantConvert<T, minf_tag, void_if<is_integral<T>>> {
    constexpr T operator()() const {return std::numeric_limits<T>::lowest();}
};

/******************************************************************************************/

}

};
