/**
 * @brief Common functors, operators, CRTP-like base classes
 *
 * @file Operators.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "TypeSupport.h"
#include "Extents.h"
#include "../Forward.h"

#include <functional>
#include <cmath>

namespace nupack {

/******************************************************************************************/

/// Default unary operation which returns the input itself
struct Identity {template <class T> constexpr remove_rref<T &&> operator() (T &&t) const {return fw<T>(t);};};
struct NoOp {template <class ...Ts> void operator() (Ts const &...) const {};};
struct AlwaysTrue {template <class ...Ts> constexpr bool operator() (Ts const &...) const {return true;};};
struct AlwaysFalse {template <class ...Ts> constexpr bool operator() (Ts const &...) const {return false;};};


/******************************************************************************************/

/// Struct allowing either unary or binary operator to be used for reducing operation
template <class Op, class F, int Mode> struct Reduce_Operator {
    F f;
    // Overload allowing f to be called with multiple arguments
    template <class T, class U, NUPACK_IF(can_call<F const &, T &&, U &&>)>
    constexpr decltype(auto) operator() (T &&t, U &&u) const {return f(fw<T>(t), fw<U>(u));}
    // Overload allowing f to be called with single arguments
    template <class T, class U, NUPACK_IF(!can_call<F const &, T &&, U &&> && Mode == 1)>
    constexpr decltype(auto) operator() (T &&t, U &&u)  const {return Op()(f(fw<T>(t)), f(fw<U>(u)));}
    /// Overload for f called only on the second argument
    template <class T, class U, NUPACK_IF(!can_call<F const &, T &&, U &&> && Mode == 0)>
    constexpr decltype(auto) operator() (T &&t, U &&u)  const {return Op()(fw<T>(t), f(fw<U>(u)));}
};

/// Make a binary operator from a unary or binary one
template <class Op, class F> constexpr auto reduce_op(F &&f) {
    return Reduce_Operator<Op, decay<decltype(f)>, 1>{fw<F>(f)};
}

/// Make a binary operator from a unary or binary one
template <class Op, class F> constexpr auto update_op(F &&f) {
    return Reduce_Operator<Op, decay<decltype(f)>, 0>{fw<F>(f)};
}

/******************************************************************************************/

/// Struct allowing either an object or function to be used for unary operation
template <class Op, class F> struct PredicateOperator {
    F f;
    // Overload allowing f to be called with multiple arguments
    template <class T, NUPACK_IF(can_call<F const &, T &&>)>
    constexpr decltype(auto) operator() (T &&t) const {return f(fw<T>(t));}
    // Overload allowing f to be called with single arguments
    template <class T, NUPACK_IF(!can_call<F const &, T &&>)>
    constexpr decltype(auto) operator() (T &&t) const {return Op()(f, fw<T>(t));}
};

/// Make a binary operator from a unary or binary one
template <class Op, class F> constexpr auto predicate_op(F &&f) {return PredicateOperator<Op, F>{fw<F>(f)};}

/******************************************************************************************/

/// Given a class with ==, defines !=
struct EqualComparable {};
template <class T, NUPACK_IF(derives_from<T, EqualComparable>)>
constexpr bool operator!=(T const &i, T const &j)  {return !(i == j);}

/// Given a class with <, define >
struct WeaklyOrdered {};
template <class T, NUPACK_IF(derives_from<T, WeaklyOrdered>)>
constexpr bool operator>(T const &i, T const &j)  {return j < i;}

/// Given a class with < and ==, defines all other operators
struct TotallyOrdered : WeaklyOrdered, EqualComparable {};
template <class T, NUPACK_IF(derives_from<T, TotallyOrdered>)>
constexpr bool operator<=(T const &i, T const &j) {return !(j < i);}
template <class T, NUPACK_IF(derives_from<T, TotallyOrdered>)>
constexpr bool operator>=(T const &i, T const &j) {return !(i < j);}

/******************************************************************************************/

/// If a class supplies an iter() const function returning an Iterable
/// supply mutating begin and end functions for it
template <class T> class ConstIterable {
    constexpr decltype(auto) get_iter() const {return static_cast<T const *>(this)->iter();}

public:
    constexpr auto begin() const {return begin_of(get_iter());}
    constexpr auto end() const {return end_of(get_iter());}

    constexpr auto size() const {return len(get_iter());}
};

namespace detail {
    template <class U> U const & with_const(U &);
    template <class U> U const & with_const(U const &);
    template <class U> U const with_const(U &&);
}

/// If a class supplies an iter() const function returning an Iterable
/// supply mutating begin and end functions for it

template <class T> class Iterable {
    constexpr decltype(auto) get_iter() {return static_cast<T *>(this)->iter();}
    /// accesses this constly, then removes the const, then gets iter(), then puts const on that thing
    constexpr decltype(auto) get_iter() const {return static_cast<decltype(detail::with_const(declref<T>().iter()))>(const_cast<T *>(static_cast<T const *>(this))->iter());}

public:
    constexpr auto begin() const {return begin_of(get_iter());}
    constexpr auto end() const {return end_of(get_iter());}
    constexpr auto begin() {return begin_of(get_iter());}
    constexpr auto end() {return end_of(get_iter());}
    constexpr auto size() const {return len(get_iter());}
};

/******************************************************************************************/

/// If a class supplies an iter() const function returning an Iterable
/// supply const begin and end functions for it
template <class T> struct ConstIndexable : ConstIterable<T> {
    template <class I> auto operator[] (I &&i) const -> decltype(declval<sink_sfinae<T, I>>().iter()[fw<I>(i)]) {
        return static_cast<T const *>(this)->iter()[fw<I>(i)];
    };
};

/******************************************************************************************/

/// If a class supplies an iter() const function returning an Iterable
/// supply mutating begin and end functions for it
template <class T> struct Indexable : Iterable<T> {
    template <class I> auto operator[] (I &&i) const -> decltype(detail::with_const(declref<sink_sfinae<T, I>>().iter()[fw<I>(i)])) {
        return add_const(const_cast<T *>(static_cast<T const *>(this))->iter()[fw<I>(i)]);
    };

    template <class I> auto operator[] (I &&i) -> decltype(declref<sink_sfinae<T, I>>().iter()[fw<I>(i)]) {
        return static_cast<T *>(this)->iter()[fw<I>(i)];
    };
};

/******************************************************************************************/

/// Helper for approximate comparison
template <class T> struct About {
    T value, epsilon, scale;

    About(T t, T e=std::sqrt(std::numeric_limits<T>::epsilon()) * 1000, T s=1) : value(t), epsilon(e), scale(s) {}

    friend bool operator==(About const &a, T const &t) {
        if (!std::isfinite(t) || !std::isfinite(a.value)) return t == a.value;
        return std::abs(t - a.value) < a.epsilon * (a.scale + std::max(std::abs(t), std::abs(a.value)));
    }

    friend bool operator==(T const &t, About const &a) {return a == t;}
    friend bool operator!=(T const &t, About const &a) {return !(a == t);}
    friend bool operator!=(About const &a, T const &t) {return !(a == t);}

    friend std::ostream & operator<<(std::ostream &os, About const &t) {return os << t.value;}
};

template <class T, class ...Us>
About<if_t<is_floating_point<T>, T, double>> about(T const &t, Us &&...us) {
    using R = if_t<is_floating_point<T>, T, double>;
    return {static_cast<R>(t), static_cast<R>(us)...};
}

/******************************************************************************************/

}
