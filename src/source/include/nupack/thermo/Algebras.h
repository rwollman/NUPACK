/**
 * @brief Dynamic program algebras for forward and backward recursions
 *
 * @file Algebras.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Rigs.h"
#include "../standard/Optional.h"

namespace nupack::thermo {

/******************************************************************************************/

template <class M, class X>
struct Expression {
    M mantissa;
    X exponent;

    template <class T=Zero>
    friend auto mantissa(Expression const & e, T hint={}) {return e.mantissa(hint);}

    template <class T=Zero>
    friend auto exponent(Expression const & e, T hint={}) {return e.exponent(hint);}
};

template <class M, class X>
auto expression(M &&mantissa, X &&exponent) {return Expression<decay<M>, decay<X>>{fw<M>(mantissa), fw<X>(exponent)};}

template <class M>
auto expression(M &&mantissa) {return expression(fw<M>(mantissa), always_zero_exp);}

NUPACK_DEFINE_TEMPLATE(is_expression, Expression, class, class);

/******************************************************************************************/

template <class Rig>
struct ForwardAlgebra {
    using rig_type = Rig;
    using is_forward = True;

    template <class T, class U>
    constexpr auto ldexp(T &&t, U &&u) const {return Rig::ldexp()(fw<T>(t), fw<U>(u));}

    template <class ...Ts>
    constexpr auto times(Ts &&...ts) const {return Rig::times()(value_of(fw<Ts>(ts))...);}

    template <class ...Ts>
    constexpr auto plus(Ts &&...ts) const {return fold(Rig::plus(), value_of(fw<Ts>(ts))...);}

    /// (Ts +...) where all the Ts are scalars
    template <class ...Ts> auto sum(Ts const &...ts) const {
        return expression([=] (auto hint) {return fold(Rig::plus(), Rig::ldexp()(mantissa(ts, hint), exponent(ts, hint))...);});
    }

    /**************************************************************************************/

    /// (Ts *...) where all the Ts are scalars
    template <class T, class ...Ts> auto product(T const &t, Ts const &...ts) const {
        return expression([=] (auto hint) {
            static_assert(is_same<decltype(hint), Zero> || is_signed<decltype(hint)>, "");
            auto new_hint = fold(Rig::plus(), hint, exponent(ts)...);
            return Rig::ldexp()(fold(Rig::times(), mantissa(t, new_hint), mantissa(ts)...), exponent(t, new_hint));
        });
    }

    /**************************************************************************************/

    // sum(range, f) sums the scalar outputs of f on the elements of the range
    template <class T, class F> auto total(T const &t, F &&f) const {
        return expression([=] (auto hint) {
            decay<decltype(mantissa(f(*begin_of(t)), hint))> ret = Rig::zero();
            for (auto i : t) {
                auto elem = f(i);
                Rig::plus_eq()(ret, Rig::ldexp()(mantissa(elem, hint), exponent(elem, hint)));
            }
            return ret;
        });
    }

    // sum(Ts[:] *...)
    template <class T, class ...Ts> auto dot(T const &t, Ts const &...ts) const {
        return expression([=] (auto hint) {
            auto map = [&](auto i) {
                return Rig::ldexp()(fold(Rig::times(), mantissa_at(t, i), mantissa_at(ts, i)...),
                                    fold(Rig::plus(), hint, exponent_at(t, i), exponent_at(ts, i)...));
            };
            return simd::map_reduce(Rig::plus_eq(), indices(t), std::move(map), Rig::sum());
        });
    }

    /**************************************************************************************/

    struct zero_t : decltype(Rig::zero()) {};
    constexpr auto zero() const {return zero_t();}

    template <class T>
    struct Conditional : Optional<T> {
        Conditional(T const &t) : Optional<T>(t) {}
        Conditional(T &&t) : Optional<T>(std::move(t)) {}
        Conditional(zero_t) {}

        template <class E=Zero>
        friend auto mantissa(Conditional const &c, E e={}) {return c ? mantissa(*c, e) : Rig::zero();}
        template <class E=Zero>
        friend auto exponent(Conditional const &c, E e={}) {return c ? exponent(*c, e) : *::nupack::zero;}
    };

    struct Maybe {
        template <class T>
        friend constexpr Conditional<decay<T>> operator&(Maybe, T &&t) {return fw<T>(t);}
    };

    constexpr auto maybe() const {return Maybe();}
};

/**************************************************************************************/

struct Recursible_Base {};
template <class F> struct Recursible : Recursible_Base {
    F expression;
    Recursible(F f) : expression(std::move(f)) {}
};
template <class F> auto recursible(F &&f) {
    static_assert(is_class<decay<F>>, "should pass in functor");
    return Recursible<decay<F>>{fw<F>(f)};
}
NUPACK_DETECT(is_recursible, void_if<(derives_from<T, Recursible_Base>)>);

/**************************************************************************************/

template <class Rig, bool Short_Circuit>
struct BackwardAlgebra {
    using is_forward = False;
    using rig_type = Rig;
    // In the following, fun is assumed to be (result, expressions...) -> bool
    // "recurse" is called with fun and a variadic number of expressions
    // The first expression will not be directly evaluated if it is recursible
    // All other expressions will be directly evaluated

    /// This returns if we should stop or not. If not short circuiting, return compile time False
    static constexpr auto short_circuit(bool b) {return if_c<Short_Circuit>(b, False());}

    // If t is a compound object, recurse into it, passing along the evaluated multiplying factors
    template <class F, class T, class ...Ts, NUPACK_IF(traits::is_recursible<T>)>
    static auto recurse(F &&fun, T const &t, Ts const &...ts) {
        return short_circuit(t.expression(fun, ts...));
    }

    // If t is not a compound object, just multiply and pass arguments to fun
    template <class F, class T, class ...Ts>
    static auto simple_recurse(F &&fun, T const &t, Ts const &...ts) {
        return short_circuit(fun([&] (auto hint) {
            auto new_hint = fold(Rig::plus(), hint, exponent(ts)...);
            return Rig::ldexp()(fold(Rig::times(), mantissa(t, new_hint), mantissa(ts)...), exponent(t, new_hint));
        }, t, ts...));
    }

    template <class F, class T, class ...Ts, NUPACK_IF(!traits::is_recursible<T>)>
    static auto recurse(F &&fun, T const &t, Ts const &...ts) {return simple_recurse(fun, t, ts...);}

    struct zero_t : decltype(Rig::zero()) {};
    constexpr auto zero() const {return zero_t();}

    template <class T>
    struct Conditional : Optional<T>, Recursible_Base {
        Conditional(T &&t) : Optional<T>(std::move(t)) {}
        Conditional(zero_t) {}

        template <class F, class ...Ts>
        bool expression(F &&f, Ts &&...ts) const {return *this && recurse(fw<F>(f), **this, fw<Ts>(ts)...);}
    };

    struct Maybe {
        template <class T>
        friend constexpr Conditional<decay<T>> operator&(Maybe, T &&t) {return fw<T>(t);}
    };

    constexpr auto maybe() const {return Maybe();}

    template <class ...Ts> auto product(Ts const &...ts) const {
        return recursible([=](auto &&fun, auto const &...factors) {return recurse(fun, ts..., factors...);});
    }

    // From variadic possibilities, recurse at most a single contribution
    template <class ...Ts> auto sum(Ts const &...ts) const {
        return recursible([=](auto &&fun, auto const &...factors) {
            bool b = false;
            NUPACK_UNPACK(b = b || recurse(fun, ts, factors...));
            return b;
        });
    }

    // From an Iterable recurse at most a single contribution
    template <class T, class F> auto total(T const &t, F &&map) const {
        return recursible([=](auto &&fun, auto const &...factors) {
            for (auto i : t) if (recurse(fun, map(i), factors...)) return true;
            return false;
        });
    }

    // From a range of products recurse a single contribution
    template <class T, class ...Ts>
    auto dot(T const &t, Ts const &...ts) const {
        return recursible([=](auto &&fun, auto const &...factors) {
            for (auto i : indices(t))
                if (recurse(fun, next(t, i), next(ts, i)..., factors...)) return true;
            return false;
        });
    }
};

template <class Rig> using SuboptAlgebra = BackwardAlgebra<Rig, false>;
template <class Rig> using SampleAlgebra = BackwardAlgebra<Rig, true>;

/******************************************************************************************/

}
