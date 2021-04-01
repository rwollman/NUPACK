/**
 * @brief C++ functor signature deduction
 * @file Signature.h
 */

#pragma once
#include "Type.h"

namespace rebind {

/******************************************************************************************/

struct NotFound {};

NotFound operator|(NotFound, NotFound);

template <class T>
T operator|(NotFound, T);

template <class T>
T operator|(T, NotFound);

/// A lightweight ordered container of types
template <class ...Ts> struct Pack {
    using pack_type = Pack;

    static constexpr auto size = sizeof...(Ts);
    using indices = decltype(std::make_index_sequence<sizeof...(Ts)>());

    template <class T>
    static constexpr bool contains = (std::is_same_v<T, Ts> || ...);

    template <std::size_t N, std::size_t ...Is>
    static constexpr auto at_sequence(std::index_sequence<Is...>) {
        return decltype((std::conditional_t<N == Is, Type<Ts>, NotFound>() | ...))();
    }

    template <std::size_t N>
    static constexpr auto at() {return at_sequence<N>(indices());}

    template <class N>
    static constexpr auto at(N) {return at_sequence<N::value>(indices());}

    template <class T, std::size_t ...Is>
    static constexpr auto position_impl(std::index_sequence<Is...>) {
        return decltype((std::conditional_t<std::is_same_v<T, Ts>,
            std::integral_constant<std::size_t, Is>, NotFound>() | ...))();
    }

    template <class T>
    static constexpr auto position(Type<T> t={}) {return position_impl<T>(indices());}

    template <std::size_t B, std::size_t ...Is>
    static constexpr auto slice_sequence(std::index_sequence<Is...>) {return Pack<typename decltype(at<B + Is>())::type...>();}

    template <std::size_t B, std::size_t E>
    static constexpr auto slice() {return slice_sequence<B>(std::make_index_sequence<E - B>());}

    template <class F, std::size_t ...Is>
    static constexpr auto indexed(F &&f, std::index_sequence<Is...>) {
        return static_cast<F &&>(f)(IndexedType<Ts>{Is}...);
    }

    template <class F>
    static constexpr auto indexed(F &&f) {return indexed(static_cast<F &&>(f), indices());}

    template <class F>
    static constexpr decltype(auto) apply(F &&f) {return f(Type<Ts>()...);}

    template <class F>
    static void for_each(F &&f) {(f(Type<Ts>()), ...);}

    using unqualified = Pack<std::remove_cv_t<std::remove_reference_t<Ts>>...>;
};

template <class ...Ts, class ...Us>
Pack<Ts..., Us...> concat(Pack<Ts...>, Pack<Us...>) {return {};}

template <template <class ...> class T, class ...Ts>
Pack<Ts...> template_as_pack_impl(T<Ts...> const &);

template <class T>
using template_as_pack = decltype(template_as_pack_impl(std::declval<T>()));

/******************************************************************************************/

template <class F, class ...Ts>
void scan_packs(Pack<Ts...>, F const &f) {f(Type<Ts>()...);}

template <class F, class P, class ...Ts, class ...Us, class ...Fs>
void scan_packs(Pack<Ts...>, Pack<Us...>, P const &p, Fs const &...fs) {
    (scan_packs(Pack<Ts..., Us>(), p, fs...), ...);
}

template <class ...Ts, class ...Fs>
void scan(Pack<Ts...> p, Fs const &...fs) {return scan_packs(Pack<>(), p, fs...);}

/******************************************************************************************/

template <class F, class=void>
struct Signature;

template <class R, class ...Ts>
struct Signature<R(Ts...)> : Pack<R, Ts...> {using return_type = R;};
template <class R, class ...Ts>
struct Signature<R(Ts...) noexcept> : Pack<R, Ts...> {using return_type = R;};

template <class R, class ...Ts>
struct Signature<R(*)(Ts...)> : Signature<R(Ts...)> {using return_type = R;};
template <class R, class ...Ts>
struct Signature<R(*)(Ts...) noexcept> : Signature<R(Ts...)> {using return_type = R;};

#define REBIND_TMP(C, Q, C2) \
    template <class R, class C, class ...Ts> \
    struct Signature<R (C::* )(Ts...) Q> : Pack<R, C2, Ts...> {using return_type = R;};

    REBIND_TMP(C, , C &);
    REBIND_TMP(C, const, C const &);
    REBIND_TMP(C, &, C &);
    REBIND_TMP(C, const &, C const &);
    REBIND_TMP(C, &&, C &&);
    REBIND_TMP(C, const &&, C const &&);
#undef REBIND_TMP

/// this is tricky...
template <class R, class C>
struct Signature<R C::*> : Pack<R &, C &> {using return_type = R &;};

/******************************************************************************************/

template <class T>
struct FunctorSignature;

#define REBIND_TMP(Q) template <class R, class C, class ...Ts> \
    struct FunctorSignature<R (C::* )(Ts...) Q> : Signature<R(Ts...)> {using return_type = R;};
    REBIND_TMP( );
    REBIND_TMP(&);
    REBIND_TMP(&&);
    REBIND_TMP(const);
    REBIND_TMP(const &);
    REBIND_TMP(const &&);
#undef REBIND_TMP

/******************************************************************************************/

template <class F>
struct Signature<F, std::void_t<decltype(&F::operator())>> : FunctorSignature<decltype(&F::operator())> {};

template <class F>
constexpr typename Signature<F>::pack_type signature(F const &) {return {};}

/******************************************************************************************/

template <class F, class=void>
struct SimplifyFunction {
    constexpr std::decay_t<F> operator()(F f) const {return f;}
};

template <class F>
struct SimplifyFunction<F, std::void_t<decltype(false ? nullptr : std::declval<F>())>> {
    constexpr auto operator()(F f) const {return false ? nullptr : f;}
};

}
