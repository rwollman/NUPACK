/**
 * @brief Definition of Pack<Ts...>, a container of types
 *
 * @file Pack.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Constants.h"

namespace nupack {

/// A lightweight ordered container of types
template <class ...Ts> struct pack {
    /// Element access by index
    template <std::size_t I>
    using type = type_at<I, Ts...>;
    /// First index of type
    template <class U>
    using index = size_constant<index_of_type<U, Ts...>>;
    /// Number of elements
    using size = size_constant<sizeof...(Ts)>;
    /// Test if it contains a type
    template <class U> using has = bool_t<is_same<U, Ts...>>;
    /// Test if the type is in the sequence but is not the first one
    template <class U> using not_first = bool_t<has<U>::value && !std::is_same<U, type<0>>::value>;
    /// Call a functor with a type_t of each type
    template <class F>
    static void for_each(F &&f) {NUPACK_UNPACK(f(type_t<Ts>()));}
    /// Call a functor with a type_t of each type while it returns true
    template <class F>
    static bool while_each(F &&f) {NUPACK_WHILE(f(type_t<Ts>())); return nupack_while_ok;}

    template <class F>
    static constexpr auto apply(F &&f) {return fw<F>(f)(type_t<Ts>()...);}

    static constexpr auto indices() {return std::make_index_sequence<sizeof...(Ts)>();}

    template <class I>
    static constexpr type_t<type_at<I::value, Ts...>> at(I) {return {};}

    template <class I>
    constexpr type_t<type_at<I::value, Ts...>> operator[](I) const {return {};}
};

template <> struct pack<> {
    using size = size_constant<0>;
    template <class U> using has = False;
    template <class U> using not_first = False;
    template <class F> static void for_each(F &&) {}
    template <class F> static constexpr bool while_each(F &&) {return true;}
    template <class F> static constexpr auto apply(F &&f) {return fw<F>(f)();}
    static constexpr auto indices() {return std::make_index_sequence<0>();}
};

template <class ...Ts> using pack_lref = pack<Ts &...>;
template <class ...Ts> using pack_rref = pack<Ts &&...>;
template <class ...Ts> using pack_cref = pack<Ts const &...>;
template <class ...Ts> using pack_const = pack<Ts const...>;

/******************************************************************************************/

template <class ...Ts, class ...Us>
constexpr pack<Ts..., Us...> operator+(pack<Ts...>, pack<Us...>) {return {};}

NUPACK_DEFINE_VARIADIC(is_pack, pack, class);

/******************************************************************************************/

namespace detail {
    template <class T> struct as_pack_t {using type = pack<T>;};
    template <class ...Ts> struct as_pack_t<std::tuple<Ts...>> {using type = pack<Ts...>;};
    template <class ...Ts> struct as_pack_t<pack<Ts...>> {using type = pack<Ts...>;};
}

template <class T> using as_pack = typename detail::as_pack_t<decay<T>>::type;

/******************************************************************************************/

struct not_found {};

template <class F, class I>
not_found find_c(pack<>, F const &, I) {return {};}

template <class T, class ...Ts, class F, class I=size_constant<0>, NUPACK_IF(decltype(declref<F const>()(type_t<T>()))::value)>
static I find_c(pack<T, Ts...>, F const &, I={}) {return {};}

template <class T, class ...Ts, class F, class I=size_constant<0>, NUPACK_IF(!decltype(declref<F const>()(type_t<T>()))::value)>
static auto find_c(pack<T, Ts...>, F const &f, I={}) {return decltype(find_c(pack<Ts...>(), f, size_constant<I::value+1>()))();}

/******************************************************************************************/

template <class F, class=void>
struct signature_t;

NUPACK_DETECT(has_signature, typename signature_t<T>::return_type);

template <class R, class ...Ts>
struct signature_t<R(Ts...)> : pack<R, Ts...> {
    using size = size_constant<sizeof...(Ts) + 1>;
    using return_type = R;
    using pack_type = pack<R, Ts...>;

    template <class F> static void for_each_arg(F &&f) {NUPACK_UNPACK(f(type_t<Ts>()));}
    template <class F> static void while_each_arg(F &&f) {NUPACK_WHILE(f(type_t<Ts>()));}
};

template <class R, class ...Ts>
struct signature_t<R(*)(Ts...)> : signature_t<R(Ts...)> {};

#define NUPACK_TMP(C, Q, C2) \
template <class R, class C, class ...Ts> \
    struct signature_t<R (C::* )(Ts...) Q> : pack<R, C2, Ts...> { \
        using return_type = R; \
        using size = size_constant<sizeof...(Ts) + 2>; \
        using pack_type = pack<R, C2, Ts...>; \
        template <class F> static void for_each_arg(F &&f) {NUPACK_UNPACK(f(type_t<Ts>()));} \
        template <class F> static void while_each_arg(F &&f) {NUPACK_WHILE(f(type_t<Ts>()));}};

    NUPACK_TMP(C, , C &);
    NUPACK_TMP(C, const, C const &);
    NUPACK_TMP(C, &, C &);
    NUPACK_TMP(C, const &, C const &);
    NUPACK_TMP(C, &&, C &&);
    NUPACK_TMP(C, const &&, C const &&);
#undef NUPACK_TMP

template <class T>
struct functor_signature;

#define NUPACK_TMP(Q) template <class R, class C, class ...Ts> \
    struct functor_signature<R (C::* )(Ts...) Q> : signature_t<R(Ts...)> {};
    NUPACK_TMP( );
    NUPACK_TMP(&);
    NUPACK_TMP(&&);
    NUPACK_TMP(const);
    NUPACK_TMP(const &);
    NUPACK_TMP(const &&);
#undef NUPACK_TMP

/******************************************************************************************/

template <class F>
struct signature_t<F, void_t<decltype(&F::operator())>> : functor_signature<decltype(&F::operator())> {};

template <class F> signature_t<F> signature() {return {};}
template <class F> signature_t<F> signature(F const &) {return {};}

/******************************************************************************************/

template <class C, class ...Ts>
auto call_operator() {
    using R = decltype(declref<C>()(declval<Ts>()...));
    using F = if_t<is_const<C>, R(decay<C>::*)(Ts...) const, R(decay<C>::*)(Ts...)>;
    return static_cast<F>(&decay<C>::operator());
}

/******************************************************************************************/

}
