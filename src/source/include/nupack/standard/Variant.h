#pragma once
#include "../reflect/Memory.h"
#include "../reflect/Print.h"
#include <variant>

namespace nupack {

/**
 * @brief Variant is std::variant
 * @tparam Ts possible types
 * @todo Switch to std::variant
 */
template <class ...Ts> using Variant = std::variant<Ts...>;

NUPACK_DEFINE_VARIADIC(is_variant, Variant, class);

/// Get the types out of a variant
template <class ...Ts>
struct detail::as_pack_t<Variant<Ts...>> {using type = pack<Ts...>;};

/// Like std::visit but the function is last
template <class V, class ...Vs, NUPACK_IF(is_variant<decay<V>>)>
decltype(auto) dispatch(V &&v, Vs &&...vs) {
    auto ts = std::forward_as_tuple(v, vs...);
    auto &&f = std::get<sizeof...(Vs)>(ts);
    return unpack(slice<0, sizeof...(Vs)>(ts), [&](auto &&...as) {
        return std::visit(fw<decltype(f)>(f), fw<decltype(as)>(as)...);
    });
}

/******************************************************************************************/

namespace detail {
    template <class A, class T>
    A get_variant(std::size_t n) {if (n) throw std::out_of_range("any out of range"); return T();}

    template <class A, class T, class U, class ...Ts>
    A get_variant(std::size_t n) {return n ? get_variant<A, U, Ts...>(n - 1) : T();}
}

/// Helper to create a variant from an enum/integer
template <class T>
struct VariantEnum {
    T index;
    NUPACK_REFLECT(VariantEnum, index);
    constexpr VariantEnum(T i={}) : index(i) {}

    template <class ...Ts>
    constexpr VariantEnum(Variant<Ts...> const &a) : index(a.which()) {}

    template <class ...Ts>
    operator Variant<Ts...>() const {return detail::get_variant<Variant<Ts...>, Ts...>(static_cast<std::size_t>(index));}
};

/******************************************************************************************/

/// std::visit for variant, only one variant argument
template <class V, class F, NUPACK_IF(is_variant<decay<V>>)>
decltype(auto) fork(V &&v, F &&f) {return std::visit(fw<F>(f), fw<V>(v));}

/// overload for std::visit in case argument is not a variant
template <class V, class F, NUPACK_IF(!is_variant<decay<V>>)>
decltype(auto) fork(V &&v, F &&f) {return fw<F>(f)(fw<V>(v));}

/******************************************************************************************/

namespace detail {
    template <class ...Ts> struct MaybeVariant {using type = Variant<Ts...>;};
    template <class T> struct MaybeVariant<T> {using type = T;};
}

template <class ...Ts> using MaybeVariant = typename detail::MaybeVariant<Ts...>::type;

template <class T, class V, NUPACK_IF(is_variant<decay<V>>)>
auto maybe_get(V &&v) {return std::get_if<T>(std::addressof(fw<V>(v)));}

template <class T, class V, NUPACK_IF(!is_variant<decay<V>>)>
auto maybe_get(V &&v) {
    static_assert(is_same<decay<V>, decay<T>>, "types do not agree");
    return std::addressof(fw<V>(v));
}

/******************************************************************************************/

/// Take a variant and return a variant based on the possible return types
template <class F, class ...Vs, class ...Ts>
auto map_variant(Variant<Vs...> const &v, F &&f, Ts &&...ts) {
    using R = Variant<decltype(fw<F>(f)(declval<Vs const &>(), declval<Ts &&>()...))...>;
    return std::visit([&](auto const &i) -> R {return fw<F>(f)(i, fw<Ts>(ts)...);}, v);
}

/// Take a variant and return a variant based on the possible return types
template <class F, class ...Vs, class ...Ts>
auto map_variant(Variant<Vs...> &v, F &&f, Ts &&...ts) {
    using R = Variant<decltype(fw<F>(f)(declval<Vs &>(), declval<Ts &&>()...))...>;
    return std::visit([&](auto &i) -> R {return fw<F>(f)(i, fw<Ts>(ts)...);}, v);
}

/// Take a variant and return a variant based on the possible return types
template <class F, class ...Vs, class ...Ts>
auto map_variant(Variant<Vs...> &&v, F &&f, Ts &&...ts) {
    using R = Variant<decltype(fw<F>(f)(declval<Vs &&>(), declval<Ts &&>()...))...>;
    return std::visit([&](auto &&i) -> R {return fw<F>(f)(std::move(i), fw<Ts>(ts)...);}, std::move(v));
}

/******************************************************************************************/

template <class ...Ts> struct memory::impl<Variant<Ts...>, void> {
    std::size_t operator()(Variant<Ts...> const &t) const {
        return sizeof(t) + dispatch(t, [](auto const &i) {return measure(i);});
    }
    void erase(Variant<Ts...> &) const {}
};

template <class ...Ts> struct io::Printer<Variant<Ts...>, void> {
    void operator()(std::ostream &os, Variant<Ts...> const &t, io::Indent id={}) const {
        dispatch(t, [&](auto const &i) {Printer<no_qual<decltype(i)>>()(os, i, id);});
    }
};

/******************************************************************************************/

}
