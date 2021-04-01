/**
 * @brief Common metaprogramming functions and traits
 *
 * @file TypeSupport.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#define BOOST_NO_IOSTREAM

#include <tuple>
#include <string>
#include <stdexcept>
#include <utility>
#include <climits>
#include <type_traits>
#include <functional>

#define NUPACK_IF(...) typename std::enable_if<__VA_ARGS__, bool>::type = 0

namespace nupack {

/******************************************************************************************/

template <class T> T declval();
template <class T> T & declref();
template <class T> static constexpr T * null_of = static_cast<T *>(nullptr);

/******************************************************************************************/

namespace detail {template <class ...Ts> struct make_void {using type = void;};}
/// void_t, takes idea from cppreference.com
template <class ...Ts> using void_t = typename detail::make_void<Ts...>::type;

/******************************************************************************************/

template <std::size_t ...Is>          using indices_t = std::index_sequence<Is...>;
template <class T>                    using no_ref  = std::remove_reference_t<T>;
template <class T>                    using no_cv   = std::remove_cv_t<T>;
template <class T>                    using no_qual = std::remove_cv_t<std::remove_reference_t<T>>;
template <bool B, class T1, class T2> using if_t    = std::conditional_t<B, T1, T2>;
template <class T>                    using decay   = std::decay_t<T>;
template <bool B>                     using int_if  = std::enable_if_t<B, int>;
template <bool B>                     using void_if = std::enable_if_t<B, void>;

/******************************************************************************************/

struct base_type_t {};

template <class T> struct type_t : base_type_t {
    using type = T;
    T operator*() const;
};

template <class T> static constexpr type_t<T> type_c{};

/******************************************************************************************/

namespace detail {
    template <class ...Ts> struct type_at_t {template <class T, class ...Us> static T get(Ts..., T t, Us...) {return t;}};

    template <std::size_t ...Is, class ...Ts>
    auto type_at_2(indices_t<Is...>, Ts...ts) {return type_at_t<decltype(Is, base_type_t())...>::get(ts...);}

    template <std::size_t I, class ...Ts>
    auto type_at(Ts...ts) {return type_at_2(std::make_index_sequence<I>(), ts...);}
}

template <std::size_t I, class ...Ts> using type_at = typename decltype(detail::type_at<I>(type_t<Ts>()...))::type;

/// Alias for std::forward()
template <class T> constexpr T && fw(no_ref<T> &t)  noexcept {return static_cast<T &&>(t);}
/// Alias for std::forward()
template <class T> constexpr T && fw(no_ref<T> &&t) noexcept {return static_cast<T &&>(t);}

/******************************************************************************************/

/// Compile time if for two expressions
template <bool B, class T, class U, NUPACK_IF(B)>
decltype(auto) if_c(T &&t, U &&) {return fw<T>(t);}

template <bool B, class T, class U, NUPACK_IF(!B)>
decltype(auto) if_c(T &&, U &&u) {return fw<U>(u);}

/******************************************************************************************/

template <class B, class T, class U, NUPACK_IF(B::value)>
constexpr T if_c(B const &, T &&t, U &&) {return fw<T>(t);}

template <class B, class T, class U, NUPACK_IF(!B::value)>
constexpr U if_c(B const &, T &&, U &&u) {return fw<U>(u);}

/******************************************************************************************/

/// Compile time if which passes a parameter into either expression
template <bool B, class T, class U, class X, NUPACK_IF(B)>
decltype(auto) if_c(X &&x, T &&t, U &&) {return fw<T>(t)(fw<X>(x));}
template <bool B, class T, class U, class X, NUPACK_IF(!B)>
decltype(auto) if_c(X &&x, T &&, U &&u) {return fw<U>(u)(fw<X>(x));}

/******************************************************************************************/

template <bool B, class T, class F, NUPACK_IF(!B)>
void eval_if(T const &, F const &) {}
template <bool B, class T, class F, NUPACK_IF(B)>
void eval_if(T &&t, F &&f) {fw<F>(f)(fw<T>(t));}

template <class B, class F, class G, NUPACK_IF(bool(B::value))>
auto eval_one(B b, F &&f, G &&) {return fw<F>(f)(b);}

template <class B, class F, class G, NUPACK_IF(!bool(B::value))>
auto eval_one(B b, F &&, G &&g) {return fw<G>(g)(b);}

/******************************************************************************************/

using std::swap;

/******************************************************************************************/

template <class T, class U> decltype(auto) sink(U &&u) {return fw<U>(u);}
template <class T, class ...Us> using sink_type = T;
template <class T, class ...Us> using sink_sfinae = if_t<true, T, void_t<Us...>>;

/******************************************************************************************/

template <class T, class S> auto construct_ptr(S const &s) -> decltype(T(s))  {return T(s);}
template <class T, class S> auto construct_ptr(S const &s) -> decltype(T(&s)) {return T(&s);}

/******************************************************************************************/

template <class T> T & remove_const(T const &t) {return const_cast<T &>(t);}
template <class T> T const & add_const(T &t) {return t;}
template <class T> auto copy(T &&t) {return fw<T>(t);}

/******************************************************************************************/

template <bool ...Bs> struct bool_pack;
template <bool ...Bs> static constexpr bool all_of_c  = std::is_same<bool_pack<true, Bs...>, bool_pack<Bs..., true>>::value;
template <bool ...Bs> static constexpr bool none_of_c = std::is_same<bool_pack<false, Bs...>, bool_pack<Bs..., false>>::value;
template <bool ...Bs> static constexpr bool any_of_c  = !none_of_c<Bs...>;

/******************************************************************************************/

template <class F, NUPACK_IF(std::is_member_pointer<F>::value)>
decltype(auto) to_functor(F f) {return [f](auto &&i) {return i.*f;};}

template <class F, NUPACK_IF(!std::is_member_pointer<no_ref<F>>::value)>
decltype(auto) to_functor(F &&f) {return fw<F>(f);}

template <std::size_t ...Is, class F>
auto call_with(F &&f) {return [call=fw<F>(f)](auto &&...is) {
    auto t = std::forward_as_tuple(is...); return call(std::get<Is>(t)...);
};}

/******************************************************************************************/

template <class T> struct constructor {
    template <class ...Ts> constexpr T operator()(Ts &&...ts) const {return T{fw<Ts>(ts)...};}
};

/******************************************************************************************/

template <class T, NUPACK_IF( std::is_pointer<T>::value)>
decltype(auto) arrow(T t) {return t;}

template <class T, NUPACK_IF(!std::is_pointer<T>::value)>
decltype(auto) arrow(T const &t) {return t.operator->();}

/******************************************************************************************/

/// Class that inherits from all template types and uses constructor arguments for the first type
template <class Base, class ...Bases> struct Mixer : Base, Bases... {
    template <class ...Ts> explicit Mixer(Ts &&...ts) : Base(fw<Ts>(ts)...) {}
};

/******************************************************************************************/

template <class To, class From, NUPACK_IF(std::is_convertible<From, To>::value)>
To convert(From &&from) {return To{std::forward<From>(from)};}

namespace types {
    template <class U> void index_of_type();

    template <class U, class T, class ...Ts, NUPACK_IF(std::is_same<U, T>::value)>
    std::integral_constant<std::size_t, 0> index_of_type();

    template <class U, class T, class ...Ts, NUPACK_IF(!std::is_same<U, T>::value)>
    auto index_of_type() {return std::integral_constant<std::size_t, decltype(index_of_type<U, Ts...>())::value + 1>();}
}

template <std::size_t N, class T>
constexpr auto at(T &&t) -> decltype(std::get<N>(std::forward<T>(t))) {return std::get<N>(std::forward<T>(t));}

template <class ...Ts> static constexpr auto index_of_type = decltype(types::index_of_type<Ts...>())::value;

static_assert(index_of_type<int, int> == 0, "");
static_assert(index_of_type<double, int, double> == 1, "");
static_assert(index_of_type<double, double, int> == 0, "");
static_assert(index_of_type<double, int, bool, double> == 2, "");

/******************************************************************************************/

namespace detail {
    template <class T, class ...Ts> struct match_t {
        template <class ...Us> constexpr auto operator()(Us ...us) const {
            return std::get<index_of_type<T, Ts...>>(std::make_tuple(us...));
        }
    };

    template <class From, class To> struct copy_qualifier_t {using type = decay<To>;};
    template <class From, class To> struct copy_qualifier_t<From &, To> {using type = decay<To> &;};
    template <class From, class To> struct copy_qualifier_t<From const &, To> {using type = decay<To> const &;};
    template <class From, class To> struct copy_qualifier_t<From const, To> {using type = decay<To> const;};

    template <class T> struct no_mutable_t {using type = T const;};
    template <class T> struct no_mutable_t<T &> {using type = T const &;};
    template <class T> struct no_mutable_t<T &&> {using type = T const &&;};
}

/******************************************************************************************/

template <class From, class To> using copy_qualifier = typename detail::copy_qualifier_t<From, To>::type;

template <class ...Ts> static constexpr auto match = detail::match_t<Ts...>();

template <class T> using no_mutable = typename detail::no_mutable_t<T>::type;

/******************************************************************************************/

template <class T> auto lval_addressof(T const &t) {return std::addressof(t);}
template <class T> auto lval_addressof(T&&) {return nullptr;}

/******************************************************************************************/

template <class T> struct caster_t {
    template <class U>
    T operator() (U &&u) const {return static_cast<T>(fw<U>(u));}
    template <class ...Us, NUPACK_IF(1 != sizeof...(Us))>
    T operator() (Us &&...us) const {return T(fw<Us>(us)...);}
};

template <class T> static constexpr auto caster = caster_t<T>{};

/******************************************************************************************/

template <std::size_t N, std::size_t I, std::size_t ...Is, NUPACK_IF(N == 0)>
std::integral_constant<std::size_t, I> index_at(indices_t<I, Is...>) {return {};}

template <std::size_t N, std::size_t I, std::size_t ...Is, NUPACK_IF(N != 0)>
auto index_at(indices_t<I, Is...>) {return index_at<N-1>(indices_t<Is...>());}

template <std::size_t N>
struct at_t {
    template <class T>
    constexpr auto operator()(T &&t) const -> decltype(std::get<N>(fw<T>(t))) {return std::get<N>(fw<T>(t));}

    template <std::size_t ...Is>
    constexpr auto operator()(indices_t<Is...> t) const {
        static_assert(N < sizeof...(Is), "Compile time index out of range");
        return index_at<N>(t);
    }
};

template <class N, class T>
constexpr decltype(auto) at_c(T &&t, N={}) {return at_t<N::value>()(fw<T>(t));}

template <std::size_t N, class T>
constexpr decltype(auto) at_c(T &&t) {return at_t<N>()(fw<T>(t));}

template <class N>
constexpr at_t<N::value> at_c(N={}) {return {};}

template <std::size_t N>
constexpr at_t<N> at_c() {return {};}

static constexpr auto first_of  = at_t<0>();
static constexpr auto second_of = at_t<1>();
static constexpr auto third_of  = at_t<2>();
static constexpr auto fourth_of = at_t<3>();
static constexpr auto fifth_of  = at_t<4>();

/******************************************************************************************/

template <class T> constexpr auto bitsof() {return CHAR_BIT * sizeof(T);}
template <class T> constexpr auto bitsof(T const &) {return CHAR_BIT * sizeof(T);}

/******************************************************************************************/

struct Ignore {
    template <class ...Ts>
    Ignore(Ts const &...ts) {}
};

}
