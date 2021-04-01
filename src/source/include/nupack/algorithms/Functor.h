/**
 * @brief Common functors
 *
 * @file Functor.h
 * @author Mark Fornace
 * @date 2018-05-18
 */
#pragma once
#include "Traits.h"

// Static lambda (bottom) attributed to http://pfultz2.com/blog/2014/09/02/static-lambda/

namespace nupack {

NUPACK_BINARY_FUNCTOR(min,         (static_cast<std::common_type_t<T, U>>(u < t ? u : t)));
NUPACK_BINARY_FUNCTOR(min_eq,      (u < t ? (t = u) : u));
NUPACK_BINARY_FUNCTOR(max,         (static_cast<std::common_type_t<T, U>>(u < t ? t : u)));
NUPACK_BINARY_FUNCTOR(max_eq,      (u > t ? (t = u) : u));
NUPACK_BINARY_FUNCTOR(assign_eq,   (fw<T>(t) = fw<U>(u)));

NUPACK_BINARY_FUNCTOR(plus,        fw<T>(t) + fw<U>(u));
NUPACK_BINARY_FUNCTOR(minus,       fw<T>(t) - fw<U>(u));
NUPACK_BINARY_FUNCTOR(times,       fw<T>(t) * fw<U>(u));
NUPACK_BINARY_FUNCTOR(divide,      fw<T>(t) / fw<U>(u));
NUPACK_BINARY_FUNCTOR(modulus,     fw<T>(t) % fw<U>(u));
NUPACK_BINARY_FUNCTOR(bitwise_or,  fw<T>(t) | fw<U>(u));
NUPACK_BINARY_FUNCTOR(bitwise_and, fw<T>(t) & fw<U>(u));
NUPACK_BINARY_FUNCTOR(bitwise_xor, fw<T>(t) ^ fw<U>(u));
NUPACK_BINARY_FUNCTOR(logical_and, fw<T>(t) && fw<U>(u));
NUPACK_BINARY_FUNCTOR(logical_or,  fw<T>(t) || fw<U>(u));
NUPACK_BINARY_FUNCTOR(equals,      fw<T>(t) == fw<U>(u));
NUPACK_BINARY_FUNCTOR(not_equals,  fw<T>(t) != fw<U>(u));
NUPACK_BINARY_FUNCTOR(less,        fw<T>(t) < fw<U>(u));
NUPACK_BINARY_FUNCTOR(greater,     fw<T>(t) > fw<U>(u));
NUPACK_BINARY_FUNCTOR(lshift,      fw<T>(t) << fw<U>(u));
NUPACK_BINARY_FUNCTOR(rshift,      fw<T>(t) >> fw<U>(u));
NUPACK_BINARY_FUNCTOR(comma,       (fw<T>(t), fw<U>(u)));
NUPACK_BINARY_FUNCTOR(subscript,   fw<T>(t)[fw<U>(u)]);

NUPACK_BINARY_FUNCTOR(less_eq,     fw<T>(t) <= fw<U>(u));
NUPACK_BINARY_FUNCTOR(greater_eq,  fw<T>(t) >= fw<U>(u));
NUPACK_BINARY_FUNCTOR(assign_to,   fw<T>(t) = fw<U>(u));
NUPACK_BINARY_FUNCTOR(plus_eq,     fw<T>(t) += fw<U>(u));
NUPACK_BINARY_FUNCTOR(minus_eq,    fw<T>(t) -= fw<U>(u));
NUPACK_BINARY_FUNCTOR(times_eq,    fw<T>(t) *= fw<U>(u));
NUPACK_BINARY_FUNCTOR(divide_eq,   fw<T>(t) /= fw<U>(u));
NUPACK_BINARY_FUNCTOR(modulus_eq,  fw<T>(t) %= fw<U>(u));
NUPACK_BINARY_FUNCTOR(lshift_eq,   fw<T>(t) <<= fw<U>(u));
NUPACK_BINARY_FUNCTOR(rshift_eq,   fw<T>(t) >>= fw<U>(u));

NUPACK_UNARY_FUNCTOR(dereference,    *fw<T>(t));
NUPACK_UNARY_FUNCTOR(increment,      ++fw<T>(t));
NUPACK_UNARY_FUNCTOR(decrement,      --fw<T>(t));
NUPACK_UNARY_FUNCTOR(post_increment, fw<T>(t)++);
NUPACK_UNARY_FUNCTOR(post_decrement, fw<T>(t)--);
NUPACK_UNARY_FUNCTOR(logical_not,    !fw<T>(t));
NUPACK_UNARY_FUNCTOR(bitwise_not,    ~fw<T>(t));
NUPACK_UNARY_FUNCTOR(unary_plus,     +fw<T>(t));
NUPACK_UNARY_FUNCTOR(unary_minus,    -fw<T>(t));
NUPACK_UNARY_FUNCTOR(address,        &fw<T>(t));
NUPACK_UNARY_FUNCTOR(address_of,     std::addressof(fw<T>(t)));
NUPACK_UNARY_FUNCTOR(iter_ptr,       std::addressof(*fw<T>(t)));
NUPACK_UNARY_FUNCTOR(size_of,        sizeof(t));
NUPACK_UNARY_FUNCTOR(throw_op,       (void)(throw fw<T>(t)));
NUPACK_UNARY_FUNCTOR(copy_op,        T(static_cast<decay<T> const &>(fw<T>(t))));
NUPACK_UNARY_FUNCTOR(move_op,        std::move(t));
NUPACK_UNARY_FUNCTOR(delete_op,      (void) delete t);

NUPACK_UNARY_FUNCTOR(is_alpha,       std::isalpha(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(is_alnum,       std::isalnum(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(is_lower,       std::islower(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(is_upper,       std::isupper(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(is_digit,       std::isdigit(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(is_space,       std::isspace(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(to_lower,       std::tolower(static_cast<unsigned char>(t)));
NUPACK_UNARY_FUNCTOR(to_upper,       std::toupper(static_cast<unsigned char>(t)));
NUPACK_BINARY_FUNCTOR(less_abs,      std::abs(fw<T>(t)) < std::abs(fw<U>(u)));
NUPACK_BINARY_FUNCTOR(greater_abs,      std::abs(fw<T>(t)) > std::abs(fw<U>(u)));

/******************************************************************************************/

/// Composition of two functions F(G(x))
template <class F, class G>
struct Composition {
    F f;
    G g;
    template <class ...Ts>
    auto operator()(Ts &&...ts) -> decltype(f(g(fw<Ts>(ts)...))) const {return f(g(fw<Ts>(ts)...));}
};

template <class F, class G>
constexpr auto compose(F &&f, G &&g) {return Composition<no_qual<F>, no_qual<G>>{fw<F>(f), fw<G>(g)};}

/******************************************************************************************/

/// call the first argument with the rest of the arguments
struct call_t {
    template <class F, class ...Ts>
    constexpr auto operator()(F &&f, Ts &&...ts) const -> decltype(fw<F>(f)(fw<Ts>(ts)...)) {return fw<F>(f)(fw<Ts>(ts)...);}
};
static constexpr auto call = call_t{};

/******************************************************************************************/

/// f(x) = x+1; it might be more optimized or more flexible depending on the situation
struct plus_one_t {
    template <class T, NUPACK_IF(can_increment<T>)>
    constexpr T operator()(T t) const {return ++t;}

    template <class T, NUPACK_IF(!can_increment<T> && can_plus_one<T>)>
    constexpr T operator()(T const &t) const {return t + 1u;}

    template <class T, NUPACK_IF(!can_increment<T> && !can_plus_one<T>)>
    constexpr T operator()(T const &t) const {return std::next(t);}
};

static constexpr auto plus_one = plus_one_t();

/// f(x) = x-1; it might be more optimized or more flexible depending on the situation
struct minus_one_t {
    template <class T, NUPACK_IF(can_decrement<T>)>
    constexpr T operator()(T t) const {return --t;}

    template <class T, NUPACK_IF(!can_increment<T> && can_minus_one<T>)>
    constexpr T operator()(T const &t) const {return t - 1u;}

    template <class T, NUPACK_IF(!can_increment<T> && !can_minus_one<T>)>
    constexpr T operator()(T const &t) const {return std::prev(t);}
};

static constexpr auto minus_one = minus_one_t();

/******************************************************************************************/

/// iterator distance, also works for integers
struct distance_t {
    template <class T, NUPACK_IF(!is_integral<T>)>
    constexpr auto operator()(T const &t, T const &u) const {return std::distance(t, u);}
    template <class T, NUPACK_IF(is_integral<T>)>
    constexpr auto operator()(T const &t, T const &u) const {return u - t;}
};

static constexpr auto distance = distance_t();
NUPACK_BINARY_FUNCTOR(rdistance, distance(fw<U>(u), fw<T>(t)));

/******************************************************************************************/

/// Overloaded functors where they must be empty, then we don't store them
template <class ...Fs> class empty_overload_t {
    template <int I, class ...Ts, NUPACK_IF(I >= sizeof...(Fs))>
    auto call(Ts const &...) const {static_assert(indices_t<I, sizeof...(Fs)>::no_overload, "No (empty) overload found");}

    template <int I, class ...Ts, NUPACK_IF(I < sizeof...(Fs)), NUPACK_IF(can_call<type_at<I, Fs...> const &, Ts &&...>)>
    constexpr decltype(auto) call(Ts &&...ts) const {return reinterpret_cast<type_at<I, Fs...> const &>(*this)(fw<Ts>(ts)...);};

    template <int I, class ...Ts, NUPACK_IF(I < sizeof...(Fs)), NUPACK_IF(!can_call<type_at<I, Fs...> const &, Ts &&...>)>
    constexpr decltype(auto) call(Ts &&...ts) const {return call<I+1>(fw<Ts>(ts)...);}
public:
    template <class ...Ts> constexpr auto operator()(Ts &&...ts) const -> decltype(call<0>(fw<Ts>(ts)...)) {return call<0>(fw<Ts>(ts)...);};
};

template <class ...Fs, NUPACK_IF(all_of_c<is_empty<Fs>...>)>
constexpr auto overload(Fs ...) {return empty_overload_t<Fs...>{};}

/******************************************************************************************/

/// Overloaded functors where they are not empty, then we have to store them
template <class ...Fs> class overload_t {
    std::tuple<Fs...> overloads;

    template <int I, class ...Ts, NUPACK_IF(I >= sizeof...(Fs))>
    auto call(Ts const &...) const {static_assert(indices_t<I, sizeof...(Fs)>::no_overload, "No overload found");}

    template <int I, class ...Ts, NUPACK_IF(I < sizeof...(Fs)), NUPACK_IF(can_call<type_at<I, Fs...> const &, Ts &&...>)>
    constexpr decltype(auto) call(Ts &&...ts) const {return std::get<I>(overloads)(fw<Ts>(ts)...);};

    template <int I, class ...Ts, NUPACK_IF(I < sizeof...(Fs)), NUPACK_IF(!can_call<type_at<I, Fs...>  const &, Ts &&...>)>
    constexpr decltype(auto) call(Ts &&...ts) const {return call<I+1>(fw<Ts>(ts)...);}
public:
    template <class ...Ts> overload_t(Ts &&...ts) : overloads{fw<Ts>(ts)...} {}
    template <class ...Ts> constexpr auto operator()(Ts &&...ts) const -> decltype(call<0>(fw<Ts>(ts)...)) {return call<0>(fw<Ts>(ts)...);};
};

template <class ...Fs, NUPACK_IF(!all_of_c<is_empty<decay<Fs>>...>)>
constexpr auto overload(Fs &&...fs) {return overload_t<decay<Fs>...>{fw<Fs>(fs)...};}

/******************************************************************************************/

/// For a given functor, discard leftmost arguments until it is callable, then call it (callable case)
template <class F, class ...Ts, NUPACK_IF(can_call<F &&, Ts &&...>)>
decltype(auto) back_call(F &&f, Ts &&...ts) {return fw<F>(f)(fw<Ts>(ts)...);}
/// For a given functor, discard leftmost arguments until it is callable, then call it (noncallable case)
template <class F, class T, class ...Ts, NUPACK_IF(!can_call<F &&, T &&, Ts &&...>)>
decltype(auto) back_call(F &&f, T &&, Ts &&...ts) {return back_call(fw<F>(f), fw<Ts>(ts)...);}

/******************************************************************************************/

template <class F>
struct empty_lambda {
    static_assert(std::is_empty<F>::value, "Static lambdas must be empty");
    template <class ...Ts> constexpr auto operator()(Ts &&...ts) const -> decltype(declref<F const &>()(fw<Ts>(ts)...)) {
        return reinterpret_cast<F const &>(*this)(fw<Ts>(ts)...);
    };
};

struct lambda_eq {template <class T> no_ref<T> * operator=(T &&t) {return &t;}};

struct lambda_wrap {
    template <class F> constexpr empty_lambda<F> operator += (F *) {return {};}
};

template <class F>
constexpr auto function_ptr(F f) -> decltype(false ? nullptr : f) {return false ? nullptr : f;}

/******************************************************************************************/

}

/// static constexpr lambda for empty lambdas (no captures)
#define NUPACK_LAMBDA(name) \
static constexpr auto const name = ::nupack::lambda_wrap() += true ? nullptr : ::nupack::lambda_eq()
