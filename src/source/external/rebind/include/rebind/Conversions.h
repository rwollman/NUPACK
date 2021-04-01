#pragma once
#include "Variable.h"

namespace rebind {

struct Default {};
struct Specialized {};
struct ADL {};

/******************************************************************************/

template <class T, class=void>
struct ImplicitConversions {
    using types = Pack<>;
};

template <class U, class T>
bool implicit_match(Variable &out, Type<U>, TargetQualifier const q, T &&t) {
    DUMP("implicit_match: ", typeid(U).name(), typeid(Type<T &&>).name(), q);
    if constexpr(std::is_convertible_v<T &&, U>)
        if (q == Value) out = {Type<U>(), static_cast<T &&>(t)};
    if constexpr(std::is_convertible_v<T &&, U &>)
        if (q == Lvalue) out = {Type<U &>(), static_cast<T &&>(t)};
    if constexpr(std::is_convertible_v<T &&, U &&>)
        if (q == Rvalue) out = {Type<U &&>(), static_cast<T &&>(t)};
    if constexpr(std::is_convertible_v<T &&, U const &>)
        if (q == Const) out = {Type<U const &>(), static_cast<T &&>(t)};
    DUMP("implicit_response result: ", out.has_value(), typeid(Type<T &&>).name(), typeid(U).name(), q);
    return out.has_value();
}

template <class U, class T>
bool recurse_implicit(Variable &out, Type<U>, TypeIndex const &idx, TargetQualifier q, T &&t);

/******************************************************************************/

template <class T>
bool implicit_response(Variable &out, TypeIndex const &idx, TargetQualifier q, T &&t) {
    DUMP("implicit_response: ", typeid(Type<T &&>).name(), idx.name(), typeid(typename ImplicitConversions<std::decay_t<T>>::types).name(), q);
    return ImplicitConversions<std::decay_t<T>>::types::apply([&](auto ...ts) {
        static_assert((!decltype(is_same(+Type<T>(), +ts))::value && ...), "Implicit conversion creates a cycle");
        return ((idx.matches(ts) && implicit_match(out, ts, q, static_cast<T &&>(t))) || ...)
            || (recurse_implicit(out, ts, idx, q, static_cast<T &&>(t)) || ...);
    });
}

template <class U, class T>
bool recurse_implicit(Variable &out, Type<U>, TypeIndex const &idx, TargetQualifier q, T &&t) {
    if constexpr(std::is_convertible_v<T &&, U &&>)
        return implicit_response(out, idx, q, static_cast<U &&>(t));
    else if constexpr(std::is_convertible_v<T &&, U &>)
        return implicit_response(out, idx, q, static_cast<U &>(t));
    else if constexpr(std::is_convertible_v<T &&, U const &>)
        return implicit_response(out, idx, q, static_cast<U const &>(t));
    else if constexpr(std::is_convertible_v<T &&, U>)
        return implicit_response(out, idx, q, static_cast<U>(t));
    return false;
}

/******************************************************************************/

/// Default response just tries implicit conversions
template <class T, TargetQualifier Q, class SFINAE>
struct Response {
    using method = Default;
    static_assert(std::is_same_v<unqualified<T>, T>);

    template <class T2>
    bool operator()(Variable &out, TypeIndex const &idx, T2 &&t) const {
        DUMP("no conversion found from source ", TypeIndex(typeid(T), Q), " to ", idx);
        return implicit_response(out, idx, Q, static_cast<T2 &&>(t));
    }
};

/******************************************************************************/

template <class T, class=void>
struct ResponseMethod {using type = Specialized;};

template <class T>
struct ResponseMethod<T, std::void_t<typename Response<T>::method>> {using type = typename Response<T>::method;};

template <class T>
using response_method = typename ResponseMethod<T>::type;

/******************************************************************************/

template <class T, TargetQualifier Q>
struct Response<T, Q, std::void_t<decltype(response(std::declval<TypeIndex>(), std::declval<T>()))>> {
    using method = ADL;
    static_assert(std::is_same_v<unqualified<T>, T>);

    template <class T2>
    bool operator()(Variable &out, TypeIndex const &idx, T2 &&t) const {
        DUMP("ADL Response: ", typeid(T), idx);
        out = response(idx, static_cast<T2 &&>(t));
        return out || implicit_response(out, idx, Q, static_cast<T2 &&>(t));
    }
};

/******************************************************************************/

void lvalue_fails(Variable const &, Dispatch &, TypeIndex);
void rvalue_fails(Variable const &, Dispatch &, TypeIndex);

/// The default behavior has no custom conversions
template <class T, class SFINAE>
struct Request {
    static_assert(std::is_same_v<unqualified<T>, T>);
    using method = Default;

    std::optional<T> operator()(Variable const &r, Dispatch &msg) const {
        return msg.error("mismatched class type", typeid(T));
    }
};

/// The default behavior has no custom conversions
template <class T, class SFINAE>
struct Request<T &, SFINAE> {
    using method = Default;

    T * operator()(Variable const &v, Dispatch &msg) const {
        lvalue_fails(v, msg, typeid(T));
        return nullptr;
    }
};

/// The default behavior tries the T -> T const & and T & -> T const & routes
template <class T, class SFINAE>
struct Request<T const &, SFINAE> {
    using method = Default;

    T const * operator()(Variable const &v, Dispatch &msg) const {
        DUMP("trying & -> const & ", typeid(T).name());
        if (auto p = v.request<T &>(msg)) return p;
        DUMP("trying temporary const & storage ", typeid(T).name());
        if (auto p = v.request<T>(msg)) return msg.store(std::move(*p));
        return msg.error("could not bind to const lvalue reference", typeid(T)), nullptr;
    }
};

/// The default behavior tries the T -> T && route
template <class T, class SFINAE>
struct Request<T &&, SFINAE> {
    using method = Default;

    T * operator()(Variable const &v, Dispatch &msg) const {
        DUMP("trying temporary && storage ", typeid(T).name());
        if (auto p = v.request<T>(msg)) return msg.store(std::move(*p));
        rvalue_fails(v, msg, typeid(T));
        return nullptr;
    }
};

/******************************************************************************/

template <class T, class=void>
struct RequestMethod {using method = Specialized;};

template <class T>
struct RequestMethod<T, std::void_t<typename Request<T>::method>> {using type = typename Request<T>::method;};

template <class T>
using request_method = typename RequestMethod<T>::type;

/******************************************************************************/

void request(int, int, int);

/// ADL version of Request
template <class T>
struct Request<T, std::void_t<decltype(request(Type<T>(), std::declval<Variable const &>(), std::declval<Dispatch &>()))>> {
    using method = ADL;

    std::optional<T> operator()(Variable const &r, Dispatch &msg) const {
        return static_cast<std::optional<T>>(request(Type<T>(), r, msg));
    }
};

/******************************************************************************/

}
