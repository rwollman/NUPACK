#pragma once
#include "Document.h"
#include "BasicTypes.h"
#include <tuple>
#include <utility>
#include <array>
#include <optional>
#include <variant>
#include <map>

namespace rebind {

/******************************************************************************/

template <class T, Qualifier Q>
struct Response<std::optional<T>, Q> {
    bool operator()(Variable &out, TypeIndex t, std::optional<T> &v) const {
        return v ? get_response<Q>(out, std::move(t), *v) : false;
    }
    bool operator()(Variable &out, TypeIndex t, std::optional<T> &&v) const {
        return v ? get_response<Q>(out, std::move(t), std::move(*v)) : false;
    }
    bool operator()(Variable &out, TypeIndex t, std::optional<T> const &v) const {
        return v ? get_response<Q>(out, std::move(t), *v) : false;
    }
};

template <class T>
struct Request<std::optional<T>> {
    std::optional<std::optional<T>> operator()(Variable const &v, Dispatch &msg) const {
        std::optional<std::optional<T>> out;
        if (!v || v.request<std::nullptr_t>()) out.emplace();
        else if (auto p = v.request<std::remove_cv_t<T>>(msg))
            out.emplace(std::move(*p));
        return out;
    }
};

/******************************************************************************/

template <class T, TargetQualifier Q>
struct Response<std::shared_ptr<T>, Q> {
    bool operator()(Variable &out, TypeIndex t, std::shared_ptr<T> const &p) const {
        DUMP("shared_ptr", t, type_index<T>(), bool(p), Q);
        return p && get_response<Q>(out, std::move(t), *p);
    }
};

template <class T>
struct Request<std::shared_ptr<T>> {
    std::optional<std::shared_ptr<T>> operator()(Variable const &v, Dispatch &msg) const {
        std::optional<std::shared_ptr<T>> out;
        if (!v || v.request<std::nullptr_t>()) out.emplace();
        else if (auto p = v.request<std::remove_cv_t<T>>(msg))
            out.emplace(std::make_shared<T>(std::move(*p)));
        return out;
    }
};

/******************************************************************************/

template <TargetQualifier Q, class ...Ts>
struct Response<std::variant<Ts...>, Q> {
    template <class V>
    bool operator()(Variable &out, TypeIndex t, V &&v) const {
        return std::visit([&](auto &&x) {return get_response<Q>(out, t, static_cast<decltype(x) &&>(x));}, static_cast<V &&>(v));
    }
};

static_assert(std::is_same_v<response_method<std::variant<int>>, Specialized>);

template <class ...Ts>
struct Request<std::variant<Ts...>> {
    template <class T>
    static bool put(std::optional<std::variant<Ts...>> &out, Variable const &v, Dispatch &msg) {
        if (auto p = v.request<T>(msg)) return out.emplace(std::move(*p)), true;
        return false;
    }

    std::optional<std::variant<Ts...>> operator()(Variable const &v, Dispatch &msg) const {
        std::optional<std::variant<Ts...>> out;
        (void) (put<Ts>(out, v, msg) || ...);
        return out;
    }
};

/******************************************************************************/

template <class V>
struct MapResponse {
    using T = std::pair<typename V::key_type, typename V::mapped_type>;

    bool operator()(Variable &o, TypeIndex t, V &&v) const {
        return range_response<T>(o, t, std::make_move_iterator(std::begin(v)), std::make_move_iterator(std::end(v)));
    }

    bool operator()(Variable &o, TypeIndex t, V const &v) const {
        return range_response<T>(o, t, std::begin(v), std::end(v));
    }
};

template <class V>
struct MapRequest {
    using T = std::pair<typename V::key_type, typename V::mapped_type>;

    std::optional<V> operator()(Variable const &v, Dispatch &msg) const {
        std::optional<V> out;
        if (auto p = v.request<Vector<T>>())
            out.emplace(std::make_move_iterator(std::begin(*p)), std::make_move_iterator(std::end(*p)));
        return out;
    }
};

template <class K, class V, class C, class A>
struct Request<std::map<K, V, C, A>> : MapRequest<std::map<K, V, C, A>> {};

template <class K, class V, class C, class A>
struct Response<std::map<K, V, C, A>> : MapResponse<std::map<K, V, C, A>> {};

/******************************************************************************/

template <class F>
struct FunctionRequest {
    std::optional<F> operator()(Variable const &v, Dispatch &msg) const {
        if (auto p = v.request<Callback<typename F::result_type>>(msg)) return F{std::move(*p)};
        return {};
    }
};

template <class R, class ...Ts>
struct Request<std::function<R(Ts...)>> : FunctionRequest<std::function<R(Ts...)>> {};

/******************************************************************************/

}
