#pragma once
#include <rebind/Document.h>
#include <armadillo>
#include <map>
#include <vector>
#include <variant>
#include <optional>
#include <nupack/reflect/Serialize.h>
#include <nupack/types/Matrix.h>

namespace rebind {

/******************************************************************************/

template <class V>
struct TupleRenderer {
    template <std::size_t I>
    static bool put(Variable &out, V const &v) {
        out.emplace(Type<decltype(std::get<I>(v))>(), std::get<I>(v));
        return true;
    }

    template <std::size_t ...Is>
    static void apply(Document &doc, std::index_sequence<Is...>) {
        doc.type(typeid(V), "std.Tuple");
        (doc.render<unqualified<std::tuple_element_t<Is, V>>>(), ...);
        doc.method(typeid(V), "[]", [](V &v, std::size_t i) {
            Variable out;
            if (i >= std::tuple_size_v<V>) throw std::out_of_range("index out of range");
            ((Is == i && put<Is>(out, v)) || ...);
            return out;
        });
        doc.method(typeid(V), "__len__", [](V const &) {return std::tuple_size_v<V>;});
    }
    void operator()(Document &doc) const {apply(doc, std::make_index_sequence<std::tuple_size_v<V>>());}
};

template <class V>
struct Renderer<V, std::enable_if_t<std::tuple_size<V>::value >= 0>> : TupleRenderer<V> {};

/******************************************************************************/

template <class It, class T>
struct Iter {
    It iter, end;
    void next() {if (iter != end) ++iter;}
    bool good() const {return iter != end;}
    T get() const {return iter != end ? *iter : throw std::out_of_range("invalid iterator");}
};

template <class It, class T>
struct Renderer<Iter<It, T>> {
    using I = Iter<It, T>;
    void operator()(Document &doc) const {
        doc.type(typeid(I), "std.Iterator");
        doc.method(typeid(I), "next", &I::next);
        doc.method(typeid(I), "good", &I::good);
        doc.method(typeid(I), "get", &I::get);
    }
};

/******************************************************************************/

template <class M>
struct MapRenderer {
    using K = typename M::key_type;
    using V = typename M::mapped_type;

    void operator()(Document &doc) const {
        using P = std::pair<typename M::key_type, typename M::mapped_type>;
        doc.render<P>();
        doc.type(typeid(M), "std.Map");
        doc.method(typeid(M), "__setitem__", [](M &m, K k, V p) {m.insert_or_assign(std::move(k), std::move(p));});
        doc.method(typeid(M), "[]", [](M &m, K const &t) -> decltype(m.at(t)) {return m.at(t);});
        doc.method(typeid(M), "__iter__", [](M &m) {return Iter<typename M::iterator, P>{m.begin(), m.end()};});
        doc.method(typeid(M), "__len__", [](M const &m) {return m.size();});
        doc.method(typeid(M), "value_type", [](M const &) {return std::type_index(typeid(typename M::value_type));});
        doc.method(typeid(M), "items", [](M const &m) {return std::vector<std::pair<K, V>>(std::begin(m), std::end(m));});
    }
};

/******************************************************************************/

template <class T, class C, class A>
struct Renderer<std::map<T, C, A>> : MapRenderer<std::map<T, C, A>> {};

/******************************************************************************/

template <class R, class ...Ts>
struct Renderer<std::function<R(Ts...)>> {
    using F = std::function<R(Ts...)>;
    void operator()(Document &doc) const {
        doc.render<unqualified<R>>();
        (doc.render<unqualified<Ts>>(), ...);
        doc.type(typeid(F), "std.Function");
        doc.method(typeid(F), "()", &F::operator());
        doc.method(typeid(F), "bool", &F::operator bool);
    }
};

/******************************************************************************/

template <class T>
struct Renderer<std::optional<T>> {
    void operator()(Document &doc) const {
        TypeIndex const t{typeid(std::optional<T>)};
        doc.type(t, "std.Optional");
        doc.method(t, "bool", [](std::optional<T> const &t) {return bool(t);});
        doc.method(t, "value", [](std::optional<T> const &t) {return t.value();});
    }
};

/******************************************************************************/

template <class ...Ts>
struct Renderer<std::variant<Ts...>> {
    void operator()(Document &doc) const {
        TypeIndex const t{typeid(std::variant<Ts...>)};
        doc.type(t, "std.Variant");
        doc.method(t, "get", [](std::variant<Ts...> const &v) {
            return std::visit([](auto const &v) -> rebind::Variable {return v;}, v);
        });
    }
};

/******************************************************************************/

template <class T>
struct ImplicitConversions<T, std::void_t<typename T::base_type>> {
    using Base = typename T::base_type;
    using types = decltype(concat(Pack<Base>(), typename ImplicitConversions<Base>::types()));
};

/******************************************************************************/

template <>
struct Renderer<bool> {
    void operator()(Document &doc) const {doc.type(typeid(bool), "std.Bool");}
};

template <class T>
struct Renderer<T, std::enable_if_t<std::is_floating_point_v<T>>> {
    void operator()(Document &doc) const {doc.type(typeid(T), "std.Float");}
};

template <class T>
struct Renderer<T, std::enable_if_t<std::is_integral_v<T>>> {
    void operator()(Document &doc) const {doc.type(typeid(T), "std.Integer");}
};

/******************************************************************************/

template <class T>
struct Renderer<T, std::enable_if_t<
    arma::is_arma_type<std::decay_t<T>>::value
    || arma::is_arma_cube_type<std::decay_t<T>>::value
    || arma::is_arma_sparse_type<std::decay_t<T>>::value
>> : NoRender {};

template <class A>
static constexpr bool is_arma_dense = std::is_same_v<A, std::decay_t<A>> && (arma::is_Mat<A>::value || arma::is_Col<A>::value || arma::is_Cube<A>::value);

template <class A>
struct Response<A, Value, std::enable_if_t<(is_arma_dense<A>)>> {
    bool operator()(Variable &out, TypeIndex t, A const &a) const {
        if (t.equals<ArrayView>())
            return out.emplace(Type<ArrayView>(), a.memptr(), ArrayLayout{nupack::la::shape(a), nupack::la::strides(&a)}), true;
        return false;
    }
};

/******************************************************************************/

template <class A>
struct Response<A, Value, std::enable_if_t<(arma::is_arma_sparse_type<A>::value)>> {
    using T = typename A::elem_type;

    bool operator()(Variable &out, TypeIndex t, A const &a) const {
        if (t.equals<rebind::Sequence>()) {
            auto &s = *out.emplace(Type<Sequence>());
            s.emplace_back(Type<ArrayView>(), a.values, a.n_nonzero + 1);
            s.emplace_back(Type<ArrayView>(), a.row_indices, a.n_nonzero + 1);
            s.emplace_back(Type<ArrayView>(), a.col_ptrs, a.n_cols + 2);
            s.emplace_back(nupack::la::shape(a));
            return true;
        }
        return false;
    }
};

/******************************************************************************/

template <class A>
struct Request<A, std::enable_if_t<(is_arma_dense<A>)>> {
    using T = typename A::elem_type;
    std::optional<A> operator()(Variable const &r, Dispatch &msg, bool copy=true) const {
        if (Debug) std::cout << "convert to arma" << std::endl;
        if (auto p = r.request<ArrayView>()) {
            if (nupack::la::depth<A> != p->layout.depth())
                return msg.error("incorrect dimensions", typeid(A), nupack::la::depth<A>, p->layout.depth());
            if (p->layout.n_elem() != 0)
                if (p->layout.shape(0) != 1 && p->layout.stride(0) != 1)
                    return msg.error("array is not column-major", typeid(A), 1, p->layout.stride(0));
            if (auto t = p->data.target<T>()) {
                static_assert(!std::is_const_v<std::remove_pointer_t<decltype(t)>>);
                return nupack::la::dense_from_data<T, nupack::la::depth<A>>(t, p->layout, copy);
            } else return msg.error("incorrect value type", typeid(T));
        } else return msg.error("not convertible to armadillo type", typeid(A));
    }
};

/******************************************************************************/

template <class A>
struct Request<A &, std::enable_if_t<(is_arma_dense<A>)>> {
    A * operator()(Variable const &r, Dispatch &msg) const {
        if (auto a = Request<A>()(r, msg, false))
            return std::addressof(msg.storage.emplace_back().emplace<A>(std::move(*a)));
        return nullptr;
    }
};

/******************************************************************************/

template <class A>
struct Request<A const &, std::enable_if_t<(is_arma_dense<A>)>> {
    A const * operator()(Variable const &r, Dispatch &msg) const {
        if (auto a = Request<A>()(r, msg, false))
            return std::addressof(msg.storage.emplace_back().emplace<A>(std::move(*a)));
        return nullptr;
    }
};

/******************************************************************************/

template <class A>
struct Request<A, std::enable_if_t<(arma::is_arma_sparse_type<A>::value)>> {
    using T = typename A::elem_type;
    std::optional<A> operator()(Variable const &r, Dispatch &msg, bool copy=true) const {
        if (Debug) std::cout << "convert to arma sparse" << std::endl;
        static_assert(nupack::la::depth<arma::Col<arma::uword>> == 1);
        static_assert(nupack::la::depth<arma::Col<T>> == 1);
        if (auto p = r.request<std::tuple<arma::Col<arma::uword>, arma::Col<arma::uword>,
                                          arma::Col<T>, std::size_t, std::size_t>>()) {
            return A(std::move(std::get<0>(*p)), std::move(std::get<1>(*p)), std::move(std::get<2>(*p)), std::get<3>(*p), std::get<4>(*p));
        } else return msg.error("not convertible to armadillo sparse type", typeid(A));
    }
};

/******************************************************************************/

template <class T, std::size_t N, class A>
struct Request<boost::container::small_vector<T, N, A>> : VectorRequest<boost::container::small_vector<T, N, A>> {};

/******************************************************************************/

template <class T, std::size_t N, class A>
struct Response<boost::container::small_vector<T, N, A>> : VectorResponse<boost::container::small_vector<T, N, A>> {};

/******************************************************************************/

template <class V>
struct VectorRenderer {
    void operator()(Document &doc) const {
        doc.render<typename V::value_type>();
        if constexpr(std::is_same_v<typename V::value_type, char>) {
            doc.type(typeid(V), "std.String");
            if constexpr(!std::is_same_v<typename V::iterator, typename V::const_iterator>) {
                doc.method(typeid(V), "append", [](V &v, typename V::value_type o) {v.push_back(std::move(o));});
            }
        } else {
            doc.type(typeid(V), "std.Vector");
            doc.method(typeid(V), "append", [](V &v, typename V::value_type o) {v.emplace_back(std::move(o));});
        }
        doc.method(typeid(V), "[]", [](V &v, std::size_t i) -> decltype(v.at(i)) {return v.at(i);});
        doc.method(typeid(V), "__len__", [](V const &v) {
            if (rebind::Debug) BEEP(v.size());
            return v.size();});
        doc.method(typeid(V), "value_type", [](V const &) {return type_index<typename V::value_type>();});
    }
};

/******************************************************************************/

template <class T, class A>
struct Renderer<std::vector<T, A>> : VectorRenderer<std::vector<T, A>> {};

template <class T, class C, class A>
struct Renderer<std::basic_string<T, C, A>> : VectorRenderer<std::basic_string<T, C, A>> {};

template <class T, class C>
struct Renderer<std::basic_string_view<T, C>> : VectorRenderer<std::basic_string_view<T, C>> {};

template <class T, std::size_t N, class A>
struct Renderer<boost::container::small_vector<T, N, A>> : VectorRenderer<boost::container::small_vector<T, N, A>> {};

/******************************************************************************/

template <>
struct Renderer<nupack::json> {
    void operator()(Document &doc) const;
};

}

/******************************************************************************/

namespace nupack {

using rebind::Document;
using rebind::Type;

/******************************************************************************/

template <class T>
void render_json(Document &doc, Type<T> t) {
    doc.method(t, "to_json", [](T const &x) {return json(x);});
    doc.method(t, "from_json", [](json const &s) {T t; s.get_to(t); return t;});
}

/******************************************************************************/

template <class T>
void render_public(Document &doc, Type<T> t) {
    if constexpr(has_names<T>)
        for_each_zip(T::names(), T::accesses(), [&](char const *n, auto f) {
            if constexpr(std::is_member_object_pointer_v<decltype(f)>) {
                doc.method(t, n-1, f);
            }
        });
    else if (!std::is_empty_v<T>) std::cout << "not rendering members of type " << typeid(T).name() << std::endl;
}

/******************************************************************************/

template <class T>
void render_comparisons(Document &doc, Type<T> t) {
    if constexpr(has_eq<T>) doc.method(t, "==", std::equal_to<T>());
    if constexpr(has_ne<T>) doc.method(t, "!=", std::not_equal_to<T>());

    if constexpr(has_lt<T>) doc.method(t, "<", std::less<T>());
    if constexpr(has_gt<T>) doc.method(t, ">", std::greater<T>());

    if constexpr(has_le<T>) doc.method(t, "<=", std::less_equal<T>());
    if constexpr(has_ge<T>) doc.method(t, ">=", std::greater_equal<T>());
}

/******************************************************************************/

template <class T>
void render_hash(Document &doc, Type<T> t) {
    doc.method(t, "__hash__", [](T const &t) {return std::hash<T>()(t);});
}

/******************************************************************************/

/// Render a view of characters
template <class V, NUPACK_IF(is_view<V> && is_character<value_type_of<V>>)>
void render(Document &doc, Type<V> t) {doc.render<value_type_of<V>>();}

template <class V, NUPACK_IF(is_view<V> && is_character<value_type_of<V>>)>
auto response(std::type_index, V const &v) {
    std::stringstream ss;
    dump_os(ss, v);
    return ss.str();
}

/****************************************************************************************/

/// Convert a view of objects to a Python list
template <class V, NUPACK_IF(is_view<V> && !is_character<value_type_of<V>>)>
void render(Document &doc, Type<V> t) {doc.render<value_type_of<V>>();}

template <class V, NUPACK_IF(is_view<V> && !is_character<value_type_of<V>>)>
std::vector<value_type_of<V>> response(std::type_index, V const &v) {
    return {begin_of(v), end_of(v)};
}

/**************************************************************************************/

template <class T>
struct Dumpable {
    std::string operator()(T const &t) const {
        std::ostringstream os;
        dump_os(os, t);
        return os.str();
    }
};

template <class T>
Dumpable<T> dumpable(Type<T> t={}) {return {};}

/**************************************************************************************/

template <class C>
rebind::Dictionary to_dictionary(C c) {
    rebind::Dictionary out;
    out.reserve(tuple_size<decltype(names_of(c))>);
    for_each_zip(names_of(c), members_of(c), [&](auto c, auto &m) {
        out.emplace_back(std::string_view(c), std::move(m));
    });
    std::sort(out.begin(), out.end(),
        [](auto const &i, auto const &j) {return i.first < j.first;});
    return out;
}

template <class C>
bool from_dictionary(rebind::Dictionary v, C &c, rebind::Dispatch &msg) {
    bool ok = true;
    for_each_zip(names_of(c), members_of(c), [&](std::string_view c, auto &m) {
        if (!ok) return;
        auto it = std::lower_bound(v.begin(), v.end(), c,
            [&](auto const &p, auto const &s) {return p.first < s;});
        if (it == v.end() || it->first != c) ok = false;
        if (auto v = it->second.template request<std::decay_t<decltype(m)>>()) m = std::move(*v);
    });
    return ok;
}

/**************************************************************************************/

template <class C>
std::optional<C> from_dictionary(rebind::Dictionary v, rebind::Dispatch &msg) {
    C c;
    if (from_dictionary(std::move(v), c, msg)) return c;
    return msg.error("member not found");
}

/**************************************************************************************/

}

#include "Core.h"
#include "Design.h"
#include "Math.h"
#include "Model.h"
#include "Thermo.h"
