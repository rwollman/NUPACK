#pragma once
#include "Variable.h"
#include "Conversions.h"
#include <cstdlib>

namespace rebind {

/******************************************************************************/

/// Built-in integer with the largest domain
#ifdef INTPTR_MAX
    using Integer = std::intptr_t;
#else
    using Integer = std::ptrdiff_t;
#endif

/// Built-in floating point with the largest domain -- long double is not used though.
using Real = double;

using Sequence = Vector<Variable>;

using Dictionary = Vector<std::pair<std::string_view, Variable>>;

/******************************************************************************/

template <class B, class E>
bool array_major(B begin, E end) {
    for (std::ptrdiff_t expected = 1; begin != end; ++begin) {
        if (begin->first < 2) continue;
        if (begin->second != expected) return false;
        expected = begin->first * begin->second;
    }
    return true;
}

struct ArrayLayout {
    /// Stores a list of shape, stride pairs
    using value_type = std::pair<std::size_t, std::ptrdiff_t>;
    Vector<value_type> contents;

    ArrayLayout() = default;

    template <class T, class U>
    ArrayLayout(T const &shape, U const &stride) {
        if (std::size(shape) != std::size(stride)) throw std::invalid_argument("ArrayLayout() strides do not match shape");
        contents.reserve(std::size(shape));
        auto st = std::begin(stride);
        for (auto const &sh : shape) contents.emplace_back(sh, *st++);
    }

    // 1 dimensional contiguous array
    ArrayLayout(std::size_t n) : contents{value_type(n, 1)} {}

    friend std::ostream & operator<<(std::ostream &os, ArrayLayout const &l) {
        os << "ArrayLayout(" << l.depth() << "):" << std::endl;
        for (auto const &p : l.contents) os << p.first << ": " << p.second << " "; os << std::endl;
        return os;
    }

    /// Return stride of given dimension
    std::ptrdiff_t stride(std::size_t i) const {return contents[i].second;}

    /// Return shape[i] for a given dimension i
    std::size_t shape(std::size_t i) const {return contents[i].first;}

    /// Synonym of shape(): Return shape[i] for a given dimension i
    std::size_t operator[](std::size_t i) const {return contents[i].first;}

    bool column_major() const {return array_major(contents.begin(), contents.end());}
    bool row_major() const {return array_major(contents.rbegin(), contents.rend());}

    /// Number of dimensions
    std::size_t depth() const {return contents.size();}

    /// Total number of elements
    std::size_t n_elem() const {
        std::size_t out = contents.empty() ? 0u : 1u;
        for (auto const &p : contents) out *= p.first;
        return out;
    }
};

/******************************************************************************/

/*
Binary data convenience wrapper for an array of POD data
 */
class ArrayData {
    void *ptr;
    std::type_info const *t;
    bool mut;

public:
    void const *pointer() const {return ptr;}
    bool mutate() const {return mut;}
    std::type_info const &type() const {return t ? *t : typeid(void);}

    ArrayData(void *p, std::type_info const *t, bool mut) : ptr(p), t(t), mut(mut) {}

    template <class T>
    ArrayData(T *t) : ArrayData(const_cast<std::remove_cv_t<T> *>(static_cast<T const *>(t)),
                                &typeid(std::remove_cv_t<T>), std::is_const_v<T>) {}

    template <class T>
    T * target() const {
        if (!mut && !std::is_const<T>::value) return nullptr;
        if (type() != typeid(std::remove_cv_t<T>)) return nullptr;
        return static_cast<T *>(ptr);
    }

    friend std::ostream & operator<<(std::ostream &os, ArrayData const &d) {
        if (!d.t) return os << "ArrayData(<empty>)";
        return os << "ArrayData(" << TypeIndex(*d.t, d.mut ? Const : Lvalue) << ")";
    }
};

/******************************************************************************/

struct ArrayView {
    ArrayData data;
    ArrayLayout layout;
};

/******************************************************************************/

template <class T>
struct Request<T *> {
    std::optional<T *> operator()(Variable const &v, Dispatch &msg) const {
        std::optional<T *> out;
        if (!v || v.request<std::nullptr_t>(msg)) out.emplace(nullptr);
        else if (auto p = v.request<T &>(msg)) out.emplace(std::addressof(*p));
        return out;
    }
};

/******************************************************************************/

template <>
struct Response<char const *> {
    bool operator()(Variable &out, TypeIndex const &t, char const *s) const {
        if (t.equals<std::string_view>()) {
            out = s ? std::string_view(s) : std::string_view();
            return true;
        } else if (t.equals<std::string>()) {
            if (s) out.emplace(Type<std::string>(), s);
            else out.emplace(Type<std::string>());
            return true;
        } else return false;
    }
};

template <>
struct Request<char const *> {
    std::optional<char const *> operator()(Variable const &v, Dispatch &msg) const {
        std::optional<char const *> out;
        if (!v || v.request<std::nullptr_t>()) out.emplace(nullptr);
        else if (auto p = v.request<std::string_view>(msg)) out.emplace(p->data());
        else if (auto p = v.request<char const &>(msg)) out.emplace(std::addressof(*p));
        return out;
    }
};

/******************************************************************************/

using Binary = std::basic_string<unsigned char>;

using BinaryView = std::basic_string_view<unsigned char>;

class BinaryData {
    unsigned char *m_begin=nullptr;
    unsigned char *m_end=nullptr;
public:
    constexpr BinaryData() = default;
    constexpr BinaryData(unsigned char *b, std::size_t n) : m_begin(b), m_end(b + n) {}
    constexpr auto begin() const {return m_begin;}
    constexpr auto data() const {return m_begin;}
    constexpr auto end() const {return m_end;}
    constexpr std::size_t size() const {return m_end - m_begin;}
    operator BinaryView() const {return {m_begin, size()};}
};

template <>
struct Response<BinaryData> {
    bool operator()(Variable &out, TypeIndex const &t, BinaryData const &v) const {
        if (t.equals<BinaryView>()) return out.emplace(Type<BinaryView>(), v.begin(), v.size()), true;
        return false;
    }
};

template <>
struct Response<BinaryView> {
    bool operator()(Variable &out, TypeIndex const &, BinaryView const &v) const {
        return false;
    }
};

template <>
struct Request<BinaryView> {
    std::optional<BinaryView> operator()(Variable const &v, Dispatch &msg) const {
        if (auto p = v.request<BinaryData>()) return BinaryView(p->data(), p->size());
        return msg.error("not convertible to binary view", typeid(BinaryView));
    }
};

template <>
struct Request<BinaryData> {
    std::optional<BinaryData> operator()(Variable const &v, Dispatch &msg) const {
        return msg.error("not convertible to binary data", typeid(BinaryData));
    }
};

/******************************************************************************/

template <class T>
struct Response<T, Value, std::enable_if_t<(std::is_integral_v<T>)>> {
    bool operator()(Variable &out, TypeIndex const &i, T t) const {
        DUMP("response from integer", typeid(T).name(), i.name());
        if (i == typeid(Integer)) return out = static_cast<Integer>(t), true;
        if (i == typeid(Real)) return out = static_cast<Real>(t), true;
        DUMP("no response from integer");
        return false;
    }
};


/*
Default Response for floating point allows conversion to Real or Integer
*/
template <class T>
struct Response<T, Value, std::enable_if_t<(std::is_floating_point_v<T>)>> {
    bool operator()(Variable &out, TypeIndex const &i, T t) const {
        if (i == typeid(Real)) return out = static_cast<Real>(t), true;
        if (i == typeid(Integer)) return out = static_cast<Integer>(t), true;
        return false;
    }
};

/*
Default Request for integer type tries to go through double precision
long double is not expected to be a useful route (it's assumed there are not multiple floating types larger than Real)
*/
template <class T>
struct Request<T, std::enable_if_t<std::is_floating_point_v<T>>> {
    std::optional<T> operator()(Variable const &v, Dispatch &msg) const {
        DUMP("convert to floating");
        if (!std::is_same_v<Real, T>) if (auto p = v.request<Real>()) return static_cast<T>(*p);
        return msg.error("not convertible to floating point", typeid(T));
    }
};


/*
Default Request for integer type tries to go through Integer
*/
template <class T>
struct Request<T, std::enable_if_t<std::is_integral_v<T>>> {
    std::optional<T> operator()(Variable const &v, Dispatch &msg) const {
        DUMP("trying convert to arithmetic", v.type(), typeid(T).name());
        if (!std::is_same_v<Integer, T>) if (auto p = v.request<Integer>()) return static_cast<T>(*p);
        DUMP("failed to convert to arithmetic", v.type(), typeid(T).name());
        return msg.error("not convertible to integer", typeid(T));
    }
};

/*
Default Request for enum permits conversion from integer types
*/
template <class T>
struct Request<T, std::enable_if_t<std::is_enum_v<T>>> {
    std::optional<T> operator()(Variable const &v, Dispatch &msg) const {
        DUMP("trying convert to enum", v.type(), typeid(T).name());
        if (auto p = v.request<std::underlying_type_t<T>>()) return static_cast<T>(*p);
        return msg.error("not convertible to enum", typeid(T));
    }
};

/*
Default Response for enum permits conversion to integer types
*/
template <class T>
struct Response<T, Value, std::enable_if_t<(std::is_enum_v<T>)>> {
    bool operator()(Variable &out, TypeIndex const &i, T t) const {
        if (i == typeid(std::underlying_type_t<T>))
            return out = static_cast<std::underlying_type_t<T>>(t), true;
        if (i == typeid(Integer))
            return out = static_cast<Integer>(t), true;
        return false;
    }
};

/*
Default Request for string tries to convert from std::string_view and std::string
*/
template <class T, class Traits, class Alloc>
struct Request<std::basic_string<T, Traits, Alloc>> {
    std::optional<std::basic_string<T, Traits, Alloc>> operator()(Variable const &v, Dispatch &msg) const {
        DUMP("trying to convert to string");
        if (auto p = v.request<std::basic_string_view<T, Traits>>())
            return std::basic_string<T, Traits, Alloc>(std::move(*p));
        if (!std::is_same_v<std::basic_string<T, Traits, Alloc>, std::basic_string<T, Traits>>)
            if (auto p = v.request<std::basic_string<T, Traits>>())
                return std::move(*p);
        return msg.error("not convertible to string", typeid(T));
    }
};

template <class T, class Traits>
struct Request<std::basic_string_view<T, Traits>> {
    std::optional<std::basic_string_view<T, Traits>> operator()(Variable const &v, Dispatch &msg) const {
        return msg.error("not convertible to string view", typeid(T));
    }
};


/******************************************************************************/

/*
    Response for CompileSequence concept -- a sequence of compile time length
*/
template <class V>
struct CompiledSequenceResponse {
    using Array = std::array<Variable, std::tuple_size_v<V>>;

    template <std::size_t ...Is>
    static Sequence sequence(V const &v, std::index_sequence<Is...>) {
        Sequence o;
        o.reserve(sizeof...(Is));
        (o.emplace_back(std::get<Is>(v)), ...);
        return o;
    }

    template <std::size_t ...Is>
    static Sequence sequence(V &&v, std::index_sequence<Is...>) {
        Sequence o;
        o.reserve(sizeof...(Is));
        (o.emplace_back(std::get<Is>(std::move(v))), ...);
        return o;
    }

    template <std::size_t ...Is>
    static Array array(V const &v, std::index_sequence<Is...>) {return {std::get<Is>(v)...};}

    template <std::size_t ...Is>
    static Array array(V &&v, std::index_sequence<Is...>) {return {std::get<Is>(std::move(v))...};}

    bool operator()(Variable &out, TypeIndex const &t, V const &v) const {
        auto idx = std::make_index_sequence<std::tuple_size_v<V>>();
        if (t == typeid(Sequence)) return out = sequence(v, idx), true;
        if (t == typeid(Array)) return out = array(v, idx), true;
        return false;
    }

    bool operator()(Variable &out, TypeIndex const &t, V &&v) const {
        auto idx = std::make_index_sequence<std::tuple_size_v<V>>();
        if (t == typeid(Sequence)) return out = sequence(std::move(v), idx), true;
        if (t == typeid(Array)) return out = array(std::move(v), idx), true;
        return false;
    }
};

template <class T>
struct Response<T, Value, std::enable_if_t<std::tuple_size<T>::value >= 0>> : CompiledSequenceResponse<T> {};

/******************************************************************************/

template <class V>
struct CompiledSequenceRequest {
    using Array = std::array<Variable, std::tuple_size_v<V>>;

    template <class ...Ts>
    static void combine(std::optional<V> &out, Ts &&...ts) {
        DUMP("trying CompiledSequenceRequest combine", bool(ts)...);
        if ((bool(ts) && ...)) out = V{*static_cast<Ts &&>(ts)...};
    }

    template <class S, std::size_t ...Is>
    static void request_each(std::optional<V> &out, S &&s, Dispatch &msg, std::index_sequence<Is...>) {
        DUMP("trying CompiledSequenceRequest request");
        combine(out, std::move(s[Is]).request(msg, Type<std::tuple_element_t<Is, V>>())...);
    }

    template <class S>
    static void request(std::optional<V> &out, S &&s, Dispatch &msg) {
        DUMP("trying CompiledSequenceRequest request");
        if (std::size(s) != std::tuple_size_v<V>) {
            msg.error("wrong sequence length", typeid(V), std::tuple_size_v<V>, s.size());
        } else {
            msg.indices.emplace_back(0);
            request_each(out, std::move(s), msg, std::make_index_sequence<std::tuple_size_v<V>>());
            msg.indices.pop_back();
        }
    }

    std::optional<V> operator()(Variable r, Dispatch &msg) const {
        std::optional<V> out;
        DUMP("trying CompiledSequenceRequest", r.type().name());
        if constexpr(!std::is_same_v<V, Array>) {
            if (auto p = r.request<std::array<Variable, std::tuple_size_v<V>>>()) {
                DUMP("trying array CompiledSequenceRequest2", r.type().name());
                request(out, std::move(*p), msg);
            }
            return out;
        }
        if (auto p = r.request<Sequence>()) {
            DUMP("trying CompiledSequenceRequest2", r.type().name());
            request(out, std::move(*p), msg);
        } else {
            DUMP("trying CompiledSequenceRequest3", r.type().name());
            msg.error("expected sequence to make compiled sequence", typeid(V));
        }
        return out;
    }
};

/// The default implementation is to accept convertible arguments or Variable of the exact typeid match
template <class T>
struct Request<T, std::enable_if_t<std::tuple_size<T>::value >= 0>> : CompiledSequenceRequest<T> {};

/******************************************************************************/

template <class R, class V>
R from_iters(V &&v) {return R(std::make_move_iterator(std::begin(v)), std::make_move_iterator(std::end(v)));}

template <class R, class V>
R from_iters(V const &v) {return R(std::begin(v), std::end(v));}

/******************************************************************************/

template <class T, class=void>
struct HasData : std::false_type {};

template <class T>
struct HasData<T, std::enable_if_t<(std::is_pointer_v<decltype(std::data(std::declval<T>()))>)>> : std::true_type {};

/******************************************************************************/

template <class T, class Iter1, class Iter2>
bool range_response(Variable &o, TypeIndex const &t, Iter1 b, Iter2 e) {
    if (t.equals<Sequence>()) {
        auto p = o.emplace(Type<Sequence>());
        p->reserve(std::distance(b, e));
        for (; b != e; ++b) p->emplace_back(Type<T>(), *b);
        return true;
    }
    if (t.equals<Vector<T>>()) return o.emplace(Type<Vector<T>>(), b, e), true;
    return false;
}

template <class V>
struct VectorResponse {
    using T = std::decay_t<typename V::value_type>;

    bool operator()(Variable &o, TypeIndex const &t, V const &v) const {
        if (range_response<T>(o, t, std::begin(v), std::end(v))) return true;
        if constexpr(HasData<V const &>::value)
            return o.emplace(Type<ArrayView>(), std::data(v), std::size(v)), true;
        return false;
    }

    bool operator()(Variable &o, TypeIndex const &t, V &v) const {
        if (range_response<T>(o, t, std::cbegin(v), std::cend(v))) return true;
        if constexpr(HasData<V &>::value)
            return o.emplace(Type<ArrayView>(), std::data(v), std::size(v)), true;
        return false;
    }

    bool operator()(Variable &o, TypeIndex const &t, V &&v) const {
        return range_response<T>(o, t, std::make_move_iterator(std::begin(v)), std::make_move_iterator(std::end(v)));
    }
};

template <class T, class A>
struct Response<std::vector<T, A>, Value, std::enable_if_t<!std::is_same_v<T, Variable>>> : VectorResponse<std::vector<T, A>> {};

/******************************************************************************/

template <class V>
struct VectorRequest {
    using T = std::decay_t<typename V::value_type>;

    template <class P>
    static std::optional<V> get(P &pack, Dispatch &msg) {
        V out;
        out.reserve(pack.size());
        msg.indices.emplace_back(0);
        for (auto &x : pack) {
            if (auto p = std::move(x).request(msg, Type<T>()))
                out.emplace_back(std::move(*p));
            else return msg.error();
            ++msg.indices.back();
        }
        msg.indices.pop_back();
        return out;
    }

    std::optional<V> operator()(Variable const &v, Dispatch &msg) const {
        if (auto p = v.request<ArrayView>()) {
            if (auto t = p->data.target<T const>()) {
                return V(t, t + p->layout.n_elem());
            }
        }
        // if (auto p = v.request<Vector<T>>()) return get(*p, msg);
        if (!std::is_same_v<V, Sequence>)
            if (auto p = v.request<Sequence>()) return get(*p, msg);
        return msg.error("expected sequence", typeid(V));
    }
};

template <class T, class A>
struct Request<std::vector<T, A>> : VectorRequest<std::vector<T, A>> {};

}
