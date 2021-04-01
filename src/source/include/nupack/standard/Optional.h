#pragma once
#include <memory>
#include <optional>

#include "../reflect/Memory.h"
#include "../common/Error.h"
#include "../reflect/Print.h"

namespace nupack {

/******************************************************************************************/

NUPACK_DEFINE_VARIADIC(is_optional, std::optional, class);

template <class T> using Optional = std::optional<T>;

/// make an optional
template <class T>
Optional<decay<T>> optional(T &&t) {return fw<T>(t);}

/// Return empty optional type which may be used in e.g. "return bool() : optional(1) : optional()"
inline auto optional() {return std::nullopt;}

/// Safe accessor for optional
template <class T, NUPACK_IF(is_optional<decay<T>> && !DebugBounds)>
decltype(auto) value_of(T &&t) {return *(fw<T>(t));}

/// Safe accessor for optional
template <class T, NUPACK_IF(is_optional<decay<T>> && DebugBounds)>
decltype(auto) value_of(T &&t) {
   try {return fw<T>(t).value();}
   catch (...) {NUPACK_BUG("Empty optional was accessed", type_name(t));}
}

/******************************************************************************************/

/// Return container with only the non-empty values of another container
template <class R, class V> auto vmap_if(V &&v) {
    R ret; reserve_space(ret, len(v));
    for (auto i = begin_of(v); i != end_of(v); ++i)
        if (*i) ret.emplace_back(static_cast<if_t<is_lref<V>, decltype(value_of(*i)), no_ref<decltype(value_of(*i))> &&>>(value_of(*i)));
    return ret;
}

/// Return container with only the non-empty values of another container
template <template <class...> class R=std::vector, class V> auto vmap_if(V &&v) {
    using T = decay<decltype(**begin_of(v))>;
    return vmap_if<std::vector<T>>(fw<V>(v));
}

/******************************************************************************************/

/// Find element in container, return optional() of that element provided it exists
template <class V, class ...Ts>
auto search(V &&v, Ts &&...ts) {
    auto it = find(v, fw<Ts>(ts)...);
    Optional<no_qual<decltype(*it)>> out;
    if (it != end_of(v)) out = *it;
    return out;
}

/******************************************************************************************/

template <class T> struct memory::impl<T, void_if<is_optional<T>>> {
    std::size_t operator()(T const &t) const {
        if (!t) return sizeof(T); // assume that the optional has stack space for the held value
        else return sizeof(T) - sizeof(value_type_of<T>) + memory::impl<value_type_of<T>>()(*t);
    }
    void erase(T &t) const {t.reset();}
};

template <class T> struct io::Printer<T, void_if<is_optional<T> && !is_streamable<T>>> {
    void operator()(std::ostream &os, T const &t, io::Indent id={}) const {
        if (!t) os << "null";
        else Printer<value_type_of<T>>()(os, *t, id);
    }
};

/******************************************************************************************/

}
