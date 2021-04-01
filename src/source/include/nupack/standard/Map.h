#pragma once
#include "../reflect/Memory.h"
#include "../reflect/Print.h"
#include "../reflect/Hash.h"

#include <map>
#include <unordered_map>

namespace nupack {

/******************************************************************************************/

template <class ...Ts> using Map = std::map<Ts...>;

template <class Key, class T, class Hash=hash<Key>, class Equal=std::equal_to<Key>,
          class Alloc=std::allocator<std::pair<Key const, T>>>
using HashMap = std::unordered_map<Key, T, Hash, Equal, Alloc>;

NUPACK_DEFINE_VARIADIC(is_map, std::map, class);
NUPACK_EXTEND_VARIADIC(is_map, std::multimap, class);
NUPACK_EXTEND_VARIADIC(is_map, std::unordered_map, class);

template <class T>
struct hash<T, void_if<(has_hash<value_type_of<T>>) && (is_map<T>)>> : RangeHash<T> {};

/******************************************************************************************/

template <class T>
struct memory::impl<T, void_if<is_map<T>>> {
    std::size_t operator()(T const &t) const {return sum(t, memory::impl<value_type_of<T>>());}
    void erase(T &t) const {T t_; t.swap(t_);}
};

/******************************************************************************************/

template <class T> struct io::PrintAsContainer<T, void_if<is_map<T>>> : PrintAsSet {};

/******************************************************************************************/

template <class M=void, class V>
auto count_map(V const &v) {
    nonvoid<M, Map<value_type_of<V>, std::size_t>> out;
    for (auto const &t : v) out->try_emplace(t, 0).first->second += 1;
    return out;
}

/******************************************************************************************/

}
