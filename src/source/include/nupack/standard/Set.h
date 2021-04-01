#pragma once
#include "../reflect/Memory.h"
#include "../reflect/Print.h"
#include "../reflect/Hash.h"

#include <set>
#include <unordered_set>

namespace nupack {

/******************************************************************************************/

NUPACK_DEFINE_VARIADIC(is_set, std::set, class);
NUPACK_EXTEND_VARIADIC(is_set, std::unordered_set, class);
NUPACK_EXTEND_VARIADIC(is_set, std::multiset, class);

template <class T>
struct hash<T, void_if<(has_hash<value_type_of<T>>) && (is_set<T>)>> : RangeHash<T> {};

/******************************************************************************************/

template <class ...Ts> using Set = std::set<Ts...>;

template <class Key, class Hash=hash<Key>, class Equal=std::equal_to<Key>, class Alloc=std::allocator<Key>>
using UnorderedSet = std::unordered_set<Key, Hash, Equal, Alloc>;

/******************************************************************************************/

template <class T>
struct memory::impl<T, void_if<is_set<T>>> {
    std::size_t operator()(T const &t) const {return sum(t, memory::impl<value_type_of<T>>());}
    void erase(T &t) const {T t_; t.swap(t_);}
};

/******************************************************************************************/

template <class T> struct io::PrintAsContainer<T, void_if<is_set<T>>> : PrintAsSet {};

template <class T>
Set<value_type_of<T>> make_set(T &&t) {return {begin_of(t), end_of(t)};}

}
