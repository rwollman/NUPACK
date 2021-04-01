#pragma once
#include <vector>
#include <boost/container/small_vector.hpp>
#include <boost/container/static_vector.hpp>

#include "../reflect/Memory.h"
#include "../reflect/Print.h"
#include "../reflect/Hash.h"
#include "../algorithms/Utility.h"

namespace nupack {

/******************************************************************************************/

template <class T, class Alloc=std::allocator<T>> using vec = std::vector<T, Alloc>;
template <class T, std::size_t s=16> using small_vec = boost::container::small_vector<T, s>;
template <class T, std::size_t s> using static_vec = boost::container::static_vector<T, s>;

//template <class T> using default_vec = if_t<sizeof(T) < 64, boost::container::small_vector<T, 16>, std::vector<T>>;
template <class T> using default_vec = boost::container::small_vector<T, 16>;

NUPACK_DEFINE_VARIADIC(is_vec, std::vector, class);
NUPACK_DEFINE_VARIADIC(is_small_vec, boost::container::small_vector, class, std::size_t, class);
NUPACK_DEFINE_VARIADIC(is_static_vec, boost::container::static_vector, class, std::size_t);

/// Copy a container to a std::vector
template <class V>
vec<value_type_of<V>> as_vec(V const &v) {return {begin_of(v), end_of(v)};}

/// Copy a container to a small vector
template <std::size_t N=16, class V>
small_vec<value_type_of<V>, N> as_small_vec(V const &v) {return {begin_of(v), end_of(v)};}

/******************************************************************************************/

template <class T>
struct hash<T, void_if<(has_hash<value_type_of<T>>) && (is_vec<T> || is_small_vec<T> || is_static_vec<T>)>> : RangeHash<T> {};

/******************************************************************************************/

/// Memory count for vectors sums across elements and currently ignores the sizeof()
template <class T> struct memory::impl<T, void_if<is_vec<T> || is_static_vec<T> || is_small_vec<T>>> {
    std::size_t operator()(T const &t) const {return sum(t, memory::impl<value_type_of<T>>());}
    void erase(T &t) const {T t_; t.swap(t_);}
};

template <class T> struct io::PrintAsContainer<T, void_if<is_vec<T> || is_static_vec<T> || is_small_vec<T>>> : PrintAsList {};

}

/******************************************************************************************/
