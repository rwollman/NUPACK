#pragma once
#include <deque>

#include "../reflect/Memory.h"
#include "../reflect/Print.h"
#include "../reflect/Hash.h"

namespace nupack {

/******************************************************************************************/

template <class ...Ts>
using deque = std::deque<Ts...>;

NUPACK_DEFINE_VARIADIC(is_deque, deque, class);

template <class T>
struct hash<T, void_if<(has_hash<value_type_of<T>>) && (is_deque<T>)>> : RangeHash<T> {};

/******************************************************************************************/

template <class T>
struct memory::impl<T, void_if<is_deque<T>>> {
    std::size_t operator()(T const &t) const {return sum(t, memory::impl<value_type_of<T>>());}
    void erase(T &t) const {T t_; t.swap(t_);}
};

template <class T> struct io::PrintAsContainer<T, void_if<is_deque<T>>> : public PrintAsList {};

/******************************************************************************************/

}
