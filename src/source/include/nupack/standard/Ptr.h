#pragma once
#include <memory>

#include "../reflect/Memory.h"
#include "../reflect/Print.h"

namespace nupack {

/******************************************************************************************/

NUPACK_DEFINE_VARIADIC(is_shared_ptr, std::shared_ptr, class);
NUPACK_DEFINE_VARIADIC(is_unique_ptr, std::unique_ptr, class);

/******************************************************************************************/

/// Check for equality of two pointers' pointees, short-circuiting if pointers are the same
template <class T, NUPACK_IF(is_pointer<T> || is_shared_ptr<T> || is_unique_ptr<T>)>
bool equal_ptr(T const &t1, T const &t2) {return (t1 == t2) || (t1 && t2 && *t1 == *t2);}

/// Ordered comparison of two pointers' pointees, short-circuiting if pointers are the same
template <class T, NUPACK_IF(is_pointer<T> || is_shared_ptr<T> || is_unique_ptr<T>)>
bool less_ptr(T const &t1, T const &t2) {
    if (t1 == t2) return false;
    return (t1 && t2) ? *t1 < *t2 : bool(t1) < bool(t2);
}

/******************************************************************************************/

/// Divide memory usage of shared_ptr by number of use_count
template <class T> struct memory::impl<T, void_if<is_shared_ptr<T>>> {
    std::size_t operator()(T const &t) const {
        if (!t) return sizeof(T);
        else return sizeof(T) + memory::impl<element_type_of<T>>()(*t) / t.use_count();
    }
    void erase(T &t) const {t.reset();}
};

template <class T> struct memory::impl<T, void_if<is_unique_ptr<T>>> {
    std::size_t operator()(T const &t) const {
        if (!t) return sizeof(T);
        else return sizeof(T) + memory::impl<element_type_of<T>>()(*t);
    }
    void erase(T &t) const {t.reset();}
};

/******************************************************************************************/

/// Object that looks lke a pointer but actually holds the value itself on the stack
template <class T> class StackPtr {
    T value;
public:
    NUPACK_REFLECT(StackPtr, value);

    T const & operator*() const {return value;}
    T & operator*() {return value;}

    template <class ...Ts> StackPtr(Ts &&...ts) : value{fw<Ts>(ts)...} {}
};

/// Object that looks lke a pointer but actually holds the value itself on the heap
template <class T> class HeapPtr {
    std::unique_ptr<T> ptr;
public:
    NUPACK_REFLECT(HeapPtr, ptr);

    T const & operator*() const {return *ptr;}
    T & operator*() {return *ptr;}

    template <class ...Ts> HeapPtr(Ts &&...ts) : ptr(std::make_unique<T>(fw<Ts>(ts)...)) {}
};

template <class T, NUPACK_IF(is_lref<T>)> auto lref_capture(T &&t) {return &t;}
template <class T, NUPACK_IF(!is_lref<T>)> auto lref_capture(T &&t) {return StackPtr<decay<T>>(fw<T>(t));}

/******************************************************************************************/

}
