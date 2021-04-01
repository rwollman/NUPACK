/**
 * @file Cast.h
 * @brief
 */
#pragma once
#include "API.h"

namespace rebind {

/******************************************************************************/

Object python_cast(Variable &&v, Object const &t, Object const &root);

/******************************************************************************/

// Unambiguous conversions from some basic C++ types to Objects

inline Object as_object(Object o) {return std::move(o);}

inline Object as_object(bool b) {return {b ? Py_True : Py_False, true};}

inline Object as_object(Integer i) {return {PyLong_FromLongLong(static_cast<long long>(i)), false};}

inline Object as_object(Real x) {return {PyFloat_FromDouble(x), false};}

inline Object as_object(std::string const &s) {return {PyUnicode_FromStringAndSize(s.data(), s.size()), false};}

inline Object as_object(std::string_view s) {return {PyUnicode_FromStringAndSize(s.data(), s.size()), false};}

inline Object as_object(BinaryView s) {return {PyByteArray_FromStringAndSize(reinterpret_cast<char const *>(s.data()), s.size()), false};}

inline Object as_object(Binary const &s) {return {PyByteArray_FromStringAndSize(reinterpret_cast<char const *>(s.data()), s.size()), false};}

/******************************************************************************/

// Initialize an object that has a direct Python wrapped equivalent
template <class T>
Object default_object(T t) {
    auto o = Object::from(PyObject_CallObject(type_object<T>(), nullptr));
    cast_object<T>(o) = std::move(t);
    return o;
}

inline Object as_object(TypeIndex t) {return default_object(std::move(t));}
inline Object as_object(Function t) {return default_object(std::move(t));}

/******************************************************************************/

/// Source driven conversion: guess the correct Python type from the source type
/// I guess this is where automatic class conversions should be done?
inline Object as_deduced_object(Variable &&ref) {
    DUMP("asking for object");
    if (!ref) return {Py_None, true};
    if (auto v = ref.request<Object>())           return std::move(*v);
    if (auto v = ref.request<Real>())             return as_object(std::move(*v));
    if (auto v = ref.request<Integer>())          return as_object(std::move(*v));
    if (auto v = ref.request<bool>())             return as_object(std::move(*v));
    if (auto v = ref.request<std::string_view>()) return as_object(std::move(*v));
    if (auto v = ref.request<std::string>())      return as_object(std::move(*v));
    if (auto v = ref.request<Function>())         return as_object(std::move(*v));
    if (auto v = ref.request<TypeIndex>())  return as_object(std::move(*v));
    if (auto v = ref.request<Binary>())           return as_object(std::move(*v));
    if (auto v = ref.request<BinaryView>())       return as_object(std::move(*v));
    if (auto v = ref.request<Sequence>())
        return map_as_tuple(std::move(*v), [](auto &&x) {return as_deduced_object(std::move(x));});
    return {};
}

/******************************************************************************/

}
