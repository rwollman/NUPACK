#pragma once
#include <utility>
#include <typeindex>
#include <string_view>
#include <iostream>
#include <functional>
#include <typeindex>

namespace rebind {

template <class T>
using unqualified = std::remove_cv_t<std::remove_reference_t<T>>;

/******************************************************************************************/

enum Qualifier : unsigned char {Value, Const, Lvalue, Rvalue};

static std::string_view QualifierNames[4] = {"value", "const", "lvalue", "rvalue"};
static std::string_view QualifierSuffixes[4] = {"", " const &", " &", " &&"};

inline std::ostream & operator<<(std::ostream &os, Qualifier q) {
    return os << QualifierNames[static_cast<unsigned char>(q)];
}

template <class T, Qualifier Q> using qualified =
    std::conditional_t<Q == Value, T,
    std::conditional_t<Q == Const, T const &,
    std::conditional_t<Q == Lvalue, T &, T &&>>>;

template <class T>
static constexpr Qualifier qualifier_of = (!std::is_reference_v<T>) ? Value
    : (std::is_rvalue_reference_v<T> ? Rvalue : (std::is_const_v<std::remove_reference_t<T>> ? Const : Lvalue));

inline constexpr Qualifier add(Qualifier a, Qualifier b) {return a == Value ? b : a;}

/******************************************************************************************/

/// Compile-time type a la boost::hana
template <class T>
struct Type  {
    using type = T;
    T operator*() const; // undefined
    constexpr Type<unqualified<T>> operator+() const {return {};}
};

template <class T, class U>
constexpr std::is_same<T, U> is_same(Type<T>, Type<U>) {return {};}

template <class T>
static constexpr Type<T> ctype = {};

template <class T>
struct IndexedType {
    std::size_t index;
    T operator*() const; // undefined
    constexpr Type<unqualified<T>> operator+() const {return {};}
};

/******************************************************************************************/

using Demangler = std::function<std::string(char const *)>;

std::string demangle(char const *);

void set_demangler(Demangler fun) noexcept;

/******************************************************************************************/

class TypeIndex {
    std::pair<std::type_info const *, Qualifier> p;
    constexpr TypeIndex(std::type_info const *t, Qualifier q=Value) noexcept : p(t, q) {}

public:

    constexpr TypeIndex() noexcept : p(nullptr, Value) {};
    constexpr TypeIndex(std::type_info const &t, Qualifier q=Value) noexcept : p(&t, q) {}

    template <class T>
    TypeIndex(Type<T>) noexcept : p(&typeid(T), qualifier_of<T>) {}

    /**************************************************************************************/

    std::type_info const & info() const noexcept {return p.first ? *p.first : typeid(void);}
    std::string name() const noexcept {return demangle(info().name());}
    constexpr Qualifier qualifier() const noexcept {return p.second;}

    operator std::type_info const &() const noexcept {return info();}
    operator std::type_index() const noexcept {return info();}

    /// For now, hash code does not incorporate the qualifier
    std::size_t hash_code() const noexcept {return info().hash_code();}

    constexpr void set_qualifier(Qualifier q) noexcept {p.second = q;}

    /// Return if the index is not empty
    constexpr explicit operator bool() const noexcept {return p.first;}

    /// Test if this type equals another one, but ignoring all qualifiers
    template <class T>
    bool matches(Type<T> t={}) const noexcept {return p.first && typeid(T) == *p.first;}

    /// Test if this type equals another one, but ignoring all qualifiers
    bool matches(TypeIndex const &t) const noexcept {return p.first == t.p.first;}

    /// Test if this type equals a type specified as a compile time argument
    template <class T>
    bool equals(Type<T> t={}) const noexcept {return p.first && typeid(T) == *p.first && qualifier_of<T> == p.second;}

    /// Add a qualifier obeying the usual C++ semantics
    constexpr TypeIndex add(Qualifier q) const noexcept {return {p.first, p.second == Value ? q : p.second};}

    /**************************************************************************************/

    constexpr bool operator==(TypeIndex const &t) const {return p == t.p;}
    constexpr bool operator!=(TypeIndex const &t) const {return p != t.p;}
    constexpr bool operator<(TypeIndex const &t) const {return p < t.p;}
    constexpr bool operator>(TypeIndex const &t) const {return p > t.p;}
    constexpr bool operator<=(TypeIndex const &t) const {return p <= t.p;}
    constexpr bool operator>=(TypeIndex const &t) const {return p >= t.p;}

    /// Return a copy without any qualifiers
    constexpr TypeIndex operator+() const {return {p.first, Value};}
};


inline std::ostream & operator<<(std::ostream &os, TypeIndex t) {
    return os << t.name() << QualifierSuffixes[static_cast<unsigned char>(t.qualifier())];
}

template <class T>
constexpr TypeIndex type_index(Type<T> t={}) {return t;}

/******************************************************************************************/

}

namespace std {

template <>
struct hash<rebind::TypeIndex> {
    size_t operator()(rebind::TypeIndex const &t) const {return t.hash_code();}
};

}
