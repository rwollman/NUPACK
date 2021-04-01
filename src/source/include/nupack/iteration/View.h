/**
 * @brief Defines View<T> class, basically a iterator pair referencing another container
 *
 * @file View.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Iterator.h"

#include "../algorithms/Numeric.h"
#include "../algorithms/Operators.h"

#include <cstring>
#include <vector>
#include <string>
#include <cwchar>

namespace nupack {

/// Non-owning view object that holds a pair of iterators and represents the range between them
template <class T> class View {
    T m_begin, m_end;

public:
    using const_iterator = T;
    using iterator = T;

    using value_type = value_type_of<T>;
    using size_type = size_type_of<T>;
    using difference_type = difference_type_of<T>;

    /******************************************************************************************/

    constexpr if_t<can_copy<T>, T, T const &> begin() const {return m_begin;}
    constexpr if_t<can_copy<T>, T, T const &> end() const {return m_end;}

    template <class B> auto set_begin(B &&b) -> decltype(m_begin = fw<B>(b)) {return m_begin = fw<B>(b);}
    template <class E> auto set_end(E &&e) -> decltype(m_end = fw<E>(e)) {return m_end = fw<E>(e);}

    /******************************************************************************************/

    static constexpr std::array<char const *, 2> names() {return {"begin", "end"};};
    auto members() {return std::tie(m_begin, m_end);};

    template <class B=zero_t, class E=zero_t, NUPACK_IF(can_construct<T, B &&> && can_construct<T, E &&>)>
    constexpr View(B &&b={}, E &&e={}) : m_begin(fw<B>(b)), m_end(fw<E>(e)) {};

    constexpr size_type size() const {return static_cast<size_type>(std::distance(begin(), end()));}
    constexpr bool empty() const {return end() == begin();}
    explicit constexpr operator bool() const {return !empty();}

    constexpr View offset(difference_type i) const {return {std::next(begin(), i), end()};}
    constexpr View offset(difference_type i, difference_type j) const {return {std::next(begin(), i), std::next(end(), j)};}

    constexpr friend View operator>>(View const &r, difference_type i) {return {r.begin(), std::next(r.end(), i)};}
    constexpr friend View operator<<(View const &r, difference_type i) {return {r.begin(), std::next(r.end(), -i)};}
    constexpr friend View operator>>(difference_type i, View const &r) {return {std::next(r.begin(), i), r.end()};}
    constexpr friend View operator<<(difference_type i, View const &r) {return {std::next(r.begin(), -i), r.end()};}

    template <int B=1, NUPACK_IF(B && can_decrement<T>)>
    auto rbegin() const {return reverse_iter(std::prev(end()));}
    template <int B=1, NUPACK_IF(B && can_decrement<T>)>
    auto rend() const {return reverse_iter(std::prev(begin()));}

    /******************************************************************************************/

    template <class V>
    void assign(V &&v) {std::copy(begin_of(v), end_of(v), begin());}

    template <class F>
    void map(F &&f) {
        size_type z = 0;
        for (auto it = begin(); it != end(); ++it, ++z) *it = f(z);
    }

    /******************************************************************************************/

    /// Explicit conversion operator to STL style containers taking iterator pairs
    template <class V, NUPACK_IF(can_construct<V, T, T>)>
    constexpr operator V() const {return V{begin(), end()};}

    /// Element access by index, if it is allowed by the iterator
    template <class I>
    auto operator[](I &&i) const -> decltype(begin()[fw<I>(i)]) {return begin()[fw<I>(i)];}
};

/******************************************************************************************/

NUPACK_DEFINE_TEMPLATE(is_view, View, class);

/******************************************************************************************/

template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_lt<value_type_of<T>, value_type_of<U>>)>
constexpr bool operator<(T const &t, U const &u) {return std::lexicographical_compare(begin_of(t), end_of(t), begin_of(u), end_of(u));}

template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_eq<value_type_of<T>, value_type_of<U>> && !can_eq<const_iterator_of<T>, const_iterator_of<U>>)>
constexpr bool operator==(T const &t, U const &u) {return std::equal(begin_of(t), end_of(t), begin_of(u), end_of(u));}

template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_eq<const_iterator_of<T>, const_iterator_of<U>>)>
constexpr bool operator==(T const &t, U const &u) {
    if (begin_of(t) == begin_of(u) && end_of(t) == end_of(u)) return true;
    return std::equal(begin_of(t), end_of(t), begin_of(u), end_of(u));
}


template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_lt<value_type_of<T>, value_type_of<U>>)>
constexpr bool operator>(T const &t, U const &u)  {return u < t;}

template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_lt<value_type_of<T>, value_type_of<U>> && can_eq<value_type_of<T>, value_type_of<U>>)>
constexpr bool operator<=(T const &t, U const &u) {return !(u < t);}

template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_lt<value_type_of<T>, value_type_of<U>> && can_eq<value_type_of<T>, value_type_of<U>>)>
constexpr bool operator>=(T const &t, U const &u) {return !(t < u);}

template <class T, class U, NUPACK_IF(is_view<T> || is_view<U>), NUPACK_IF(can_eq<value_type_of<T>, value_type_of<U>>)>
constexpr bool operator!=(T const &t, U const &u)  {return !(t == u);}

/******************************************************************************************/

struct view_t {
    /// Make a view from two iterators
    template <class T> View<T> operator()(T const &b, T const &e) const {return View<T>(b, e);}
    /// Make a view from a container, start index, and end index
    template <class V, class B, class E>
    constexpr auto operator()(V && v, B b, E e) const {return (*this)(begin_of(v) + b, begin_of(v) + e);}
    /// Make a view from a container
    template <class V> constexpr auto operator()(V const & v) const {return (*this)(begin_of(v), end_of(v));}
};

static constexpr auto view = view_t();

/******************************************************************************************/

/// Iterator halfway through a container
template <class V> auto midpoint(V &&v) {return std::next(begin_of(v), len(v) / 2);}
/// Make a view from a C pointer and number of elements
template <class T> auto ptr_view(T t, std::ptrdiff_t n) {return view(t, t + n);}
/// Make a reversed view of a container
template <class V> constexpr auto reversed(V &&v) {
    return view(reverse_iter(std::prev(end_of(v))), reverse_iter(std::prev(begin_of(v))));
}
/// Make a view which puts a mapping onto the original elements
template <class V, class F> auto indirect_view(V &&v, F &&f) {
    auto func = to_functor(fw<F>(f));
    using iter = IndirectIter<decay<decltype(begin_of(v))>, decltype(func)>;
    return view(iter(begin_of(v), func), iter(end_of(v), func));
}
/// Make a view of keys of a container with std::pair-ish elements
template <class V> auto key_view(V &&v) {return indirect_view(fw<V>(v), first_of);}
/// Make a view of items of a container with std::pair-ish elements
template <class V> auto item_view(V &&v) {return indirect_view(fw<V>(v), second_of);}
/// Make two views for the first and second halves of a container
template <class V> auto bisect(V const &v) {
    return make_array(view(begin_of(v), midpoint(v)), view(midpoint(v), end_of(v)));
}
/// Make a view for the last n elements of a container
template <class V> auto last_n(V const &v, size_type_of<V> i) {return view(end_of(v) - i, end_of(v));}

/******************************************************************************************/

template <class V> auto moved_view(V &v) {return view(move_begin(v), move_end(v));}

/******************************************************************************************/

template <class V>
struct indexer {
    V v;
    template <class I> auto operator()(I &&i) const -> decltype(v[fw<I>(i)]) {return v[fw<I>(i)];}
    template <class I> auto operator()(I &&i) -> decltype(v[fw<I>(i)]) {return v[fw<I>(i)];}
};

template <class I, class V>
auto indexed_view(I &&i, V &&v) {return indirect_view(fw<I>(i), indexer<remove_rref<V>>{v});}

/******************************************************************************************/

template <class T> auto null_pos(T t)  {while (*t) ++t; return t;}
inline auto null_pos(char *t)          {return t + std::strlen(t);}
inline auto null_pos(char const *t)    {return t + std::strlen(t);}
inline auto null_pos(wchar_t *t)       {return t + std::wcslen(t);}
inline auto null_pos(wchar_t const *t) {return t + std::wcslen(t);}

/******************************************************************************************/

template <class R=void, class V, NUPACK_IF(!is_view<V>)>
V const & convert_from_view(V const &v) {return v;}

template <class R=void, class V, NUPACK_IF(is_view<V>)>
auto convert_from_view(V const &v) {return nonvoid<R, std::vector<value_type_of<V>>>{v};}

/******************************************************************************************/

}
