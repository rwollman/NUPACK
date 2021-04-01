/**
 * @brief Iterator utilities including a iterator wrapper adaptor
 *
 * @file Iterator.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../algorithms/Traits.h"

#include <iterator>
#include <memory>
#include <array>

namespace nupack { namespace detail {

/******************************************************************************************/

/// Base class containing value_type if it is available
template <class T, class=void> struct value_type_base {};
template <class T> struct value_type_base<T, void_t<value_type_of<T>>> {using value_type = value_type_of<T>;};
/// Base class containing reference_type if it is available
template <class T, class=void> struct reference_type_base {};
template <class T> struct reference_type_base<T, void_t<reference_type_of<T>>> {using reference = reference_type_of<T>;};
/// Base class containing pointer if it is available
template <class T, class=void> struct pointer_type_base {};
template <class T> struct pointer_type_base<T, void_t<pointer_type_of<T>>> {using pointer = pointer_type_of<T>;};
/// Base class containing difference_type if it is available
template <class T, class=void> struct difference_type_base {};
template <class T> struct difference_type_base<T, void_t<difference_type_of<T>>> {using difference_type = difference_type_of<T>;};
/// Base class containing size_type if it is available
template <class T, class=void> struct size_type_base {};
template <class T> struct size_type_base<T, void_t<size_type_of<T>>> {using size_type = size_type_of<T>;};
/// Base class containing iterator_category if it is available
template <class T, class=void> struct iterator_category_base {};
template <class T> struct iterator_category_base<T, void_t<iterator_category_of<T>>> {using iterator_category = iterator_category_of<T>;};
template <class T> struct iterator_category_base<T, void_if<is_arithmetic<T> || is_pointer<T>>> {using iterator_category = std::random_access_iterator_tag;};

/******************************************************************************************/

}

/// The WrapIter object contains a given iterator, providing all deducible traits and methods
/// It is meant to be a base class for overriding those traits and methods selectively
template <class C, class T> struct WrapIter :
    detail::value_type_base<T>, detail::reference_type_base<T>,
    detail::pointer_type_base<T>, detail::iterator_category_base<T>,
    detail::size_type_base<T>, detail::difference_type_base<T> {

protected:
    T iter;
    constexpr decltype(auto) crtp() {return static_cast<C &>(*this);}
    constexpr decltype(auto) crtp() const {return static_cast<C const &>(*this);}

public:
    using simple_type = std::true_type;

    constexpr WrapIter(T t) : iter{std::move(t)} {}
    T original() const {return iter;}

    /**************************************************************************************/

    template <int B=1, NUPACK_IF(B && can_increment<T>)> constexpr decltype(auto) operator++() {++iter; return crtp();}
    template <int B=1, NUPACK_IF(B && can_decrement<T>)> constexpr decltype(auto) operator--() {--iter; return crtp();}

    template <class N> constexpr auto operator+=(N &&n) -> decltype(iter += n, declref<C>()) {iter += n; return crtp();}
    template <class N> constexpr auto operator-=(N &&n) -> decltype(iter -= n, declref<C>()) {iter -= n; return crtp();}

    template <int B=1, NUPACK_IF(B && can_subtract<T>)>
    constexpr auto operator-(C const &c2) const {return iter - c2.iter;}

    template <int B=1, NUPACK_IF(B && has_eq<T>)> constexpr auto operator==(C const &o) const {return iter == o.iter;}
    template <int B=1, NUPACK_IF(B && has_lt<T>)> constexpr auto operator< (C const &o) const {return iter < o.iter;}

    template <int B=1, NUPACK_IF(B && can_arrow<T>)> constexpr decltype(auto) operator->() const {return iter;}
    template <int B=1, NUPACK_IF(B && can_dereference<T>)> constexpr decltype(auto) operator*() const {return *iter;}

    template <class N> constexpr auto operator[](N &&n) const -> decltype(iter[n]) {return iter[n];}

    /**************************************************************************************/

    template <class C_=C, NUPACK_IF(can_copy<C_> && can_increment<C_>)>
    constexpr auto operator++(int) {C_ copy{crtp()}; ++crtp(); return copy;}

    template <class C_=C, NUPACK_IF(can_copy<C_> && can_decrement<C_>)>
    constexpr auto operator--(int) {C_ copy{crtp()}; --crtp(); return copy;}

    template <class N> friend constexpr auto operator+(C c, N &&n) -> no_ref<decltype(c += fw<N>(n))> {return c += fw<N>(n);}
    template <class N> friend constexpr auto operator+(N &&n, C c) -> no_ref<decltype(c += fw<N>(n))> {return c += fw<N>(n);}
    template <class N> friend constexpr auto operator-(C c, N &&n) -> no_ref<decltype(c -= fw<N>(n))> {return c -= fw<N>(n);}
    template <class N> friend constexpr auto operator-(N &&n, C c) -> no_ref<decltype(c -= fw<N>(n))> {return c -= fw<N>(n);}

    template <int B=1, NUPACK_IF(B && has_lt<T>)> constexpr auto operator>(C const &o) const {return o < crtp();}
    template <int B=1, NUPACK_IF(B && has_eq<T>)> constexpr auto operator!=(C const &o) const {return !(crtp() == o);}
    template <int B=1, NUPACK_IF(B && has_lt<T> && has_eq<T>)> constexpr auto operator<=(C const &o) const {return !(crtp() > o);}
    template <int B=1, NUPACK_IF(B && has_lt<T> && has_eq<T>)> constexpr auto operator>=(C const &o) const {return !(crtp() < o);}
};

/******************************************************************************************/

/// Iterator which prevents access to non-const dereferenced value
template <class T> struct ConstIter : WrapIter<ConstIter<T>, T> {
    constexpr auto const & operator*() const {return *(this->iter);}
    constexpr auto const * operator->() const {return arrow(this->iter);}
    //template <class N>
    //constexpr auto const & operator[](N &&n) const -> decltype((*this).iter) {return (*this + n).iter;}
};

template <class U, class T> struct CastIter : WrapIter<CastIter<U, T>, T> {
    using value_type = U;
    using reference = U &&;
    constexpr U operator*() const {return static_cast<U>(*(this->iter));}
    using base_type = WrapIter<CastIter<U, T>, T>;
    using base_type::base_type;
};

template <class U, class T>
auto cast_iter(T &&t) {return CastIter<U, no_qual<T>>(t);}

/******************************************************************************************/

/// Iterator which is the same as its dereferenced value
template <class T> struct ValueIter : WrapIter<ValueIter<T>, T> {
    using value_type = T;
    using pointer = T;
    using reference	= T const &;
    using base_type = WrapIter<ValueIter<T>, T>;
    using base_type::base_type;

    constexpr auto const & operator*() const {return this->iter;}
    constexpr auto const * operator->() const {return &(this->iter);}

    template <class N> constexpr auto operator[](N &&n) const -> decltype((*this + n).iter) {return (*this + n).iter;}
};

/******************************************************************************************/

/// Iterator which yields same object over and over whenever dereferenced
template <class T> struct CopiesIter : WrapIter<CopiesIter<T>, std::size_t> {
    T const *value;
    using reference = T const &;
    using value_type = T;
    using pointer = T const *;

    constexpr CopiesIter(T const *v, std::size_t n) : WrapIter<CopiesIter<T>, std::size_t>(n), value(v) {}

    constexpr reference operator*() const {return *value;}
    constexpr decltype(auto) operator->() const {return &(*(*this));}
    constexpr reference operator[](std::size_t) const {return *(*this);}
};

/******************************************************************************************/

/// IndirectIter wraps an iterator with a *const* mapping function
template <class T, class F> class IndirectIter : public WrapIter<IndirectIter<T, F>, T> {

    template <class P, NUPACK_IF(is_lref<P>)>
    static auto reference_wrap(P p) {return &p;};

    template <class P, NUPACK_IF(is_rref<P> || !is_ref<P>)>
    static auto reference_wrap(P p) {return std::make_unique<no_ref<P>>(std::move(p));};

public:
    using base_type = typename IndirectIter::WrapIter;
    using base_type::iter;
    F const map;

    using pointer = pointer_type_of<T>;
    using reference = to_reference<decltype(declref<F const>()(*declref<T const>()))>;
    using value_type = decay<reference>;

    constexpr IndirectIter(T t, F f) : base_type(std::move(t)), map{std::move(f)} {}
    constexpr IndirectIter & operator=(IndirectIter it) {iter = std::move(it.iter); return *this;}
    constexpr IndirectIter(IndirectIter const &) = default;
    constexpr IndirectIter(IndirectIter &&) = default;

    constexpr decltype(auto) operator*() const {return map(*iter);}
    constexpr decltype(auto) operator->() const {return reference_wrap<decltype(map(*iter))>(map(*iter));}
    constexpr decltype(auto) operator()() const {return *iter;}
    constexpr decltype(auto) operator+() const {return iter;}

    template <class N>
    constexpr auto operator[](N &&n) const -> decltype(map(iter[n])) {return map(iter[n]);}
};

/******************************************************************************************/

/// ReverseIter is like std::reverse_iterator with a little bit of added deduction
template <class T> struct ReverseIter : WrapIter<ReverseIter<T>, T> {
    using base_type = typename ReverseIter::WrapIter;
    using base_type::iter;

    using value_type = value_type_of<T>;
    using pointer = pointer_type_of<T>;
    using reference	= reference_type_of<T>;

    using base_type::base_type;

    template <int B=1, NUPACK_IF(B && can_decrement<T>)> constexpr decltype(auto) operator++() {--iter; return *this;}
    template <int B=1, NUPACK_IF(B && can_increment<T>)> constexpr decltype(auto) operator--() {++iter; return *this;}

    template <int B=1, NUPACK_IF(B && (can_copy<T>) && (can_increment<T>))>
    constexpr auto operator++(int) {ReverseIter copy{*this}; --iter; return copy;}

    template <int B=1, NUPACK_IF(B && (can_copy<T>) && (can_decrement<T>))>
    constexpr auto operator--(int) {ReverseIter copy{*this}; ++iter; return copy;}

    template <class N> constexpr auto operator+=(N &&n) -> decltype(iter -= n, declref<ReverseIter>()) {iter -= n; return *this;}
    template <class N> constexpr auto operator-=(N &&n) -> decltype(iter += n, declref<ReverseIter>()) {iter += n; return *this;}

    template <int B=1, NUPACK_IF(B && can_subtract<T>)>
    constexpr auto operator-(ReverseIter const &c2) const {return c2.iter - iter;}

    template <int B=1, NUPACK_IF(B && has_lt<T>)> constexpr auto operator<(ReverseIter const &o) const {return o.iter < iter;}

    template <class N> constexpr auto operator[](N &&n) const -> decltype(iter[-n]) {return iter[-n];}
};

template <class T>
auto reverse_iter(T &&t) {return ReverseIter<no_qual<T>>{t};}

/******************************************************************************************/

}
