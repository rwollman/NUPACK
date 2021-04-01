/**
 * @brief Range<T> template for Python-like range(), also for iterators() view
 *
 * @file Range.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "View.h"

#include "../algorithms/Traits.h"
#include "../algorithms/Extents.h"
#include "../algorithms/Functor.h"
#include "../algorithms/Numeric.h"

namespace arma {
    class span;
}

namespace nupack {

/******************************************************************************************/

/// CRTP providing common functionality for Range and ReverseRange
template <class R, class I>
class RangeBase : TotallyOrdered {
    using diff = typename ValueIter<I>::difference_type;
public:
    constexpr I start() const {return *static_cast<R const &>(*this).begin();}
    constexpr I stop() const {return *static_cast<R const &>(*this).end();}

    constexpr bool operator==(RangeBase const &o) const {
        return start() == o.start() && stop() == o.stop();
    }
    constexpr bool operator<(RangeBase const &o) const {
        return start() < o.start() || (start() == o.start() && stop() < o.stop());
    }

    std::array<I, 2> save_repr() const {return {start(), stop()};}
    void load_repr(std::array<I, 2> const &a) {static_cast<R &>(*this) = R(a[0], a[1]);}

    friend std::ostream & operator<<(std::ostream &os, R const &s) {
        os << '[' << s.start() << ':' << s.stop();
        if (s.stride() != 1) os << ':' << s.stride();
        return os << ')';
    }

    constexpr R shift(diff i, diff j) const {
        return {static_cast<R const &>(*this).begin()[i], static_cast<R const &>(*this).end()[j]};
    }
    constexpr R shift(diff i) const {
        return {static_cast<R const &>(*this).begin()[i], static_cast<R const &>(*this).end()[i]};
    }
};

/******************************************************************************************/

template <class I> struct ReverseRange;

/**
 * @brief Container such that dereferencing the iterator just returns the iterator
 * Used so that range-based for loop can be used for indices and iterators
 * @tparam I type of the underlying iterator or integer
 */
template <class I> struct Range : View<ValueIter<I>>, RangeBase<Range<I>, I> {
    static_assert(can_copy<ValueIter<I>>, "must be copyable");
    using base_type = typename Range::View;
    using base_type::begin; using base_type::end;

    template <class B, class E>
    constexpr Range(B b, E e) : base_type(static_cast<I>(b), static_cast<I>(e)) {}

    template <class E>
    explicit constexpr Range(E e) : base_type(static_cast<I>(0u), static_cast<I>(e)) {}

    constexpr Range() : base_type(static_cast<I>(0u), static_cast<I>(0u)) {}

    /// Shift a span
    friend constexpr Range operator+(I i, Range const &s) {return {*s.begin()+i, *s.end()+i};}
    friend constexpr Range operator+(Range const &s, I i) {return {*s.begin()+i, *s.end()+i};}
    friend constexpr Range operator-(Range const &s, I i) {return {*s.begin()-i, *s.end()-i};}

    /// Shift and reverse a Range
    constexpr Range reversed(I n) const {return {n-*end(), n-*begin()};}
    friend constexpr Range operator-(I i, Range const &s) {return s.reversed(i);}

    constexpr ReverseRange<I> operator~() const {return {*std::prev(end()), *std::prev(begin())};}

    operator arma::span() const; // include Matrix.h for this definition

    constexpr int stride() const {return +1;}
};

/******************************************************************************************/

using span = Range<uint>;
static_assert(can_copy<span>, "");
using cspan = span const;

NUPACK_DEFINE_TYPE(is_span, span);

/******************************************************************************************/

/// Range which goes backward, created by calling ~ on Range
template <class I> struct ReverseRange : View<ReverseIter<ValueIter<I>>>, RangeBase<Range<I>, I> {
    static_assert(can_copy<ReverseIter<ValueIter<I>>>, "must be copyable");
    using base_type = typename ReverseRange::View;
    using base_type::begin; using base_type::end;

    constexpr ReverseRange() : base_type{I{}, I{}} {}
    constexpr ReverseRange(I i, I j) : base_type{i, j} {}

    constexpr Range<I> operator~() const {return {*std::prev(end()), *std::prev(begin())};}
    constexpr int stride() const {return -1;}
};

NUPACK_EXTEND_TEMPLATE(is_view, Range, class);

/******************************************************************************************/

/// Checks endpoint by !=
template <class T=void, class I> constexpr auto range(I i) {
    using S = nonvoid<T, I>;
    return Range<S>(static_cast<S>(zero), static_cast<S>(i));
}

/// Checks endpoint by <
template <class T=void, class I> constexpr auto lrange(I i) {
    using S = nonvoid<T, I>;
    return Range<S>(static_cast<S>(zero), std::max(static_cast<S>(zero), static_cast<S>(i)));
}

/// Checks endpoint by !=
template <class T=void, class I1, class I2> constexpr auto range(I1 i1, I2 i2) {
    using I = nonvoid<T, std::common_type_t<I1, I2>>;
    return Range<I>(static_cast<I>(i1), static_cast<I>(i2));
}

/// Checks endpoint by <
template <class T=void, class I1, class I2> constexpr auto lrange(I1 i1, I2 i2) {
    using I = nonvoid<T, std::common_type_t<I1, I2>>;
    return Range<I>(static_cast<I>(i1), std::max(static_cast<I>(i1), static_cast<I>(i2)));
}

/******************************************************************************************/

/// Return a view into a container between indices given by a span
template <class V, class T>
auto subview(V &&v, Range<T> const &s) {return view(fw<V>(v), s.start(), s.stop());}

/// Return a view of a container whose value_type is the iterator itself
template <class V> constexpr auto iterators(V &&v) {return range(begin_of(v), end_of(v));}

/// Return a view of a container whose value_type is the iterator itself, begin offset
template <class V, class B> constexpr auto iterators(V &&v, B b) {
    return range(std::next(begin_of(v), b), end_of(v));

}
/// Return a view of a container whose value_type is the iterator itself, begin and end offsets
template <class V, class B, class E> constexpr auto iterators(V &&v, B b, E e) {
    return range(std::next(begin_of(v), b), std::next(end_of(v), e));
}

/// Return a view on the indices of a container i.e. [0 : its length)
template <class I=void, class V> constexpr auto indices(V const &v) {
    using S = nonvoid<I, decay<decltype(len(v))>>;
    return range(static_cast<S>(zero), static_cast<S>(len(v)));
}

/// Return a view of copies of the same object
template <class T> constexpr auto copies(T const &t, std::size_t n) {
    return view(CopiesIter<T>{std::addressof(t), 0u}, CopiesIter<T>{std::addressof(t), n});
}

/******************************************************************************************/

}
