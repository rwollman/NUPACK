/**
 * @brief Tensor adapters and functions to avoid partition function overflow
 *
 * @file Overflow.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "SIMD.h"
#include "Tensor.h"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/adapted/std_pair.hpp>


namespace nupack {

NUPACK_DETECT(is_scalar_range, void_if<is_scalar<value_type_of<T>>>);
NUPACK_DETECT(is_compound_range, void_if<!is_scalar<value_type_of<T>>>);

/******************************************************************************************/

template <class T, class U>
using zip_iterator = boost::iterators::zip_iterator<std::pair<T, U>>;

/******************************************************************************************/

template <class T>
auto first_iter(T it) -> decltype(it.get_iterator_tuple().first) {return it.get_iterator_tuple().first;}

/******************************************************************************************/

template <class T>
auto second_iter(T it) -> decltype(it.get_iterator_tuple().second) {return it.get_iterator_tuple().second;}

/******************************************************************************************/

template <class T> using exponent_t = int_of_size<sizeof(T)>;
template <class T> using overflow = std::pair<T, exponent_t<T>>;

NUPACK_DEFINE_TEMPLATE(is_overflow, overflow, class);
NUPACK_TEMPLATE_FALLBACK(mantissa_t, typename T::first_type, T);

template <class T>
struct numeric_limits<overflow<T>> {
    static constexpr auto max_exponent = std::numeric_limits<exponent_t<T>>::max();
};

/******************************************************************************************/

namespace thermo {

/******************************************************************************************/

/// Tensor storage for overflow. The mantissa and exponent are stored separately
template <class T>
struct TensorBase<overflow<T>> {
    using man_type = vec<T, simd::allocator<T>>;
    using exp_type = vec<exponent_t<T>, simd::allocator<exponent_t<T>>>;
    std::pair<man_type, exp_type> storage;
    using iterator = zip_iterator<typename man_type::iterator, typename exp_type::iterator>;
    using const_iterator = zip_iterator<typename man_type::const_iterator, typename exp_type::const_iterator>;

    NUPACK_REFLECT(TensorBase, storage);

    using value_type = overflow<T>;
    using data_type = decltype(storage);

    TensorBase() = default;
    constexpr TensorBase(iseq i, T t={}) :
        storage(std::piecewise_construct, std::make_tuple(i, t), std::make_tuple(i, 0)) {}

    TensorBase(data_type &&d) : storage(std::move(d)) {}

    template <class B, class E, class O>
    static void read_span(B b, E e, O o) {
        if constexpr(std::is_scalar_v<std::decay_t<decltype(*b)>>) {
            std::transform(b, e, o, simd::ifrexp);
        } else {
            std::copy(first_iter(b), first_iter(e), first_iter(o));
            std::copy(second_iter(b), second_iter(e), second_iter(o));
        }
    }

    // move mantissa matrix
    explicit TensorBase(TensorBase<T> &&u) : storage(std::move(u.storage), std::size(u.storage)) {}

    // copy mantissa matrix, initialize exponents to 0
    template <class U, NUPACK_IF(std::is_scalar_v<U>)>
    explicit TensorBase(TensorBase<U> const &u) : storage(std::size(u.storage), std::size(u.storage)) {
        std::transform(u.storage.begin(), u.storage.end(), begin(), simd::ifrexp);
    }

    // copy mantissa and exponent matrices
    template <class U, NUPACK_IF(!std::is_scalar_v<U>)>
    explicit TensorBase(TensorBase<U> const &u) : storage(view(u.storage.first), view(u.storage.second)) {}

    template <class ...Ts> void resize(Ts ...ts) {storage.first.resize(ts...); storage.second.resize(ts...);}
    auto size() const {return storage.first.size();}
    constexpr auto strides() const {return std::array<iseq, 1>{1};}

    void fill(T t) {::nupack::fill(storage.first, t); ::nupack::fill(storage.second, 0);}

    void fill(T t, span s) { // fill but just for range
        std::fill(begin_of(storage.first) + s.start(), begin_of(storage.first) + s.stop(), t);
        if (t != 0) std::fill(begin_of(storage.second) + s.start(), begin_of(storage.second) + s.stop(), 0);
    }

    template <class I, class C>
    static void copy_out(I begin, I end, C &out) {
        std::copy(first_iter(begin), first_iter(end), std::back_inserter(out.first));
        std::copy(second_iter(begin), second_iter(end), std::back_inserter(out.second));
    }

    auto begin() {return iterator(std::make_pair(begin_of(storage.first), begin_of(storage.second)));}
    auto end() {return iterator(std::make_pair(end_of(storage.first), end_of(storage.second)));}
    auto begin() const {return const_iterator(std::make_pair(begin_of(storage.first), begin_of(storage.second)));}
    auto end() const {return const_iterator(std::make_pair(end_of(storage.first), end_of(storage.second)));}

    auto data() const {return begin();}
};

/******************************************************************************************/

template <class V, class I, NUPACK_IF(is_scalar_range<V> && is_integral<I>)>
auto mantissa_at(V &&v, I i) -> decltype(v[i]) {return v[i];}

template <class V, class I, NUPACK_IF(is_scalar_range<V> && is_integral<I>)>
auto mantissa_at(V &&v, I i) -> decltype(v(i)) {return v(i);}

template <class V, int N, NUPACK_IF(is_scalar_range<V>)>
decltype(auto) mantissa_at(V &&v, simd::Chunk<N> i) {
    return simd::load_pack<N>(std::addressof(*(begin_of(v) + i.value)));
}

template <class V, class I, NUPACK_IF(is_scalar_range<V>)>
Zero exponent_at(V &&v, I i) {return {};}

/******************************************************************************************/

template <class V, class I, NUPACK_IF(is_compound_range<V> && is_integral<I>)>
decltype(auto) mantissa_at(V &&v, I i) {return first_iter(begin_of(v))[i];}

template <class V, class I, NUPACK_IF(is_compound_range<V> && is_integral<I>)>
decltype(auto) exponent_at(V &&v, I i) {return second_iter(begin_of(v))[i];}

template <class V, int N, NUPACK_IF(is_compound_range<V>)>
decltype(auto) mantissa_at(V &&v, simd::Chunk<N> i) {
    return simd::load_pack<N>(&first_iter(begin_of(v))[i.value]);
}

template <class V, int N, NUPACK_IF(is_compound_range<V>)>
decltype(auto) exponent_at(V &&v, simd::Chunk<N> i) {
    return simd::load_pack<N>(&second_iter(begin_of(v))[i.value]);
}

/******************************************************************************************/

template <class P, class T=Zero>
auto mantissa(P &&p, T={}) -> decltype(fw<P>(p).first) {return fw<P>(p).first;}

template <class P, class T=Zero>
auto mantissa(P &&p, T={}) -> decltype(p->first) {return p->first;}

template <class P>
auto exponent(P &&p) -> decltype(p.second) {return p.second;}

template <class P>
auto exponent(P &&p) -> decltype(p->second) {return p->second;}

template <class P, class T=Zero>
auto exponent(P &&p, T hint) -> decltype(p.second + hint) {return p.second + hint;}

template <class P, class T=Zero>
auto exponent(P &&p, T hint) -> decltype(p->second + hint) {return p->second + hint;}

/******************************************************************************************/

template <class P, class T=Zero, NUPACK_IF(!is_class<decay<decltype(*std::declval<P>())>>)>
decltype(auto) mantissa(P &&p, T={}) {return *p;}

template <class P, class T=Zero, NUPACK_IF(!is_class<decay<P>>)>
remove_rref<P &&> mantissa(P &&p, T={}) {return fw<P>(p);}

template <class P, class T=Zero, NUPACK_IF(!is_class<decay<decltype(*std::declval<P>())>>)>
T exponent(P const &, T hint={}) {return hint;}

template <class P, class T=Zero, NUPACK_IF(!is_class<P>)>
T exponent(P const &, T hint={}) {return hint;}

/******************************************************************************************/

template <class V, NUPACK_IF(is_compound_range<V>)>
auto mantissa_view(V &&v) {return view(first_iter(begin_of(v)), first_iter(end_of(v)));}

template <class V, NUPACK_IF(is_compound_range<V>)>
auto exponent_view(V &&v) {return view(second_iter(begin_of(v)), second_iter(end_of(v)));}

template <class V, NUPACK_IF(is_scalar_range<V>)>
auto mantissa_view(V &&v) {return v;}

/******************************************************************************************/

template <class T>
static constexpr int overflow_bits = int(CHAR_BIT * sizeof(T)) / (is_overflow<T> ? -2 : 1);

template <bool B, class T> using oflow = if_t<B, overflow<T>, T>;

/******************************************************************************************/


}}

/******************************************************************************************/

namespace nupack::simd {

/******************************************************************************************/

/// Perform a map operation with SIMD, where the output is modified in place
template <class O, class F>
void map(O &&out, int i, int stop, F &&f) noexcept {

    if constexpr(std::is_scalar_v<value_type_of<O>>) {
        constexpr auto Z = pack_size<value_type_of<O>>;
        for (; i + Z <= stop; i += Z)
            bs::store(f(Chunk<Z>(i)).first, std::addressof(*(begin_of(out) + i)));
        for (; i < stop; ++i) {
            *(begin_of(out) + i) = f(i).first;
            static_assert(std::is_same_v<std::decay_t<decltype(f(i).second)>, Zero>);
        }
    } else {
        auto m = begin_of(thermo::mantissa_view(out));
        auto e = begin_of(thermo::exponent_view(out));
        constexpr auto Z = pack_size<value_type_of<decltype(m)>>;
        for (; i + Z <= stop; i += Z) {
            auto p = f(Chunk<Z>(i));
            bs::store(std::move(p.first), std::addressof(*(m + i)));
            bs::store(std::move(p.second), std::addressof(*(e + i)));
        }
        for (; i < stop; ++i) {
            std::tie(*(m + i), *(e + i)) = f(i);
        }
    }
}

template <class O, class F>
void map(O &&out, span s, F &&f) noexcept {
    return map(fw<O>(out), s.start(), s.stop(), fw<F>(f));
}

/******************************************************************************************/

}
