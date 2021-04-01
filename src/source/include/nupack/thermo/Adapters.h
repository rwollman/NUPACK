/**
 * @brief Adapters to make Tensors show up as column major, row major, symmetric
 *
 * @file Adapters.h
 * @author Mark Fornace
 * @date 2018-06-01
 */
#pragma once
#include "Tensor.h"
#include "../standard/Optional.h"
#include "../types/Complex.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

/// Base class for 2D tensor that must be square
template <class T> struct SquareBase : Tensor<T, 2> {
    using base_type = Tensor<T, 2>;
    using base_type::base_type;

    SquareBase(Tensor<T, 2> t) : base_type(std::move(t)) {}

    // template <class ...Ts>
    // SquareBase(Ts &&...ts) : base_type(static_cast<Ts &&>(ts)...) {
    //     NUPACK_REQUIRE(base_type::shape()[0], ==, base_type::shape()[1]);
    // }

    auto subsquare(span s) {return base_type::operator()(s, s);}
    auto subsquare(span s) const {return base_type::operator()(s, s);}

    void copy_square(span i, span j) {
        copy_tensor_block(*this, i, i, *this, j, j);
    }

    /// preserve value (i, j) but set the exponent to the maximum of (i+1, j), (i, j-1), (i, j)
    void reset_exponent(uint i, uint j) {
        if constexpr(!std::is_scalar_v<std::decay_t<T>>) {
            auto &&y = *base_type::operator()(i, j);
            auto e = y.second + max(0, simd::ifrexp(y.first).second);
            if (i < j) e = max(e, max(exponent(base_type::operator()(i+1, j)), exponent(base_type::operator()(i, j-1))));
            if (i > j) e = max(e, max(exponent(base_type::operator()(i-1, j)), exponent(base_type::operator()(i, j+1))));
            mantissa(y) = simd::ldexp(mantissa(y), exponent(y) - e);
            exponent(y) = e;
        }
    }

    template <class A, class F>
    bool set_value(bool ij, uint i, uint j, A const &, F &&rule) {
        bool err = false;
        if constexpr(std::is_scalar_v<std::decay_t<T>>) {
            *base_type::operator()(i, j) = A::rig_type::element_value(err, static_cast<F &&>(rule), Zero());
        } else {
            typename std::decay_t<T>::second_type e0;
            if (i == j) e0 = 0;
            else if (ij) e0 = max(exponent(base_type::operator()(i+1, j)), exponent(base_type::operator()(i, j-1)));
            else e0 = max(exponent(base_type::operator()(i-1, j)), exponent(base_type::operator()(i, j+1)));
            *base_type::operator()(i, j) = A::rig_type::element_value(err, static_cast<F &&>(rule), e0);
        }
        return err;
    }

    base_type const & unglued() const {return *this;}
};

/******************************************************************************************/

/// Can be spanned on second index
template <class T> struct Upper : SquareBase<T> {
    using base_type = SquareBase<T>;
    using base_type::base_type;

    Upper(Tensor<T, 2> t) : base_type(std::move(t)) {}

    template <class U>
    Upper(Upper<U> const &u) : base_type(u) {
        if constexpr(!std::is_scalar_v<std::decay_t<T>>) {
            for (auto o : lrange(1, base_type::shape()[0]))
                for (auto i : range(base_type::shape()[0] - o))
                    base_type::reset_exponent(i, i + o);
        }
    }

    template <class A, class F>
    bool set(iseq i, iseq j, A const &a, F &&rule) {
        return base_type::set_value(true, i, j, a, static_cast<F &&>(rule));
    }
    auto write(span i, span j, bool) const {return base_type::write(i, j);}

    template <class M>
    void read(span is, span js, M const &m) {
        base_type::read(is, js, m);
        if constexpr(!std::is_scalar_v<std::decay_t<T>>) {
            for (auto o : range(max(is.stop(), js.start()) - is.stop(), js.stop() - is.start()))
                for (auto i : range(is.start(), js.stop() - o))
                    base_type::reset_exponent(i, i + o);
        }
    }
};

/******************************************************************************************/

constexpr auto reversed_strides(bool b) {return !b;} // C major to F major
template <class A> constexpr A reversed_strides(A a) {return {a[1], a[0]};}

/******************************************************************************************/

/// Can be spanned on first index
template <class T> struct Lower : SquareBase<T> {
    using base_type = SquareBase<T>;
    using base_type::base_type;

    auto strides() const {return reversed_strides(SquareBase<T>::strides());}

    Lower(Tensor<T, 2> t) : base_type(std::move(t)) {}

    template <class U>
    Lower(Lower<U> const &u) : base_type(u) {
        if constexpr(!std::is_scalar_v<std::decay_t<T>>)
            for (auto o : lrange(1, base_type::shape()[0]))
                for (auto i : range(base_type::shape()[0] - o))
                    base_type::reset_exponent(i + o, i);
    }

    decltype(auto) operator() (iseq i, iseq j) const {return base_type::operator()(j, i);}
    decltype(auto) operator() (iseq i, iseq j) {return base_type::operator()(j, i);}
    decltype(auto) operator() (span i, iseq j) const {return base_type::operator()(j, i);}

    template <class A, class F>
    bool set(iseq i, iseq j, A const &a, F &&rule) {
        return base_type::set_value(false, j, i, a, static_cast<F &&>(rule));
    }

    auto write(span i, span j, bool) const {return base_type::write(j, i);}

    template <class M>
    void read(span is, span js, M const &m) {
        base_type::read(js, is, m);
        if constexpr(!std::is_scalar_v<std::decay_t<T>>) {
            for (auto o : range(max(is.stop(), js.start()) - is.stop(), js.stop() - is.start()))
                for (auto i : range(is.start(), js.stop() - o))
                    base_type::reset_exponent(i + o, i);
        }
    }
};

/******************************************************************************************/

/// Can be spanned on first or second index, but must be kept symmetric
template <class T> struct Symmetric : SquareBase<T> {
    using base_type = SquareBase<T>;
    using base_type::base_type;

    Symmetric(Tensor<T, 2> t) : base_type(std::move(t)) {}

    template <class U>
    Symmetric(Symmetric<U> const &u) : base_type(u) {
        if constexpr(!std::is_scalar_v<std::decay_t<T>>)
            for (auto o : lrange(1, base_type::shape()[0]))
                for (auto i : range(base_type::shape()[0] - o)) {
                    base_type::reset_exponent(i, i + o);
                    *base_type::operator()(i + o, i) = *base_type::operator()(i, i + o);
                }
    }

    decltype(auto) operator() (iseq i, iseq j) const {return base_type::operator()(i, j);}
    decltype(auto) operator() (iseq i, iseq j) {return base_type::operator()(i, j);}
    decltype(auto) operator() (span i, iseq j) const {return base_type::operator()(j, i);}
    decltype(auto) operator() (iseq i, span j) const {return base_type::operator()(i, j);}

    template <class A, class F>
    bool set(iseq i, iseq j, A const &a, F &&rule) {
        bool out = base_type::set_value(true, i, j, a, static_cast<F &&>(rule));
        *base_type::operator()(j, i) = *base_type::operator()(i, j);
        return out;
    }

    auto write(span i, span j, bool) const {
        if (Debug) {
            for (auto x : i) for (auto y : j)
                NUPACK_DREQUIRE(*base_type::operator()(y, x), ==, *base_type::operator()(x, y));
        }
        return base_type::write(i, j);
    }

    template <class M>
    void read(span is, span js, M const &m) {
        base_type::read(is, js, m);
        if constexpr(!std::is_scalar_v<std::decay_t<T>>) {
            for (auto o : range(max(is.stop(), js.start()) - is.stop(), js.stop() - is.start()))
                for (auto i : range(is.start(), js.stop() - o))
                    base_type::reset_exponent(i, i + o);
        }
        for (auto i : is) for (auto j : js) *base_type::operator()(j, i) = *base_type::operator()(i, j);
    }
};

/******************************************************************************************/

NUPACK_EXTEND_VARIADIC(is_like_tensor, Lower, class);
NUPACK_EXTEND_VARIADIC(is_like_tensor, Upper, class);
NUPACK_EXTEND_VARIADIC(is_like_tensor, Symmetric, class);

/******************************************************************************************/

template <class T>
struct XTensor : MemberComparable {
    using value_type = no_qual<T>;
    using tensor_type = Tensor<value_type, 2>;
    using x_type = std::array<tensor_type, 3>;
    using V = small_vec<x_type>;
    using slice_type = if_t<is_cref<T>, View<const_iterator_of<V>>,
                       if_t<is_lref<T>, View<iterator_of<V>>, V>>;
    slice_type slices;
    small_vec<std::size_t> prefixes;
    NUPACK_REFLECT(XTensor, slices, prefixes);

    /// Only write anything if the X recursion is not complete
    auto write(span i, span, bool complete) const {
        return complete ? std::nullopt : std::make_optional(slices.at(sequence_index(i)));
    }

    template <class M>
    void read(span i, span, M const &X) {
        if (!X) return;
        auto &s = at(slices, sequence_index(i));
        s[0] = (*X)[0];
        s[1] = (*X)[1];
        s[2] = (*X)[2];
        if constexpr(!std::is_scalar_v<std::decay_t<T>>) {
            for (auto &&x : s) for (auto &&y : x) {
                auto p = simd::ifrexp(y.first);
                if (p.second > 0) {y.first = p.first; y.second += p.second;}
            }
        }
    }

    XTensor() = default;

    explicit XTensor(slice_type xs, small_vec<std::size_t> p) : slices(std::move(xs)), prefixes(std::move(p)) {}

    template <bool B=true, class U, NUPACK_IF(B && !is_ref<T>)>
    XTensor(Complex const &s, U u) {
        prefixes.emplace_back(0);
        slices.resize(s.n_strands());
        for (auto const i : range(s.n_strands()))
            prefixes.emplace_back(prefixes.back() + s.length(i));
    }

    /// Initialize memory for a single subblock
    template <bool B=true, class V, class T2,  NUPACK_IF(B && is_ref<T>)>
    void initialize(V const &seq, T2 const &zero) {
        NUPACK_REQUIRE(len(slices), ==, 1);
        auto const n = seq.length(0);
        // Maximum length of interior loop (in # of unpaired bases)
        auto const size = unsigned_minus(seq.n_strands() == 1 ? n : n + seq.length(seq.n_strands() - 1), 4);
        for (auto &s : slices[0]) {s.resize(n, size); s.fill(zero);}
    }

    template <class U, NUPACK_IF(!is_ref<T> && !is_ref<U>)>
    XTensor(XTensor<U> const &x) : slices(vmap<slice_type>(x.slices, [](auto const &s) {return x_type{s[0], s[1], s[2]};})), prefixes(x.prefixes) {
        if constexpr(std::is_scalar_v<std::decay_t<U>> && !std::is_scalar_v<std::decay_t<T>>) {
            for (auto &s : slices) for (auto &&x : s) for (auto &&y : x) {
                auto p = simd::ifrexp(y.first);
                NUPACK_DASSERT(std::isfinite(p.first) && std::isfinite(y.first), p, y);
                if (p.second > 0) {y.first = p.first; y.second += p.second;}
            }
        }
    }

    // get a span of sequences from a span of bases
    auto sequence_index(span s) const {
        auto i = find_index(prefixes, s.start());
        NUPACK_REQUIRE(i, <, len(slices), "span indices do not line up in X");
        NUPACK_ASSERT(contains(prefixes, prefixes[i] + len(s)), prefixes, s);
        return i;
    }

    auto subsquare(span s) {
        auto const i = sequence_index(s);
        return XTensor<T &>(view(slices, i, i+1), view(prefixes, i, i+2));
    }

    auto subsquare(span s) const {
        auto const i = sequence_index(s);
        return XTensor<T const &>(view(slices, i, i+1), view(prefixes, i, i+2));
    }

    void copy_square(span i, span j) const {} // don't think any copies are needed since X doesn't factor into higher calculations

    auto & operator[](iseq i) {return at(at(slices, 0), i);}
    auto const & operator[](iseq i) const {return at(at(slices, 0), i);}

    template <class A, class F>
    bool set(iseq i, iseq, A, F &&rule) { // j is implied already
        cspan S(0, (*this)[0].shape()[1]);
        return rule((*this)[0](i, S), add_const((*this)[2])(i+1, S));
    }

    void increment() {
        swap((*this)[2], (*this)[1]);
        swap((*this)[1], (*this)[0]);
    }
};

/******************************************************************************************/

/// Defines blocking behavior for a collection of vectors
/// To access a given vector i at element j, write vecs(i, j)
template <class T, std::size_t N>
struct Rows : Tensor<T, 2> {
    using base_type = Tensor<T, 2>;
    constexpr span all() const {return {0, N};}

    template <class ...Ts, NUPACK_IF(sizeof...(Ts) != 2)>
    Rows(Ts &&...ts) : base_type(fw<Ts>(ts)...) {NUPACK_REQUIRE(len(*this), ==, N);}

    template <bool B=true, NUPACK_IF(B && (is_same<T, decay<T>>))>
    Rows(iseq n, T t) : Rows(N, n, t) {}

    auto write(span i, span, bool) const {return base_type::write(all(), i);}
    template <class M> void read(span i, span, M const &m) {
        base_type::read(all(), i, m);
    }

    auto subsquare(span i) {return base_type::operator()(all(), i);}
    auto subsquare(span i) const {return base_type::operator()(all(), i);}

    void copy_square(span i, span j) {
        copy_tensor_block(*this, all(), i, *this, all(), j);
    }
};

/******************************************************************************************/

}}
