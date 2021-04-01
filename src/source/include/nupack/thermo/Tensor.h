#pragma once
#include "SIMD.h"
#include "../iteration/Range.h"
#include "../reflect/Reflection.h"
#include "../common/Error.h"
#include "../standard/Vec.h"
#include <ios>
#include <iomanip>

namespace nupack {

namespace thermo {

template <class T, NUPACK_IF(can_call<T>)>
decltype(auto) value_of(T &&t) noexcept {return t();}

template <class T, NUPACK_IF(!can_call<T> && !can_dereference<T>)>
decltype(auto) value_of(T &&t) noexcept {return fw<T>(t);}

template <class T, NUPACK_IF(can_dereference<T>)>
decltype(auto) value_of(T t) noexcept {return *t;}

/******************************************************************************************/

template <class T, std::size_t N> class Tensor;
NUPACK_DEFINE_TEMPLATE(is_tensor, Tensor, class, std::size_t);
NUPACK_DEFINE_TEMPLATE(is_like_tensor, Tensor, class, std::size_t);

/******************************************************************************************/

// template <class T, class U, class=void>
// struct ConvertTensor {
//     // template <class O, class M>
//     // static O make(M &&m) {return cast_container<O>(fw<M>(m));}

//     template <class B, class E, class O>
//     static void copy(B b, E e, O out) noexcept {std::copy(b, e, out);}
// };

// /// copy values from a container of type I to an output of value_type O
// template <class T, class U, class B, class E, class Out>
// void copy_tensor(B b, E e, Out out) {ConvertTensor<T, U>::copy(b, e, out);}

// /// make an O from
// template <class O, class T, class U, class I>
// O make_tensor(I &&i) {return ConvertTensor<T, U>::template make<O>(fw<I>(i));}

/******************************************************************************************/

/// Base class for a non subview tensor holds the data in a vec
template <class T>
struct TensorBase {
    using data_type = vec<T, simd::allocator<T>>;
    using iterator = decltype(declref<data_type>().begin());
    using const_iterator = decltype(declref<data_type const>().begin());
    using element_type = T;
    using value_type = T;

    data_type storage;
    NUPACK_REFLECT(TensorBase, storage);

    TensorBase() = default;
    explicit TensorBase(iseq n, T t={}) : storage(n, t) {}
    explicit TensorBase(data_type &&data) noexcept : storage(std::move(data)) {}

    template <class U>
    explicit TensorBase(TensorBase<U> const &u) : storage(view(u.storage)) {}

    template <class V>
    void assign(V const &v) noexcept {storage.assign(begin_of(v), end_of(v));}

    template <class B, class E, class O>
    static void read_span(B b, E e, O o) {std::copy(b, e, o);}

    auto data() const noexcept {return storage.begin();}

    auto begin() const noexcept {return storage.begin();}
    auto begin() noexcept {return storage.begin();}

    auto end() noexcept {return storage.end();}
    auto end() const noexcept {return storage.end();}

    auto rbegin() const noexcept {return storage.rbegin();}
    auto rend() const noexcept {return storage.rend();}

    template <class I, class C>
    static void copy_out(I const &b, I const &e, C &out) {std::copy(b, e, std::back_inserter(out));}

    void resize(iseq n) {storage.resize(n);}

    void fill(T t) noexcept {contiguous_fill(storage.data(), storage.data() + storage.size(), t);}

    void fill(T t, span s) noexcept {contiguous_fill(storage.data()+s.start(), storage.data()+s.stop(), t);}
};

/******************************************************************************************/

/// a struct to allow x(i) access instead of x[i] for iterators
template <class T> struct Offset {
    T iter;
    template <class I> decltype(auto) operator()(I &&i) const {return iter[fw<I>(i)];}
};

/// make a contiguous view starting at some iterator
template <class It>
auto offset_view(It it, span s) noexcept {return view(it + s.start(), it + s.stop());}

/******************************************************************************************/

/// 1 dimensional tensor
/// Tensor<overflow<T>, 1>
template <class T>
class Tensor<T, 1> : public TensorBase<T>, MemberComparable {
public:
    using base_type = TensorBase<T>;
    Tensor() = default;
    explicit Tensor(iseq n) : TensorBase<T>(n) {}

    template <class U> Tensor(Tensor<U, 1> const &u) noexcept : TensorBase<T>(u) {}
    template <class U> Tensor(Tensor<U, 1> &&u) noexcept : TensorBase<T>(std::move(u)) {}

    auto operator()(iseq i) const noexcept {return base_type::begin() + i;}
    auto operator()(iseq i) noexcept {return base_type::begin() + i;}

    auto operator()(span i) const noexcept {return offset_view(base_type::begin(), i);}
    auto operator()(span i) noexcept {return offset_view(base_type::begin(), i);}

    friend std::ostream & operator<<(std::ostream &os, Tensor const &t) noexcept {
        dump_os(os, view(t.begin(), t.end()));
        return os;
    }

    void resize(iseq n) {base_type::resize(n);}

    iseq size() const noexcept {return base_type::end() - base_type::begin();}
    std::array<iseq, 1> shape() const noexcept {return {size()};}
    constexpr auto strides() const noexcept {return std::array<iseq, 1>{1};}
};

/******************************************************************************************/

/// 2-dimensional tensor subview
template <class T>
class Tensor<T &, 2> {
    using iterator = decltype(declref<copy_qualifier<T, Tensor<decay<T>, 2>>>().begin());
    using const_iterator = decltype(declref<Tensor<decay<T>, 2> const>().begin());
    using base_type = TensorBase<T>;
    std::array<iseq, 2> dims;
    iterator m_data;
    iseq length;
public:
    using element_type = T &;
    using value_type = decay<T>;

    auto iter(iseq i, iseq j) const noexcept {return m_data + i * length + j;}

    constexpr auto strides() const noexcept {return std::array<iseq, 2>{length, 1};}
    auto data() const noexcept {return m_data;}

    NUPACK_REFLECT(Tensor, dims, length, m_data);

    Tensor(iterator b, iseq i, iseq j, iseq l) : dims{{i, j}}, m_data(b), length(l) {}
    Tensor(copy_qualifier<T, Tensor<T, 2>> &all) noexcept : dims(all.shape()), m_data(begin_of(all)), length(all.shape()[1]) {}

    auto operator()(iseq i, iseq j) const noexcept {return const_iterator(iter(i, j));}
    auto operator()(iseq i, iseq j) noexcept {return iter(i, j);}

    auto operator()(iseq i, span j) noexcept {return offset_view(iter(i, 0), j);}
    auto operator()(iseq i, span j) const noexcept {return offset_view(const_iterator(iter(i, 0)), j);}

    auto operator()(span i, span j) const noexcept {return Tensor<T &, 2>(iter(i.start(), j.start()), len(i), len(j), length);}
    /// TODO: improve formatting here
    friend std::ostream & operator<<(std::ostream &os, Tensor const &t) {
        std::ios::fmtflags f(os.flags());
        os << std::setprecision(6);
        if (product(t.dims)) for (auto i : range(t.dims[0])) {
            os << "\n";
            for (auto j : range(t.dims[1])) dump_os(os, std::setw(11), *t(i, j), ", ");
        }
        os.flags(f);
        return os;
    }
    /// Make a copy of the subview contents
    auto write(span i, span j) const {NUPACK_ERROR("should not be used");}
    /// Read in values from another matrix - R should be output of write()
    template <class R>
    void read(span i, span j, R const &r) {NUPACK_ERROR("should not be used");}

    auto shape() const noexcept {return dims;}
    auto size() const noexcept {return dims[0];}
};

template <class T, class U>
void copy_tensor_block(Tensor<T, 2> const &from, span i, span j, Tensor<U, 2> &to, span k, span l) {
    NUPACK_REQUIRE(i.stop(), <=, from.shape()[0]);
    NUPACK_REQUIRE(j.stop(), <=, from.shape()[1]);
    NUPACK_REQUIRE(k.stop(), <=, to.shape()[0]);
    NUPACK_REQUIRE(l.stop(), <=, to.shape()[1]);
    NUPACK_REQUIRE(i.start(), <=, i.stop());
    NUPACK_REQUIRE(j.start(), <=, j.stop());
    NUPACK_REQUIRE(k.start(), <=, k.stop());
    NUPACK_REQUIRE(l.start(), <=, l.stop());
    zip(i, k, [&](auto a, auto b) {std::copy(from.iter(a, j.start()), from.iter(a, j.stop()), to.iter(b, l.start()));});
}

/******************************************************************************************/

/// 2D tensor (not subview)
template <class T>
class Tensor<T const, 2>;

template <class T>
class Tensor<T, 2> : public TensorBase<T>, MemberComparable {
    using base_type = TensorBase<T>;
    std::array<iseq, 2> dims = {0, 0};


public:
    NUPACK_EXTEND_REFLECT(Tensor, base_type, dims);

    auto iter(iseq i, iseq j) noexcept {return base_type::begin() + i * dims[1] + j;}
    auto iter(iseq i, iseq j) const noexcept {return base_type::begin() + i * dims[1] + j;}

    Tensor() = default;

    template <class U>
    Tensor(Tensor<U, 2> const &u) : base_type(u), dims(u.shape()) {}

    template <class U>
    explicit Tensor(iseq m, iseq n, U t) : TensorBase<T>(m * n, t), dims{{m, n}} {}

    explicit Tensor(iseq m, iseq n) : TensorBase<T>(m * n), dims{{m, n}} {}

    explicit Tensor(typename TensorBase<T>::data_type &&data, iseq m, iseq n) noexcept : TensorBase<T>(std::move(data)), dims{{m, n}} {}

    auto operator()(iseq i, iseq j) const noexcept {return base_type::begin() + (i * dims[1] + j);}
    auto operator()(iseq i, iseq j) noexcept {return base_type::begin() + (i * dims[1] + j);}

    auto operator()(iseq i, span j) const noexcept {return offset_view(iter(i, 0), j);}
    auto operator()(iseq i, span j) noexcept {return offset_view(iter(i, 0), j);}

    auto operator()(span i, span j) const noexcept {return Tensor<T const &, 2>(iter(i.start(), j.start()), len(i), len(j), dims[1]);}
    auto operator()(span i, span j) noexcept {return Tensor<T &, 2>(iter(i.start(), j.start()), len(i), len(j), dims[1]);}
    /// TODO: improve formatting here
    friend std::ostream & operator<<(std::ostream &os, Tensor const &t) {
        std::ios::fmtflags f(os.flags());
        os << std::setprecision(6);
        if (product(t.dims)) for (auto i : range(t.dims[0])) {
            os << "\n";
            for (auto j : range(t.dims[1])) dump_os(os, std::setw(11), *t(i, j), ", ");
        }
        os.flags(f);
        return os;
    }
    /// Make a copy of some subblock of the matrix
    auto write(span i, span j) const {
        NUPACK_REQUIRE(i.start(), <=, i.stop());
        NUPACK_REQUIRE(j.start(), <=, j.stop());
        NUPACK_REQUIRE(i.stop(), <=, dims[0]);
        NUPACK_REQUIRE(j.stop(), <=, dims[1]);

        typename Tensor<T, 2>::data_type data;
        for (auto a : i) base_type::copy_out(iter(a, j.start()), iter(a, j.stop()), data);
        return Tensor<T, 2>{std::move(data), len(i), len(j)};
    }

    /// Read values into some subblock of the matrix
    template <class M>
    void read(span i, span j, M const &m) {
        NUPACK_REQUIRE(len(i) * len(j), ==, product(m.shape()));
        auto in = begin_of(m);
        for (auto a : i) {base_type::read_span(in, in + len(j), iter(a, j.start())); in += len(j);}
        NUPACK_ASSERT(in == end_of(m), i, j);
    }

    template <class V>
    void fill_subcolumn(iseq i, span j, V v) noexcept {base_type::fill(v, {i * dims[1] + j.start(), i * dims[1] + j.stop()});}

    void resize(iseq m, iseq n) {
        dims = {m, n};
        base_type::resize(m * n);
    }

    auto shape() const noexcept {return dims;}
    constexpr auto strides() const noexcept {return std::array<iseq, 2>{dims[1], 1};}
    auto size() const noexcept {return dims[0];}

    template <class P>
    auto indices_of(P && p) const noexcept {return std::array<iseq, 2>{};}

    auto indices_of(typename base_type::const_iterator t) const noexcept {
        return unpack_as<std::array<iseq, 2>>(div<false>(t - base_type::data(), dims[1]));
    }

    auto indices_of(typename base_type::iterator t) const noexcept {
        return unpack_as<std::array<iseq, 2>>(div<false>(t - base_type::data(), dims[1]));
    }

    template <class P>
    constexpr bool has(P const &p) const noexcept {return false;}

    auto has(typename base_type::const_iterator p) const noexcept {return base_type::data() <= p && p < base_type::data() + product(dims);}
    auto has(typename base_type::iterator p) const noexcept {return base_type::data() <= p && p < base_type::data() + product(dims);}
};


/******************************************************************************************/

}}
