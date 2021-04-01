/** \file Matrix.h
 * @brief Declares matrix types to be used, and abstracts Eigen/Armadillo partially
 */
#pragma once
#include "../standard/Vec.h"
#include "../standard/Array.h"
#include "../reflect/Print.h"
#include "../reflect/Memory.h"
#include "../common/Error.h"

#include "../iteration/Range.h"
#include "../iteration/Patterns.h"
#include "../iteration/View.h"
#include "../algorithms/Operators.h"
#include "../algorithms/Utility.h"

#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

namespace nupack {

template <class T>
Range<T>::operator arma::span() const {return arma::span(*begin(), *end()-static_cast<T>(1u));}

template <class T> struct detail::value_type_of_t<arma::SpMat<T>> {using type = T;};

NUPACK_DETECT(has_eval, decltype(declval<T>().eval()));

/// Is an armadillo type
template <class T> static constexpr bool is_arma = arma::is_arma_type<decay<T>>::value
                                                || arma::is_arma_cube_type<decay<T>>::value
                                                || arma::is_arma_sparse_type<decay<T>>::value;

/******************************************************************************************/

namespace la {

using namespace arma;
namespace traits {}
using namespace la::traits;

namespace detail {
    template <class T, usize N> struct Dense;
    template <class T> struct Dense<T, 1> {using type = arma::Col<T>;};
    template <class T> struct Dense<T, 2> {using type = arma::Mat<T>;};
    template <class T> struct Dense<T, 3> {using type = arma::Cube<T>;};

    template <class T, class D, std::size_t ...Is>
    typename Dense<T, sizeof...(Is)>::type dense_from_data(T *t, D const &dim, bool copy, bool strict, std::index_sequence<Is...>) {
        return {t, static_cast<arma::uword>(dim[Is])..., copy, strict};
    }
}

/// Get an armadillo type for the element type and number of dimensions
template <class T, usize N>
using Dense = typename detail::Dense<T, N>::type;

/// Make an armadillo tensor from underlying mutable data
template <class T, std::size_t N, class D>
Dense<T, N> dense_from_data(T *t, D const &dim, bool copy=true, bool strict=false) {
    return detail::dense_from_data(t, dim, copy, strict, std::make_index_sequence<N>());
}

/// Make an armadillo tensor by copying immutable data
template <class T, std::size_t N, class D>
Dense<T, N> dense_from_data(T const *t, D const &dim, bool copy=true, bool strict=false) {
    return detail::dense_from_data(const_cast<T *>(t), dim, true, false, std::make_index_sequence<N>());
}

/******************************************************************************************/

/// Get the dimensions of an armdillo type as a tuple
template <class T> std::array<usize, 1> shape(arma::Col<T> const &t) {return {t.n_rows};}
template <class T> std::array<usize, 1> shape(arma::Row<T> const &t) {return {t.n_cols};}
template <class T> std::array<usize, 2> shape(arma::Mat<T> const &t) {return {t.n_rows, t.n_cols};}
template <class T> std::array<usize, 2> shape(arma::SpMat<T> const &t) {return {t.n_rows, t.n_cols};}
template <class T> std::array<usize, 3> shape(arma::Cube<T> const &t) {return {t.n_rows, t.n_cols, t.n_slices};}

template <class A, NUPACK_IF(arma::is_Col_fixed_only<A>::value)>
constexpr std::array<usize, 1> shape() {return {A::n_rows};}

template <class A, NUPACK_IF(arma::is_Row_fixed_only<A>::value)>
constexpr std::array<usize, 1> shape() {return {A::n_cols};}

template <class A, NUPACK_IF(arma::is_Mat_fixed_only<A>::value)>
constexpr std::array<usize, 2> shape() {return {A::n_rows, A::n_cols};}

template <class A>
static constexpr auto depth = std::tuple_size<decltype(shape(declval<A>()))>::value;

/******************************************************************************************/

/// Get the strides of an armdillo type: 1D types have stride 1
template <class A, NUPACK_IF(arma::is_Col<A>::value || arma::is_Row<A>::value)>
constexpr std::array<usize, 1> strides(A const *) {return {1u};}
/// Get the strides of an armdillo type: Mat is stored column-major
template <class A, NUPACK_IF(arma::is_Mat_fixed_only<A>::value)>
constexpr std::array<usize, 2> strides(A const *) {return {1u, A::n_rows};}
/// Get the strides of an armdillo type: Mat is stored column-major
template <class T>
std::array<usize, 2> strides(arma::Mat<T> const *t) {return {1u, t->n_rows};}
/// Get the strides of an armdillo type: Cube is stored as list of matrices
template <class T>
std::array<usize, 3> strides(arma::Cube<T> const *t) {return {1u, t->n_rows, t->n_rows * t->n_cols};}

/******************************************************************************************/

template <class T> constexpr auto re(T &&t) -> decltype(fw<T>(t).real()) {return fw<T>(t).real();}
template <class T> constexpr auto re(T &&t) -> decltype(arma::real(fw<T>(t))) {return arma::real(fw<T>(t));}

/******************************************************************************************/

/// Evaluate expression template if possible
template <class T, NUPACK_IF(has_eval<T &&>)>
decltype(auto) eval(T &&t) {return fw<T>(t).eval();}

template <class T, NUPACK_IF(!has_eval<T &&>)>
decltype(auto) eval(T &&t) {return fw<T>(t);}

/// Is an expression already evaluated?
template <class T> static constexpr bool is_eval = is_same<no_ref<decltype(eval(declval<T>()))>, T>;

/******************************************************************************************/

/// Type from evaluating an expression template
template <class T>
using eval_result = no_qual<decltype(eval(declval<if_t<std::is_array<T>::value, decay<T>, T>>()))>;

/// Is a sparse matrix object
template <class T>
static constexpr bool is_sparse = arma::is_arma_sparse_type<eval_result<T>>::value;

/// Is a dense matrix object
template <class T>
static constexpr bool is_dense = is_arma<eval_result<T>> && !is_sparse<T>;

/// Make a matrix of zeros
template <class M, class ...Ts, NUPACK_IF(is_arma<M>)>
M zeros(Ts &&...ts) {return M(fw<Ts>(ts)..., arma::fill::zeros);}

/******************************************************************************************/

/// Resize an armadillo object
template <class M, class ...Ts, NUPACK_IF(is_arma<M>)>
void resize(M &m, Ts &&...ts) {m.resize(fw<Ts>(ts)...);}

/******************************************************************************************/

template <class T> using Col = arma::Col<T>;
template <class T> using Row = arma::Row<T>;
template <class T> using Mat = arma::Mat<T>;
template <class T> using SpMat = arma::SpMat<T>;
template <class T> using Cube = arma::Cube<T>;

using real_col = Col<real>;
using real_row = Row<real>;
using real_mat = Mat<real>;
using real_csc = SpMat<real>;

/******************************************************************************************/

NUPACK_DETECT(is_col, void_if<arma::resolves_to_colvector<decay<T>>::value>);
NUPACK_DETECT(is_row, void_if<arma::resolves_to_rowvector<decay<T>>::value>);

/******************************************************************************************/

template <class T, NUPACK_IF(is_col<T>)> decltype(auto) to_col(T const &t) {return t;}
template <class T, NUPACK_IF(is_row<T>)> decltype(auto) to_col(T const &t) {return trans(t);}

template <class T, NUPACK_IF(is_col<T>)> decltype(auto) flip(T &&t) {return flipud(fw<T>(t));}
template <class T, NUPACK_IF(is_row<T>)> decltype(auto) flip(T &&t) {return fliplr(fw<T>(t));}

/******************************************************************************************/

// vector vector
template <class T, class U, class F, NUPACK_IF(!(is_col<T> && is_row<U>) && !(is_row<T> && is_col<U>))>
decltype(auto) align(T &&t, U &&u, F &&f) {return f(fw<T>(t), fw<U>(u));}
// vector vector.t()
template <class T, class U, class F, NUPACK_IF((is_col<T> && is_row<U>) || (is_row<T> && is_col<U>))>
decltype(auto) align(T &&t, U &&u, F &&f) {return f(fw<T>(t), fw<U>(u).t());}

/******************************************************************************************/

/// Elementwise product
template <class T, class U, NUPACK_IF(is_scalar<T> && is_scalar<U>)>
auto schur(T const &t, U const &u) {return t * u;}
/// Elementwise product
template <class T, class U, NUPACK_IF((is_scalar<T> || is_scalar<U>) && (is_arma<T> || is_arma<U>))>
auto schur(T const &t, U const &u) {return t * u;}
/// Elementwise product, have to use % here
template <class T, class U, NUPACK_IF(is_arma<T> && is_arma<U>)>
auto schur(T const &t, U const &u) {return align(t, u, modulus);}

/******************************************************************************************/

template <class M> auto fill_zero(M &m) -> decltype(m.zeros()) {return m.zeros();}
template <class M> auto n_rows(M const &m) -> decltype(m.n_rows) {return m.n_rows;}
template <class M> auto n_cols(M const &m) -> decltype(m.n_cols) {return m.n_cols;}

/******************************************************************************************/

// Generalized dot product. The first and last arguments are columns. The middle arguments are matrices.
template <class T, class ...Ts, NUPACK_IF(sizeof...(Ts) >= 1)>
auto dot(T const &t, Ts &&... ts) {return as_scalar(t.t() * (ts * ...));}

template <class T, NUPACK_IF(is_scalar<T>)>
T trans(T const &t){return t;}

/// Outer product of two vectors, any orientation accepted
template <class T1, class T2, NUPACK_IF(is_row<T1> && is_row<T2>)>
auto outer(T1 const &t1, T2 const &t2) {return eval(trans(t1) * t2);}
/// Outer product of two vectors, any orientation accepted
template <class T1, class T2, NUPACK_IF(is_row<T1> && is_col<T2>)>
auto outer(T1 const &t1, T2 const &t2) {return eval(trans(t1) * trans(t2));}
/// Outer product of two vectors, any orientation accepted
template <class T1, class T2, NUPACK_IF(is_col<T1> && is_row<T2>)>
auto outer(T1 const &t1, T2 const &t2) {return eval(t1 * t2);}
/// Outer product of two vectors, any orientation accepted
template <class T1, class T2, NUPACK_IF(is_col<T1> && is_col<T2>)>
auto outer(T1 const &t1, T2 const &t2) {return eval(t1 * trans(t2));}

/******************************************************************************************/

// Call function with arguments of matrix, vector, vector like f(o(:, j), t(:), u(j));
template <class O, class T, class U, class F>
void visit_outer(O &&o, T &&t, U &&u, F &&f) {
    for (usize j = 0; j != o.n_cols; ++j) f(o.col(j), t, u(j));
}

// o(i, j) += t(i) * t(u)
template <class O, class T, class U>
void add_outer(O &o, T const &t, U const &u) {
    visit_outer(o, t, u, [](auto &&o, auto const &t, auto u) {o += t * u;});
}

/******************************************************************************************/

template <class M> auto at(M &m, usize i, usize j) -> decltype(m(i,j)) {return m(i,j);}

/******************************************************************************************/

template <class M> auto esum(M const &m) -> decltype(arma::accu(m)) {return arma::accu(m);}

/******************************************************************************************/

template <class M> auto esum(M const &m) -> decltype(m.sum()) {return m.sum();}

template <class M> auto msum(M const &m, usize i) -> decltype(arma::sum(m, i)) {return arma::sum(m, i);}
template <class M> auto mmin(M const &m, usize i) -> decltype(arma::min(m, i)) {return arma::min(m, i);}
template <class M> auto mmax(M const &m, usize i) -> decltype(arma::max(m, i)) {return arma::max(m, i);}

/******************************************************************************************/

template <class M> auto eabs(M const &m) -> decltype(m.cwiseAbs()) {return m.cwiseAbs();}

/******************************************************************************************/

template <class M> auto eabs(M const &m) -> decltype(arma::abs(m)) {return arma::abs(m);}

/******************************************************************************************/

template <class Mat, class F> auto for_cols(Mat const &M, F &&f) {
    for (size_type_of<Mat> j = 0; j != n_cols(M); ++j) f(j);
}

template <class Mat, class F> auto for_rows(Mat const &M, F &&f) {
    for (size_type_of<Mat> i = 0; i != n_rows(M); ++i) f(i);
}

/******************************************************************************************/

template <class Mat> auto matrix_chi_squared(Mat const &x, Mat const &y) {
    Mat ret = (x - y) * (x - y) / (x + y);
    replace_if(ret, is_nan, 0);
    return esum(std::move(ret));
}

/******************************************************************************************/

template <class T>
Col<T> raveled(Mat<T> &m) {return {m.memptr(), m.n_rows * m.n_cols, false};}

template <class T>
arma::umat nonzero_indices(SpMat<T> const &A) {
    arma::umat out(2, A.n_nonzero);
    arma::uword i = 0;
    for (auto it = A.begin(); it != A.end(); ++it, ++i) {
        out(0, i) = it.row();
        out(1, i) = it.col();
    }
    return out;
}

template <class T>
Col<T> sparse_view(SpMat<T> &m) {return Col<T>(const_cast<T *>(m.values), m.n_nonzero, false);}

/******************************************************************************************/

template <class V, class M, class F>
void sparse_map(V const &idx, M &&m, F &&f) {
    for (auto i : range(idx.n_cols)) f(i, m(idx(0, i), idx(1, i)));
}

template <class V, class F>
void sparse_map(V const &idx, F &&f) {
    for (auto i : range(idx.n_cols)) f(i, idx(0, i), idx(1, i));
}

/******************************************************************************************/

template <class V>
Mat<value_type_of<value_type_of<V>>> stack_columns(V const &v) {
    if (!len(v)) return {};
    Mat<value_type_of<value_type_of<V>>> out(len(front(v)), len(v));
    izip(v, [&](auto i, auto const &c) {out.col(i) = c;});
    return out;
}

/******************************************************************************************/

struct Solver {
    string kind;

public:
    NUPACK_REFLECT(Solver, kind);

    void wait() const {}
    Solver(string k="superlu") : kind(k) {}

    template <class X, class A, class B, NUPACK_IF(is_sparse<A>)>
    void operator()(X &&x, A &&a, B &&b, real tol=1.e-14) const {x = spsolve(fw<A>(a), fw<B>(b), kind.c_str());}

    template <class X, class A, class B, NUPACK_IF(!is_sparse<A>)>
    void operator()(X &&x, A &&a, B &&b, real tol=1.e-14) const {x = solve(fw<A>(a), fw<B>(b));}
};

/******************************************************************************************/

}

using la::Row;
using la::Col;
using la::Mat;
using la::SpMat;
using la::Cube;
using la::real_mat;
using la::real_row;
using la::real_col;
using la::real_csc;

/******************************************************************************************/

NUPACK_DETECT(has_n_elem, decltype(declval<T>().n_elem));

template <class T> struct extents::length<T, void_if<is_arma<T> && !has_size<T> && has_n_elem<T>>> {
    constexpr auto operator()(T const &t) const {return t.n_elem;}
};

template <class T> struct ConvertConstant<T, size_constant<0>, void_if<is_arma<T>>> {
    template <class ...Ts>
    constexpr T operator()(Ts &&...ts) const {return T(fw<Ts>(ts)..., arma::fill::zeros);}

    static constexpr bool ok = true;
};

static_assert(ConvertConstant<arma::Mat<double>::fixed<4,4>, size_constant<0>, void>::ok, "");

template <class T> struct io::SingleLine<T, void_if<is_arma<T>>> : False {
    constexpr bool operator()(T const &) const {return false;}
};

template <class T> struct memory::impl<T, void_if<is_arma<T> && is_same<la::eval_result<T>, T> && can_construct<T>>> {
    std::size_t operator()(T const &t) const {return sizeof(T) + t.n_elem * sizeof(typename T::elem_type);}
    void erase(T &t) const {T t_; swap(t, t_);}
};

template <class T> struct memory::impl<T, void_if<is_arma<T> && !(is_same<la::eval_result<T>, T> && can_construct<T>)>> {
    std::size_t operator()(T const &) const {return sizeof(T);}
    void erase(T &) const {}
};

/******************************************************************************************/

}
