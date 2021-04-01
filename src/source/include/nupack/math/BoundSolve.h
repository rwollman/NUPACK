
#pragma once
#include "../types/Matrix.h"

namespace nupack {

/******************************************************************************************/

struct ScalarBound {
    double minimum, maximum, regularization;
    NUPACK_REFLECT(ScalarBound, minimum, maximum, regularization);

    ScalarBound(double min=0, double max=*inf, double reg=0) : minimum(min), maximum(max), regularization(reg) {
        if (maximum < minimum) throw std::invalid_argument("minimum must not be greater than maximum.");
    }

    template <class T>
    auto operator()(T const &value, Ignore={}, Ignore={}) const {
        return std::clamp<T>(value, minimum, maximum);
    }

    template <class B>
    auto regularize(B const &b) const {return regularization * arma::norm(b);}
};

/******************************************************************************************/

template <class T>
struct VectorBound {
    Mat<T> bounds; // (2, N)
    Col<T> regularization;
    NUPACK_REFLECT(VectorBound, bounds, regularization);

    VectorBound() = default;
    VectorBound(Mat<T> b, Col<T> reg) : bounds(std::move(b)), regularization(std::move(reg)) {
        if (bounds.n_rows != 2) throw std::invalid_argument("Number of rows must be 2");
        NUPACK_ALL_EQUAL("Inconsistent dimension", bounds.n_cols, regularization.n_rows);
    }

    auto operator()(T const &value, uint i, Ignore={}) const {
        return std::clamp(value, bounds(0, i), bounds(1, i));
    }

    template <class B>
    auto regularize(B const &b) const {return regularization % b;}
};


/******************************************************************************************/

struct AlternatingOptions {
    std::size_t iters;
    double tolerance;
    bool warm_start;
    NUPACK_REFLECT(AlternatingOptions, iters, tolerance, warm_start);

    AlternatingOptions(std::size_t n=5000, double tol=1e-8, bool warm=false) :
        iters(n), tolerance(tol), warm_start(warm) {
            if (iters == 0) throw std::invalid_argument("number of iterations should be greater than 0");
            if (tolerance <= 0) throw std::invalid_argument("tolerance must be positive.");
        }


    bool operator()(Ignore, Ignore, double obj0, double obj) const {
        // print(obj0, obj, options.tolerance, (obj0 - obj) / options.tolerance);
        return !(obj0 - obj > sq(tolerance) * std::abs(obj));
    }
};


/******************************************************************************************/

template <class T>
struct AlternatingResult {
    uint unconverged = 0, iters = 0;
    T objective = 0;
    NUPACK_REFLECT(AlternatingResult, unconverged, iters, objective);
};

// solve A X = G
// then
// A^T A X = A^T G
// (A^T A) A^T = A^T G

/******************************************************************************************/

template <class A>
struct ClampSolver {
    using value_type = value_type_of<std::decay_t<A>>;
    Col<value_type> u, m, a_diag;
    A a;
    AlternatingOptions options;

    ClampSolver(A a0, AlternatingOptions const &ops) : a_diag(a0.diag()), a(static_cast<A &&>(a0)), options{ops} {}

    // solve $\min_x x^T A x - 2 b^T x$ subject to clamp(x)
    // if D is the user domain, clamp(k, x0) should return \min_{x \in D} (x - x0)^2
    template <class X, class B, class F>
    std::tuple<value_type, uint, bool> operator()(X &&x, B const &b, F const &clamp) {
        NUPACK_ALL_EQUAL("Inconsistent dimensions", a.n_rows, a.n_cols, b.n_rows, x.n_rows);
        if (options.warm_start && x.n_rows == u.n_rows) x = u; // use previous guess

        m = b - a.t() * x - clamp.regularize(b);
        real objective = -(arma::dot(m, x) + arma::dot(b, x)); // x^T A x - 2 x^T b
        uint z = 0;
        bool conv = false;
        for (; z != options.iters; ++z) {
            u = x;
            auto obj = objective;
            for (uint k = 0; k != a.n_cols; ++k) {
                if (a_diag[k] != 0) { // else x(k) will just be left at 0
                    value_type const t = clamp(x[k] + m[k] / a_diag[k], k); // the latter is the analytic unconstrained argmin
                    if (t != x[k]) m += (x[k] - t) * a.col(k); // old - new
                    x[k] = t;
                } else if (x[k] != 0) { // fix up residual if x(k) not already at 0 (e.g. from initialization)
                    m += x[k] * a.col(k);
                    x[k] = 0;
                }
            }
            objective = -(arma::dot(m, x) + arma::dot(b, x));
            if ((conv = options(u, x, obj, objective))) break;
        }
        return {objective, z, !conv};
    }
};

/******************************************************************************************/

template <class A, class C>
struct LogNormalSolver : ClampSolver<A> {
    using base = ClampSolver<A>;
    using base::a, base::a_diag, base::options, base::convergence, base::m, base::u;
    using value_type = typename base::value_type;
    C c;
    Col<value_type> d, n, c_diag;

    LogNormalSolver(A a0, C c0, AlternatingOptions const &ops) : base(static_cast<A &&>(a0), ops), c(static_cast<C &&>(c0)), c_diag(c0.diag()) {}

    // solve $\min_x (e^x)^T A e^x - 2 b^T e^x + x^T C x - 2 d^T x$
    template <class X, class B, class D>
    std::tuple<value_type, uint, bool> operator()(X &&x, B const &b, D const &d) {
        NUPACK_ALL_EQUAL("Inconsistent dimensions", a.n_rows, a.n_cols, b.n_rows, x.n_rows);
        if (options.warm_start) x = u; // use previous guess

        m = b - a.t() * arma::exp(x);
        n = d - c.t() * x;

        uint z = 0;
        for (; z == 0 || (z < options.iters && convergence(u, x)); ++z) {
            u = x;
            for (uint k = 0; k != a.n_cols; ++k) {
                value_type mk = m[k], nk = n[k], xk, ek;
                value_type x0 = mk > 0 ? std::log(mk / a_diag[k]) : 0, e0 = std::exp(x0);

                // Newton solve in 1 dimension
                do {
                    xk = x0 + (e0 * (b - a * e0) - c * x0) / (c + e0 * (2 * a * e0 - b));
                    ek = std::exp(xk);
                    // scalar updates to residuals
                    mk += (ek - e0) * a_diag[k];
                    nk += (xk - x0) * c_diag[k];

                    x0 = xk;
                    e0 = ek;
                } while (xk - xk > 1e-6);

                // vector updates to stored residuals
                m += (ek - std::exp(x[k])) * a.col(k);
                n += (xk - x[k]) * c.col(k);
            }
        }
        m = arma::exp(x);
        return {arma::dot(a.t() * m, m) - 2 * arma::dot(b, m)
              + arma::dot(c.t() * x, x) - 2 * arma::dot(d, x), z != options.iters};
    }
};

/******************************************************************************************/

/*
 * NNLS modified to take A^T * A, A^T * B instead of A, B, and B is matrix instead of vector
 * The solution x is modified in place
 * Returns number of unconverged columns of B and the objective (x^T A x - 2 x^T b) summed over columns of B
 * For a least squares problem, the error is just the objective + b^T b
 * The residual norm^2 is added to in place if provided
 * For least squares, the residual norm is x^T AA x - 2 BA x, so you must add ||B||^2 if you want the true error
 * Returns the number of unconverged points and the total residual squared norm
 * Description: sequential Coordinate-wise algorithm for non-negative least square regression A x = b, s^t. x >= 0
 * Modified from: https://github.com/linxihui/Misc/blob/master/Practice/NMF/nnls.cpp
 * Reference: http://cmp.felk.cvut.cz/ftp/articles/franc/Franc-TR-2005-06.pdf
 */
template <class F, class A, class T>
AlternatingResult<T> bound_solve(Mat<T> &x, A const &a, Mat<T> const &b, F const &bound, AlternatingOptions const &ops, Col<T> *norm2=nullptr) {
	NUPACK_ALL_EQUAL("Inconsistent dimensions", b.n_cols, x.n_cols);
    ClampSolver<A const &> solve(a, ops);
    if (norm2 && norm2->n_rows != b.n_cols) {norm2->set_size(b.n_cols); norm2->zeros();}

    AlternatingResult<T> out{0, 0, 0};

    // For each column of B
    for (uint i = 0; i != b.n_cols; ++i) {
        auto [err, iters, unconv] = solve(x.col(i), b.col(i), bound);// [&](auto x, auto k) {return bound(x, k, i);});
        if (norm2) (*norm2)(i) += err;
        out.objective += err;
        if (unconv) ++out.unconverged;
        out.iters += iters;
    }
    return out;
}

/******************************************************************************************/

// // 2 bounds, each bound either non-existent, 0, or finite. gives 3 * 3 - 1 = 8 possibilities
// template <class A, class T>
// AlternatingResult<T> bound_solve(Mat<T> &x, A const &a, Mat<T> const &b, Bounds const &bound, AlternatingOptions const &ops, Col<T> *norm2=nullptr) {
//     NUPACK_REQUIRE(ops.minimum, <=, ops.maximum);
//     // Optimization for unconstrained case
//     if (ops.minimum == real(minf) && ops.maximum == real(inf)) {
//         if constexpr(la::is_sparse<A>) arma::spsolve(x, a, b);
//         else arma::solve(x, a, b);
//         return {0, 0, sq(arma::norm(a * x - b))};
//     }
//     return bound_solve(x, a, b, [min=T(ops.minimum), max=T(ops.maximum)](Ignore, Ignore, T const &x) {
//         return std::clamp(x, min, max);
//     }, ops, norm2);
// }

/******************************************************************************************/

// Non-negative least squares - just solve A^T A x = A^T b instead
template <class A, class B, class F>
auto bound_least_squares(Mat<A> const &a, Mat<B> const &b, F const &bound, AlternatingOptions const &ops={}) {
	using V = std::common_type_t<A, B>;
    Mat<V> x(a.n_cols, b.n_cols, arma::fill::zeros);
    auto res = bound_solve(x, Mat<V>(a.t() * a), Mat<V>(a.t() * std::move(b)), bound, ops);
    // auto err = arma::norm(b - a * x, "fro");
    res.objective += la::accu(b % b);
    return std::make_pair(std::move(x), std::move(res));
}

/******************************************************************************************/

// Cichocki & Phan: Algorithms for Nonnegative matrix and tensor factorization (2008)
template <class Y, class T>
auto als(Y const &y, Mat<T> a, ScalarBound const &ops) {
    Mat<T> x = arma::clamp(arma::solve(a, y), ops.minimum, ops.maximum);
    a = arma::clamp(arma::solve(x, y.t()), ops.minimum, ops.maximum);
    return std::make_pair(std::move(a), std::move(x));
}

template <class Y>
auto als(Y const &y, std::size_t n, ScalarBound const &ops) {
    return als(y, arma::randu<Mat<value_type_of<Y>>>(y.n_rows, n), ops);
}

/******************************************************************************************/

// Cichocki & Phan: Algorithms for Nonnegative matrix and tensor factorization (2008)
template <class T>
auto hals_nmf(Mat<T> const &y, std::size_t m, AlternatingOptions const &ops) {
    std::size_t iter;
    Mat<T> W, V, P, Q, A, B;

    for (iter = 0; iter != ops.iters; ++iter) {
        W = y.t() * A;
        V = A.t() * A;
        for (auto j : range(m)) {
            B.col(j) += W.col(j) - B * V.col(j);
            // clamp
        }
        P = y.t() * B;
        Q = B.t() * B;
        for (auto j : range(m)) {
            A.col(j) += W.col(j) - A * V.col(j);
            // clamp
            A.col(j) /= arma::norm(A.col(j));
        }
    }
    Col<T> w = arma::norm(B);
    B.each_col() /= w;
    return std::make_pair(A, B);
}


// // Cichocki & Phan: Algorithms for Nonnegative matrix and tensor factorization (2008)
// template <class T>
// auto symmetric_nmf(Mat<T> const &B, std::size_t m, AlternatingOptions const &ops) {
//     std::size_t iter;
//     Mat<T> W, V, A(B.n_cols, m, arma::fill::ones);

//     for (iter = 0; iter != ops.iters; ++iter) {
//         W = B.t() * A;
//         V = A.t() * A;
//         for (auto j : range(m)) {
//             A.col(j) += W.col(j) - A * V.col(j);
//             // clamp
//             c
//         }
//     }
//     return A;
// }

// tr[A A.T - B]

/******************************************************************************************/

// Cichocki & Phan: Algorithms for Nonnegative matrix and tensor factorization (2008)
template <class Y>
auto hals_ntf(Y const &y, std::size_t m, AlternatingOptions const &ops) {
    constexpr uint N = la::depth<Y>;
    using T = value_type_of<Y>;
    std::array<Mat<T>, N> U;

    Mat<T> UU, T2, T3, T1 = U[0].t() * U[0];
    for (auto n : range(1, N)) T1 %= U[n].t() * U[n];

    std::size_t iter = 0;
    for (Col<T> gamma; iter != ops.iters; ++iter) {
        UU = U[N].t() * U[N];
        gamma = UU.diag();
        for (uint n = N - 1; n--;) {
            if (!n) gamma.ones();
            NUPACK_ERROR("not implemented");
            // T2 = contract over all indices not m and not n in Y, U;
            T3 = T1 / (U[n].t() * U[n]);
            for (auto j : range(m)) {
                // U[n].col(j) = arma::clamp(gamma(j) * U[n].col(j) + T2.col(j) - U[n] * T3.col(j), ops.minimum, ops.maximum);
                if (n) U[n].col(j) /= arma::abs(arma::norm(U[n].col(j)));
            }
            T1 = T3 * (U[n].t() * U[n]);
        }
    }
    return U;
}

/******************************************************************************************/

}
