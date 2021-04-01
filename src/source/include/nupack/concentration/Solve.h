/**
 * @brief Equilibrium toolkit algorithms - usually it will suffice to include Equilibrate.h
 *
 * @file Concentration.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Equilibrate.h"
#include "../common/Random.h"
#include "../math/BoundSolve.h"
#include "../iteration/Search.h"
#include "../reflect/SerializeMatrix.h"
#include "../reflect/Serialize.h"
#include <spdlog/spdlog.h>
#include <fmt/ostream.h>

namespace nupack::concentration {

/******************************************************************************************/

// where x is a vector and H is symmetric, return A.T * diag(x) * A
// this is probably the rate limiting step, complexity is (# strands)^2 (# complexes)
template <class M, class V>
void symmetric_mdm(M &H, M const &A, V const &x) {
    H.set_size(A.n_cols, A.n_cols);
    for (auto i : range(H.n_rows)) for (auto j : range(i+1))
        H.at(j, i) = H.at(i, j) = accu(A.col(i) % x % A.col(j));
}


/******************************************************************************************/

// x_{c, I} = \sum_j (A_{j, c} - A_{j, c_i}) \lambda_j + q_c - q_{c_i}  \\
// y_{c, I} = A_{c, I} e^{x_c, I}  \\
// \log b - \log a - q = \hat{A}^T \lambda + \log ( 1^T Y ) + \frac{ (\lambda - \lambda^0)^T (A^T Y - \hat{A} d(1^T Y)) } { 1^T Y  }

template <class T, class M>
Output<T> equilibrate_lse(M const &A, Col<T> const &logb, Col<T> const &logq, Options ops) {
    NUPACK_ASSERT(logb.is_finite(), logb, "Strand concentrations should be finite");
    NUPACK_ASSERT(logq.is_finite(), logq, "Partition functions should be finite");
    NUPACK_REQUIRE(A.min(), >=, 0, A, "Coefficient matrix should be non-negative");

    Output<T> out;
    auto &y = out.dual_solution = la::solve(A, la::log(A * la::exp(logb)) - logq);
    auto &x = out.solution;
    la::uvec c(logb.n_cols);
    M X(A.n_rows, A.n_cols), Y, Alog = la::log(A);
    Mat<T> G;
    Col<T> rhs, error, xc, shift;

    for (out.iters = 0; out.iters != ops.max_iters; ++out.iters) {
        x = A * y + logq;
        c = la::index_max(Alog.each_col() + x, 0).t();
        xc = -(A.rows(c) * y + logq(c));
        X.each_col() = x;
        X.each_row() += xc.t();
        Y = A % la::exp(X);

        rhs = la::sum(Y, 0).t();
        error = logb + xc - la::log(rhs);
        // print(out.iters, la::max(la::abs(error)), json(c));
        if (la::max(la::abs(error)) < ops.tolerance) {out.converged = true; break;}

        G = Y.t() * A;
        rhs %= error;
        y += la::solve(G, rhs); // la::solve_opts::refine
        x = A * y + logq;
    }
    x = la::exp(x);
    return out;
}

/******************************************************************************************/

// This one converges nicely but isn't very guarded against exponentiation overflow
template <class T, class M>
Output<T> equilibrate_cd(M const &A, Col<T> const &logb, Col<T> const &logq, Options ops) {
    Output<T> out;
    auto &y = out.dual_solution = la::solve(A, la::log(A * la::exp(logb)) - logq);

    Mat<T> H, V, AV;
    Col<T> e, vy, logx = A * y + logq;
    auto &x = out.solution = la::exp(logx);

    Col<T> const x0 = la::solve(A.t(), la::exp(logb));

    for (out.iters = 0; out.iters != ops.max_iters; ++out.iters) {
        if (out.iters % logb.n_rows == 0) { // O(m n^2), done every n is O(m n)
            symmetric_mdm(H, A, x);
            NUPACK_ASSERT(la::eig_sym(e, V, H), "eigendecomposition failed", H, x);
            // V.clean(1e-12); // could remove small entries just in case there's noise
            AV = A * V;
            vy = V.t() * y; // thus A y == AV vy
        }

        for (auto i : indices(logb)) { // O(m n) overall
            T const s = la::dot(AV.col(i), x), // O(m)
                   s0 = la::dot(AV.col(i), x0), // O(m)
                    h = la::accu(AV.col(i) % AV.col(i) % x); // O(m)
            if (h == 0) continue;
            T shift = (s0 - s) / h;
            NUPACK_ASSERT(std::isfinite(shift), shift, y(i), s, s0, s0 - s, h)

            if (std::abs(shift) > 16) { // exp(16) ~ 1e7, no reason magnitude should need to fluctuate that much
                // print("bad condition", shift, s, s0, s0 - s, h);
                shift = std::copysign(T(16), shift);
            }

            vy(i) += shift;
            logx += shift * AV.col(i); // O(m)
            x = la::exp(logx); // O(m)
            NUPACK_ASSERT(x.is_finite(), x, logx, vy, shift, s, s0, h, AV.col(i));
        }

        y = V * vy;
        if (la::max(la::abs(la::log(A.t() * x) - logb)) < ops.tolerance) {out.converged = true; break;}
    }

    return out;
}

/******************************************************************************************/

/// Using linear solver, gradient, Hessian, delta, and minimum delta return the dogleg direction to go to
template <class V, class M, class T>
V find_direction(V const &grad, M const &Hess, T const delta, T min_delta) {
    // Probably could optimize the order of Newton, Cauchy evaluations below
    // Calculate Newton step with a SPD, singular accepting solver
    V newt = -grad;
    bool const newton = la::solve(newt, Hess, newt,
#if (ARMA_VERSION_MAJOR >= 9) && (ARMA_VERSION_MINOR >= 500)
        la::solve_opts::fast + la::solve_opts::allow_ugly + la::solve_opts::likely_sympd
#else
// #       warning "Armadillo 9.500 or newer is recommended."
        la::solve_opts::fast + la::solve_opts::allow_ugly
#endif
    );
    T const newt_norm = norm(newt);

    // Take Newton if we are inside the minimum trust region
    if (newton && newt_norm > 0 && (delta < min_delta || newt_norm < delta))
        return newt;

    // Calculate Cauchy step
    V cauchy = grad / norm(grad);
    cauchy = grad * (-1 / dot(cauchy, Hess * cauchy));
    T const cauchy_norm = norm(cauchy);

    // Take Cauchy if Newton failed or we are outside the trust region
    if (!newton || newt_norm == 0 || !std::isfinite(sq(newt_norm)) || cauchy_norm > delta)
        return sqrt(delta / cauchy_norm) * cauchy;

    // Dogleg - take this if we are in intermediate region
    T const newt_cauchy = dot(newt, cauchy);

    auto q = quadratic_solve<T>(sq(newt_norm) + sq(cauchy_norm) - 2 * newt_cauchy,
                                2 * (newt_cauchy - sq(cauchy_norm)),
                                sq(cauchy_norm) - sq(delta));

    auto beta = std::min(q.first, q.second, less_abs); // choose correct root for mixing coefficient
    NUPACK_REQUIRE(abs(beta), <=, 1, beta, newt_norm, cauchy_norm, newt_cauchy);

    if (beta < 0) return (beta + 1) * cauchy;
    else return (1 - beta) * cauchy + beta * newt;
}

/******************************************************************************************/

template <class V, class P, class O>
struct DualSystem {
    using gradient = V;
    using hessian = Mat<value_type_of<V>>;

    template <class V_, class P_, class O_>
    DualSystem(V_ &&v, P_ &&p, O_ &&o) : dual(fw<V_>(v)), to_primal(fw<P_>(p)), objective_function(fw<O_>(o)) {
        to_primal(primal, dual);
        objective = objective_function(primal, dual);
    }

    V primal, dual;
    real objective;

    P to_primal; // (primal, dual) -> primal
    O objective_function;  // (primal, dual) -> float

    NUPACK_REFLECT(DualSystem, primal, dual, objective, to_primal, objective_function);

    // Take another system and set self to its dual value plus a shift
    template <class T> void update(DualSystem const &s, T &&shift) {
        dual = s.dual + fw<T>(shift);
        to_primal(primal, dual);
        objective = objective_function(primal, dual);
    }

    void swap(DualSystem &x) {
        swap_all(std::tie(primal, dual, objective), std::tie(x.primal, x.dual, x.objective));
    }
};

template <class ...Ts> auto dual_system(Ts &&...ts) {return DualSystem<decay<Ts>...>(fw<Ts>(ts)...);}

/******************************************************************************************/

/// Converge a system using a trust-region method
/// Provide initial system and delta, functions to compute gradient, offset, direction, convergence, and radius adjustment
template <class Sys, class G, class H, class D, class C, class A>
Sys trust_region(Sys s, G &&gradient, H &&hessian, D &&direction, C &&condition, A &&adjust_delta, real delta) {
    auto s2 = s;
    typename Sys::gradient grad, p;
    typename Sys::hessian Hess;
    for (std::size_t iter = 0; true; ++iter) {
        throw_if_signal();
        gradient(grad, s);
        hessian(Hess, s);
        if (condition(s, grad, Hess)) break;
        p = direction(grad, Hess, delta);
        s2.update(s, p);

        // rho measures actual improvement divided by expected improvement (if exact, = 1)
        auto const expected = dot(grad, p) - dot(p, Hess * p) / 2; // should almost always be < 0
        auto const rho = (s2.objective - s.objective) / expected; // better for rho to be >0

        NUPACK_ASSERT(!std::isnan(rho) && std::isfinite(delta), "error in trust region solver",
            iter, delta, rho, p, grad, Hess, s.objective, s.primal, s.dual, s2.objective, s2.primal, s2.dual, expected);

        // Possibly accept the new values. Also, adjust the trust region size
        delta = adjust_delta(delta, rho);
        if (s2.objective <= s.objective) s.swap(s2);

        // Check validity of current solution
        NUPACK_ASSERT(s.primal.is_finite() && s.dual.is_finite() && std::isfinite(s.objective),
            "trust region solver encountered non-finite value",
            iter, delta, rho, p, grad, Hess, s.objective, s.primal, s.dual, s2.objective, s2.primal, s2.dual, expected);
    }
    return s;
}

/******************************************************************************************/

template <class M, class V>
V initial_dual_guess(Method in, M const &A, V const &x0, V const &q, V const &rhs) {
    V c(len(x0));

    // Initial guess that is used by Concentration -- uses q plus a uniform guess for x0
    switch (in) {
        // Initial guess using input concentrations
        case Method::given: c = la::log(x0);
        // Initial guess using absolute value of least squares
        case Method::absls: c = la::log(A * la::abs(solve(A.t() * A, rhs)));
        // Initial guess using NNLS
        case Method::nnls: c = la::log(std::get<0>(bound_least_squares(la::eval(A.t()),
            static_cast<M const &>(rhs), ScalarBound(0, real(inf)), AlternatingOptions(5000))));
        default: c.fill(1);
    }

    // Get rid of any NaN by multiplying by min^2/max, where min and max are of finite values of c
    real bump = 1;
    for (auto const x : c) if (std::isfinite(x)) bump = std::min<real>(x, bump);
    for (auto &x : c) if (!std::isfinite(x)) x = bump;

    // Weight lower free energy complexes more for a least squares initial guess
    V const weight = arma::exp(q - q.max());
    M const AwA = A.t() * la::diagmat(weight) * A;

    return la::solve(AwA, A.t() * (weight % (c - q)));
}


/**
 * @brief Solve equilibrium concentrations
 * The objective is O(y) = 1^T exp(A y + q) + y^T A^T x0
 * @param A Coefficient matrix (complexes, strands)
 * @param x0 initial concentrations of complexes
 * @param q log partition functions
 * @param ops solving options
 * @return V equilibrated concentrations
 */
template <class T>
Output<T> equilibrate_gradient(Mat<T> A, Col<T> const &x0, Col<T> const &q, Options const &ops) {
    static auto logger = spdlog::get("concentration");
    NUPACK_ALL_EQUAL("Inconsistent number of complexes", x0.n_rows, q.n_rows, A.n_rows);

    if (logger) logger->info("equilibrate started.");
    if (logger) logger->info("A rows: {}, A columns: {}", A.n_rows, A.n_cols);

    if (!len(x0)) return {x0, 0, 0, true};

    Mat<T> orig_A, orth_A;
    if (ops.orthogonalize) {
        orig_A = A;
        orth_A = orth(A.t());
        A = A * orth_A;
    }

    // Direction finder with minimum radius delta
    auto direction = [dmin=ops.delta_min](auto const &...ts) {
        auto p = find_direction(ts..., dmin);
        // It is unknown why this can happen, but it does in some irreproducible design cases
        for (auto &x : p) if (!std::isfinite(x)) x = 0;
        return p;
    };

    usize n = 0; // number of condition checks (= 1 + number of iterations)
    real error;
    bool good = false;
    // Negative total concentrations of each strand type in the chosen basis
    Col<T> const rhs = -(A.t() * x0);
    // Negative total concentrations of each strand type in the normal basis
    Col<T> const normalization = 1 / ((ops.orthogonalize ? orig_A : A).t() * x0);
    NUPACK_ASSERT(normalization.is_finite() && normalization.min() > 0, normalization);
    // Convergence criterion based on gradient norm
    auto condition = [&] (auto const &sys, auto const &grad, auto const &...) {
        if (!ops.orthogonalize) error = la::max(la::abs(grad) % normalization);
        else error = la::max(la::abs(orth_A * grad) % normalization);
        good = error < ops.tolerance;
        return (ops.max_iters < ++n) || good;
    };
    // Function to calculate primal vector from previous primal and dual
    auto primal = [&](auto &x, auto const &y) {
        // In the future, 1e100 is a bit arbitrary.
        // This is done to prevent infinity from occurring during intermediate solution
        // It would be better to constrain the dual vector y to be in a reasonable regime.
        x = la::clamp(la::exp(A * y + q), std::numeric_limits<real>::min(), 1e100);
        NUPACK_ASSERT(x.is_finite(), x, A, y, q, A * y + q, n);
    };
    // Function to calculate objective from primal and dual
    auto objective = [&](auto const &x, auto const &y) {return accu(x) + dot(y, rhs);};

    auto sys = dual_system(initial_dual_guess(ops.method, A, x0, q, rhs), primal, objective);

    // Function to calculate gradient in dual space from system
    auto gradient = [&](auto &grad, auto const &s) {
        grad = A.t() * s.primal + rhs;
        NUPACK_ASSERT(grad.is_finite(), s.primal, n);
    };
    // Function to calculate gradient in dual space from system
    // this is probably the rate limiting step, complexity is (# strands)^2 (# complexes)
    auto hessian = [&](auto &H, auto const &s) {
        symmetric_mdm(H, A, s.primal);
        NUPACK_ASSERT(H.is_finite(), s.primal, n);
    }; // same as H = A.t() * la::diagmat(s.primal) * A;
    // Function to adjust trust region radius and return if the new values should be accepted
    auto adjust_delta = [=](auto const delta, auto const rho) {
        if (delta <= ops.delta_min) return delta;
        if (rho > 0.75) return std::min(2 * delta, ops.delta_max);
        if (rho < 0.25) return delta / 4;
        return delta;
    };
    // Run trust region with initial radius 1000
    auto ret = trust_region(std::move(sys), gradient, hessian, direction, condition, adjust_delta, ops.delta_max);

    if (!good && logger)
        logger->error("A:\n{}\northo A:\n{}\nx0:\n{}\ng:\n{}\nx:\n{}\n", orig_A, A, x0, q, ret.primal);

    if (logger) logger->info("equilibrate finished. Number of iterations: {}", n - 1);

    // Undo the threshold that was put on in primal()
    ret.primal.replace(std::numeric_limits<real>::min(), 0);

    if (ops.orthogonalize) ret.dual = orth_A * ret.dual;
    return {std::move(ret.primal), std::move(ret.dual), ret.objective, error, n - 1, good};
}

/******************************************************************************************/

template <class T>
Output<T> equilibrate_finite(Mat<T> const &A, Col<T> const &logb, Col<T> const &logq, Options const &ops) {
    switch (ops.method) {
        case Method::cd: return equilibrate_cd<T>(A, logb, logq, ops);
        case Method::fit: return equilibrate_lse<T>(A, logb, logq, ops);
        default: return equilibrate_gradient(A, Col<T>(la::solve(A.t(), la::exp(logb))), logq, ops);
    }
}

/******************************************************************************************/

}
