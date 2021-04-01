#include <nupack/concentration/Solve.h>

namespace nupack::concentration {

/******************************************************************************************/

template <class V, class T, class F=Identity>
auto maxmap(V &&v, T const &init, F &&f={}) {
    std::decay_t<decltype(f(*begin_of(v)))> out = init;
    for (auto &&x : v) {auto tmp = f(x); if (tmp > out) out = tmp;}
    return out;
}

/******************************************************************************************/

struct LexicographicalCompare {
    template <class T, class U>
    bool operator()(T const &t, U const &u) const {
        return std::lexicographical_compare(begin_of(t), end_of(t), begin_of(u), end_of(u));
    }
};

Output<real> equilibrate(Mat<real> const &A, Col<real> logb, Col<real> const &logq, Options const &ops) {
    if (logq.has_nan()) throw std::domain_error("Input log Q contains NaN");
    if (logq.max() == real(*inf)) throw std::domain_error("Input log Q contains +inf");

    // Reduce complexes of the same composition
    std::map<Col<real>, real, LexicographicalCompare> unique;
    izip(logq, [&](auto i, auto logq) {
        if (logq == -real(*inf)) return;
        Col<real> a = A.row(i).t();
        auto &p = unique.try_emplace(std::move(a), -real(*inf)).first->second;
        p = log_sum_exp(p, logq);
    });

    Mat<real> A2(unique.size(), A.n_cols);
    Col<real> logq2(unique.size());
    izip(unique, [&](auto i, auto const &p) {
        A2.row(i) = p.first.t();
        logq2(i) = p.second;
    });

    auto out = equilibrate_finite<real>(A2, logb, logq2, ops);
    out.solution = la::exp(A * out.dual_solution + logq);
    return out;
}

/******************************************************************************************/

Output<real> solve_complexes(vec<small_vec<uint>> const &indices, Col<real> logq, Col<real> x0, Options const &ops, bool rotational_correction, bool as_strands) {
    NUPACK_ALL_EQUAL("Inconsistent number of complexes", len(indices), len(logq));
    Mat<real> A(len(logq), maxmap(indices, 0, [](auto const &p) {return maxmap(p, 0) + 1;}), arma::fill::zeros);

    izip(indices, [&](auto c, auto const &x) {
        for (auto s : x) A(c, s) += 1;
        if (rotational_correction) logq(c) -= std::log(real(rotational_symmetry(x)));
    });

    NUPACK_REQUIRE(x0.min(), >, 0, "All concentrations should be positive");

    if (as_strands) {
        NUPACK_REQUIRE(len(x0), ==, A.n_cols, "Incorrect number of concentrations given");
        x0 = la::log(std::move(x0));
    } else {
        NUPACK_REQUIRE(len(x0), ==, A.n_rows, "Incorrect number of concentrations given");
        x0 = la::log(A.t() * x0);
    }

    return equilibrate(std::move(A), std::move(x0), std::move(logq), ops);
}

/******************************************************************************************/

}
