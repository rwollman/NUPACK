/**
 * @brief Equilibrium toolkit functionality
 *
 * @file Equilibrate.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../types/Matrix.h"

namespace nupack::concentration {

enum class Method : uint {fit, cd, uniform, given, absls, nnls};

struct Options {
    usize max_iters = 1e4;
    real tolerance = 1.e-8;
    real delta_min = 1.e-12;
    real delta_max = 1000.0;
    bool orthogonalize = true; //< takes care of cases where # strands > # complexes
    Method method = Method::cd; //< usually the non-uniform method is better
    NUPACK_REFLECT(Options, max_iters, tolerance, delta_min, delta_max, method, orthogonalize);
};

template <class T>
struct Output {
    Col<T> solution, dual_solution;
    real objective, error;
    usize iters = 0;
    bool converged = false;
    NUPACK_REFLECT(Output, solution, dual_solution, objective, error, iters, converged);
};

/**
 * @brief Solve equilibrium concentrations
 * @param A Coefficient matrix (complexes, strands)
 * @param logb initial LOG of STRAND concentrations
 * @param q log partition functions
 * @param ops solving options
 * @return V equilibrated concentrations
 */
Output<real> equilibrate(Mat<real> const &A, Col<real> logb, Col<real> const &q, Options const &ops={});

/**
 * @brief Solve equilibrium concentrations for complexes
 * indices: list of ordered indices of strands
 * logq: list of log partition functions (distinguishable)
 * x0: strand concentrations
 */
Output<real> solve_complexes(vec<small_vec<uint>> const &indices, Col<real> logq, Col<real> x0, Options const &ops, bool rotational_correction=true, bool as_strands=true);

}
