#pragma once
#include "DesignComponents.h"
#include "OutputResult.h"

namespace nupack { namespace newdesign {


/**
 * @brief A thin wrapper associating a Sequence with its evaluated/estimated
 *     Defect (at any level in the design).
 * @details Allows functions in the Designer to return a Sequence and its
 *     Defect together.
 */
struct Result {
    Sequence sequence;
    vec<Defect> defects;
    vec<real> weights;
    mutable SingleResult evaluated;

    NUPACK_REFLECT(Result, sequence, defects, weights, evaluated);

    Result() = default;

    /** @brief forward to first Defect::total */
    Defect defect(uint i=0) const;
    Defect weighted_defect(uint i) const;
    vec<Defect> weighted_defects() const;

    vec<real> totals() const;
    vec<real> weighted_totals() const;

    real total() const;
    real weighted_total() const;

    SingleResult const & full_evaluation(Designer const &) const;

    bool operator<(Result const &other) const;
    bool operator>(Result const &other) const;
    bool operator==(Result const &other) const;
    bool operator!=(Result const &other) const;
};

static Result const inf_result{{}, {defect_vec(1, std::make_pair(uint(0), std::numeric_limits<real>::infinity()))}, {1.0}, {}};

vec<uint> first_defect_sample(Result const &res, uint num=1);
vec<uint> stochastic_hierarchical_sample(Result const &res, uint num=1);
vec<uint> scalarized_sample(Result const &res, uint num=1);
vec<uint> uniform_sample(Result const &res, uint num=1);


/**
 * @brief a structure to hold all of the best encountered sequences and defects at all levels of the design algorithm
 *
 */
template <class T>
struct DesignState {
    T dfault;

    T full;
    vec<T> forest {};
    T leaf_opt;
    T leaf_mut;

    DesignState() = default;
    DesignState(T in) : dfault(in), full(in), leaf_opt(in), leaf_mut(in) {}

    void reset(T &t) { t = dfault; }
    void reset(vec<T> &t) { t.clear(); }

    template <class ...Ts, NUPACK_IF(sizeof...(Ts) >= 1 && all_of<(is_t<Ts, T> || is_t<Ts, vec<T>>)...>)>
    void reset(Ts &...ts) { (reset(ts), ...); }
    void reset() { reset(full, leaf_opt, leaf_mut, forest); }
    void resize_forest(uint new_size) { forest.resize(new_size, dfault); }

    NUPACK_REFLECT(DesignState, dfault, full, forest, leaf_opt, leaf_mut);
};


using ResultState = DesignState<Result>;

} // newdesign
} // nupack
