/** \file ComplexSet.h
 * @brief Contains ComplexSet, a class representing all strands/complexes present in a State
 */
#pragma once

#include "../model/Move.h"

#include "../iteration/Search.h"
#include "../algorithms/Operators.h"
#include "../algorithms/Utility.h"

#include "../types/Fenwick.h"
#include "../types/IO.h"

namespace nupack {

/******************************************************************************************/

struct ComplexSet : ConstIndexable<ComplexSet> {
    using Index = iseq;
    using Indices = vec<iseq>;
    using Mat = BaseMat<real>;

    struct StrandData {
        Index x, pos, loop;
        NUPACK_REFLECT(StrandData, x, pos, loop);
        friend std::ostream & operator<<(std::ostream &os, StrandData const &t) {
            dump_os(os, std::make_tuple(t.x, t.pos, t.loop));
            return os;
        };
    };

    /**************************************************************************************/

    /// Strand # -> (Complex #, index within complex, loop #)
    vec<StrandData> strand_map;
    /// Complex # -> Ordered strand indices in that complex
    vec<Indices> complex_indices;
    /// Loop join data
    vec<Mat> join_rates;
    /// Complex join data
    Fenwick<Mat> complex_rates, x_rates_sq;
    //bool disable_joins;

    NUPACK_REFLECT(ComplexSet, strand_map, complex_indices, join_rates, complex_rates, x_rates_sq);

    /**************************************************************************************/

    /// Total join rate for state
    template <class RF> real join_rate(RF const &) const;
    /// Check that indices of strand_map and complex_indices are consistent
    bool check() const;

    /**************************************************************************************/

    /// Choose a join move
    template <class RF>
    JoinMove get_join_move(real, RF const &) const;
    /// Choose a join move
    JoinMove get_join_move_nondimensional(real) const;
    /// Register that strands i and j have joined
    void register_join(Index i, Index j);
    /// Register that strands i and j have split
    void register_split(Index i, Index j);

    /**************************************************************************************/

    ComplexSet() : complex_rates(Mat(zero)), x_rates_sq(Mat(zero)) {}
    ComplexSet(std::size_t n) : strand_map(n), join_rates(n, Mat(zero)), complex_rates(Mat(zero)), x_rates_sq(Mat(zero)) {}

    /// Reserve space for "s" number of strands
    //void reserve(Index s) {strand_map.resize(s); join_rates.resize(s, );}
    /// Add complex by giving its strand indices
    template <class... Args> void emplace_back(Args&&...);
    /// Reorder indices so strand "s" is first in its complex (unused)
    void rotate(Index s);
    /// Update join rates for strand "i"
    void update_join_rates(Index i, BaseMat<real> const &);
    /// Update loop index of strand "i" to "o"
    void set_loop_index(Index i, Index o) {strand_map[i].loop = o;}

    auto const & iter() const {return complex_indices;}

    /**************************************************************************************/

    friend std::ostream & operator<<(std::ostream &os, ComplexSet const &c) {dump_os(os, c.complex_indices); return os;}
};


/******************************************************************************************/

template <class RF>
real ComplexSet::join_rate(RF const &rf) const {
    if (complex_indices.size() < 2) return 0.0;
    auto mat = la::eval(times_transpose(complex_rates.total()) - x_rates_sq.total());
    return 0.5 * accu(mat) * (rf.molarity * rf.bimolecular_scaling);
}

/******************************************************************************************/

template <class... Args>
void ComplexSet::emplace_back(Args&&... args) {
    complex_indices.emplace_back(fw<Args>(args)...);
    for (auto i : indices<iseq>(complex_indices.back()))
        strand_map[complex_indices.back()[i]] = {static_cast<iseq>(complex_indices.size() - 1), i, static_cast<iseq>(-1)};
    complex_rates.emplace_back(zero);
    x_rates_sq.emplace_back(zero);
}

/******************************************************************************************/

template <class RF>
JoinMove ComplexSet::get_join_move(double r, RF const &rf) const {
    NUPACK_REQUIRE(r, <, join_rate(rf));
    NUPACK_REQUIRE(r, >=, 0);
    return get_join_move_nondimensional(2 * r / (rf.molarity * rf.bimolecular_scaling));
}

/******************************************************************************************/

/// Call a callback each join location in complex "x" of state "w" of base "b"
// that can pair to an external base "c"
template <class State, class X, class F>
void for_join_sites_in_complex(State &&w, X const &x, Base b, Base c, F &&f) {
    // for each exterior loop in the complex
    fork(w.rate_function, [&](auto const &rf) {
        for_exterior_loops_in_complex(w, x, [&](auto &&o) {
            // for each matching base in the exterior loop
            for_join_locs_in_loop(o, b, c, w.model, rf, [&](auto &&m) {
                // call back with the loop index and base location
                f(o.index(), m);
            });
        });
    });
}

/******************************************************************************************/

/// Call a callback on each ComplexJoin Move possible between complexes x and y in the given State
template <class State, class X, class F>
void for_joins_between(State const &w, X const &x, X const &y, F &&f) {
    // For every base pair (e.g. AT GC CG TA)
    for (auto b : CanonicalBases) for (auto c : CanonicalBases) if (w.model.pairable.can_close(b, c))
        // For every base in the first complex of nucleotide "b" which can pair to "c"
        for_join_sites_in_complex(w, x, b, c, [&](auto const &o1, auto const &m1){
            // For every base in the second complex of nucleotide "c" which can pair to "b"
            for_join_sites_in_complex(w, y, c, b, [&](auto const &o2, auto const &m2){
                // Call back with the join move between these bases
                f(ComplexJoinMove{static_cast<iseq>(o1), static_cast<iseq>(o2), m1, m2});
            });
        });
}

/******************************************************************************************/

/// Call a callback on each ComplexJoin Move possible between complexes in the State
template <class State, class F>
void for_all_joins(State const &w, F &&f) {
    // iterate over all complexes (x, y) where x < y
    for_ordered_pairs(begin_of(w.complexes), end_of(w.complexes),
        [&](auto const &x, auto const &y){
            // callback for each join move between the two complexes
             for_joins_between(w, *x, *y, f);
        }
    );
}

/******************************************************************************************/

}
