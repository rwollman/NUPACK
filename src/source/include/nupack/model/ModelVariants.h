#pragma once
#include <cmath>
#include "../types/Sequence.h"
#include "../standard/Variant.h"
#include "StackProgram.h"
#include "ParameterSet.h"

namespace nupack {

/******************************************************************************************/

/// simpler Rigs for use outside dynamic programs
struct RigPF {
    template <class ...Ts>
    static auto plus(Ts const &...ts) {return fold(::nupack::plus, ts...);}

    template <class ...Ts>
    static auto times(Ts const &...ts) {return fold(::nupack::times, ts...);}

    static real zero() {return 0;}
    static real one() {return 1;}

    static auto boltz(real beta, real t) {return boltzmann_factor(beta, t);}
};


struct RigMFE {
    template <class ...Ts>
    static auto plus(Ts const &...ts) {return fold(::nupack::min, ts...);}

    template <class ...Ts>
    static auto times(Ts const &...ts) {return fold(::nupack::plus, ts...);}

    static real zero() {return *inf;}
    static real one() {return 0;}

    static auto boltz(real, real t) {return t;}
};

/******************************************************************************************/

enum class Ensemble : std::uint_fast8_t {nostacking, stacking, min, all, none};

extern std::array<char const *, 5> EnsembleNames;


static std::array<Ensemble, 5> const AllEnsembles
    = {Ensemble::nostacking, Ensemble::stacking, Ensemble::min, Ensemble::all, Ensemble::none};

Ensemble as_ensemble(string_view s);

/******************************************************************************************/

/// Historical dangles=some
struct MinDangles : Empty {
    template <class T>
    constexpr T reduce(T e1, T e2, usize s=3) const {return s == 3 ? std::min(e1, e2) : e1 + e2;}
};

/// Historical dangles=all
struct AllDangles : Empty {
    template <class T>
    constexpr T reduce(T e1, T e2, usize s=3) const {return e1 + e2;}
};

/// New dangles=none
/// Same as historical dangles=none except closing dangle issue
struct NoStacking : Empty {
    template <class T>
    constexpr T reduce(T, T, usize s=3) const {return *zero;}
};

/// New dangles=coax
struct Stacking : Empty {
    template <class T>
    constexpr T reduce(T, T, usize s=3) const {NUPACK_ERROR("should not use."); return *zero;}
};

using EnsembleType = Variant<NoStacking, Stacking, MinDangles, AllDangles>;

inline EnsembleType ensemble_variant(Ensemble e) {
    switch (e) {
        case Ensemble::nostacking: return NoStacking();
        case Ensemble::stacking: return Stacking();
        case Ensemble::none: return NoStacking();
        case Ensemble::min: return MinDangles();
        case Ensemble::all: return AllDangles();
    }
}

/******************************************************************************************/

template <class Dangles, class P>
struct DangleFunction : Dangles {
    P const *p;
    DangleFunction(P const &pset) : p(&pset) {}
    auto energy5(Base i, Base j, Base k) const {return p->data(dangle5, i, j, k);}
    auto energy3(Base i, Base j, Base k) const {return p->data(dangle3, i, j, k);}
};

template <class P>
struct DangleFunction<NoStacking, P> : NoStacking {
    DangleFunction(P const &pset) {}
    value_type_of<P> energy5(Base, Base, Base) const {return *zero;}
    value_type_of<P> energy3(Base, Base, Base) const {return *zero;}
};

template <class D, class P>
DangleFunction<D, P> dangle_function(D, P const &p) {return {p};}

/******************************************************************************************/

template <class M, class V>
value_type_of<M> stacking_energy(NoStacking, M const &, V const &, int const) {return *zero;}

template <class M, class V>
value_type_of<M> stacking_energy(Stacking, M const &model, V const &v, int const nick) {
    auto pf = stacking_sum<real>(v, nick, model);
    return inverse_boltzmann(model.beta, pf);
}

// A list of sequences. nick indicates which one the strand break is before
// For example, the starting loop always has a nick of 0
template <class D, class M, class V, NUPACK_IF(is_same<D, MinDangles, AllDangles>)>
value_type_of<M> stacking_energy(D dangle, M const &model, V const &v, int const nick) {
    using T = value_type_of<M>;
    if (len(v) == 1) return *zero;
    auto const si = begin_of(v), sf = end_of(v) - 1;

    auto get = [&, sn=begin_of(v) + nick](auto s, auto t, auto u) {
        T e5, e3;
        if (sn != t) e5 = model.dG(dangle5, back(*s), front(*t), front(*t, 1));
        if (sn != u) e3 = model.dG(dangle3, back_index(*t, 1), back(*t), front(*u));
        if (sn == t) return e3;
        if (sn == u) return e5;
        return dangle.reduce(e5, e3, len(*t));
    };

    T en = *zero;

    if (len(*si) != 2) en += get(sf, si, si + 1);
    for (auto s = begin_of(v) + 1; s != end_of(v) - 1; ++s)
        if (len(*s) != 2) en += get(s - 1, s, s + 1);
    if (len(*sf) != 2) en += get(sf - 1, sf, si);

    return en;
}

/******************************************************************************************/


/// Dangle at the beginning of an edge, return 0 if s == sn
template <class D, class V, class It>
real safe_dangle5(D const &dangle, V const &v, It s) {
    if (len(*s) <= 2 || front(*s) == Base('_')) return *zero;
    auto const lo = cyclic_prev(v, s);
    return dangle.energy5(back(*lo), front(*s), front(*s, 1));
}

/******************************************************************************************/

/// Dangle at the end_of of an edge, return 0 if s == sn - 1
template <class D, class V, class It>
real safe_dangle3(D const &dangle, V const &v, It s) {
    if (len(*s) <= 2 || back(*s) == Base('_')) return *zero;
    auto const up = cyclic_next(v, s);
    return dangle.energy3(back_index(*s, 1), back(*s), front(*up));
}

/******************************************************************************************/

}
