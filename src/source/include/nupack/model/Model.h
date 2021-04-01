#pragma once

#include "ModelVariants.h"
#include "ParameterSet.h"
#include "StackEnumeration.h"
#include "../algorithms/Utility.h"
#include "../algorithms/Numeric.h"
#include "../standard/Ptr.h"
#include "../standard/Optional.h"
#include "../standard/Map.h"
#include "../types/Complex.h"

namespace nupack {

/******************************************************************************************/

/// Condition descriptors, these are fine as reals since they aren't often used
struct ModelConditions {
    real temperature = DefaultTemperature;
    real na_molarity = 1.0;
    real mg_molarity = 0.0;

    NUPACK_REFLECT(ModelConditions, temperature, na_molarity, mg_molarity);
    using is_member_ordered = True;
};

/******************************************************************************************/

template <class T>
struct Model : MemberOrdered {
    using value_type = T;
    Model() = default;
    explicit Model(Ensemble, ParameterFile const &p={}, ModelConditions const &cs={}, Optional<WobblePairing> gu={});

    /**************************************************************************************/

    NUPACK_REFLECT(Model, parameters, beta, conditions, possible_pairs, has_terminal_penalty, pairable, ensemble);

    std::array<small_vec<Base, 4>, 4> possible_pairs;
    ParameterSet<T> parameters;
    ModelConditions conditions;
    T beta;
    Ensemble ensemble;
    Pairable pairable;
    bool has_terminal_penalty = false;

    bool valid() const {return bool(parameters.data.array);}

    template <class ...Is>
    decltype(auto) dG(Is ...is) const {NUPACK_DASSERT(valid(), "Empty model"); return parameters.data(is...);}

    template <class U>
    Model(Model<U> const &o) {
        static_assert(!is_same<U, T>, "Should use normal copy constructor");
        members_of(*this) = members_of(o);
    }

    auto boltz(T e) const {return boltzmann_factor(beta, e);};

    /**************************************************************************************/

    template <class S> T hairpin_energy(const S &) const;
    template <class S1, class S2> T interior_energy(const S1 &, const S2 &) const;
    template <class V> T linear_multi_energy(V const &) const;
    template <class V> T multi_energy(V const &) const;
    template <class V> T exterior_energy(V const &, int const) const;
    template <class V> T loop_energy(V const &, int) const;
    template <class V> T terminal_penalty_sum(V const &v) const;

    /**************************************************************************************/

    T interior_size_energy(int) const;
    T interior_asymmetry(int, int) const;
    T interior_mismatch(Base, Base, Base, Base) const;
    T terminal_mismatch(Base, Base, Base, Base) const;

    T join_penalty() const {return dG(::nupack::join_penalty);}
    T multi_init() const {return dG(::nupack::multi_init);}
    T multi_base() const {return dG(::nupack::multi_base);}
    T multi_pair() const {return dG(::nupack::multi_pair);}
    auto pairs(Base i) const {return possible_pairs[i];}
    T terminal_penalty(Base i, Base j) const {return dG(::nupack::terminal_penalty, i, j);}
    T coaxial_stack_energy(Base, Base, Base, Base) const;

    /**************************************************************************************/

    template <class F>
    auto dangle_switch(F &&f) const {
        switch (ensemble) {
            case Ensemble::nostacking: return fw<F>(f)(dangle_function(NoStacking(), parameters));
            case Ensemble::stacking: return fw<F>(f)(dangle_function(Stacking(), parameters));
            case Ensemble::none: return fw<F>(f)(dangle_function(NoStacking(), parameters));
            case Ensemble::min: return fw<F>(f)(dangle_function(MinDangles(), parameters));
            case Ensemble::all: return fw<F>(f)(dangle_function(AllDangles(), parameters));
        }
    }

    auto ensemble_type() const {return ensemble_variant(ensemble);}

    // T dangle5(Base i, Base j, Base k) const {return dG(::nupack::dangle5, i, j, k);}
    // T dangle3(Base i, Base j, Base k) const {return dG(::nupack::dangle3, i, j, k);}
};

NUPACK_DEFINE_TEMPLATE(is_model, Model, class);

/******************************************************************************************/

template <class T>
Model<T>::Model(Ensemble e, ParameterFile const &p, ModelConditions const &cs, Optional<WobblePairing> gu)
    : conditions(cs), parameters({p, "dG", dna_salt_correction(cs.temperature,
        cs.na_molarity, cs.mg_molarity), cs.temperature}), beta(1.0 / (Kb * cs.temperature)), ensemble(e) {

    // Use parameter wobble setting if not specified
    if (gu) pairable.wobble_pairing = (*gu == WobblePairing::on);
    else pairable.wobble_pairing = parameters.default_wobble_pairing;

    pairable.wobble_closing = pairable.wobble_pairing && (e == Ensemble::nostacking || e == Ensemble::stacking);

    for (auto i : CanonicalBases) for (auto j : CanonicalBases) if (pairable(i, j)) {
        possible_pairs[i].emplace_back(j);
        has_terminal_penalty |= terminal_penalty(i, j) != 0;
    }
}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::loop_energy(V const &v, int nick) const {
    if (nick != -1) return exterior_energy(v, nick);
    else if (len(v) == 1) return hairpin_energy(v[0]);
    else if (len(v) == 2) return interior_energy(v[0], v[1]);
    else return multi_energy(v);
}

/******************************************************************************************/

template <class T>
T Model<T>::interior_size_energy(int s) const {
    NUPACK_REQUIRE(s, >, 0);
    if (s <= 30) return dG(interior_size, s - 1); // Interior with >2, <30 total
    return static_cast<T>(dG(interior_size, interior_size.back()) + log((s) / 30.0) * dG(log_loop_penalty)); // Big interior loop
}

/******************************************************************************************/

template <class T>
T Model<T>::interior_asymmetry(int n1, int n2) const {
    auto ninio_number = std::min((decltype(n2)) 4, std::min(n2, n1)) - 1;
    auto asymmetry = std::abs(n1 - n2);
    return std::min<T>(asymmetry * dG(ninio, ninio_number), dG(ninio, ninio.back()));
}

/******************************************************************************************/

template <class T> /// Interior mismatch energy for (b1, b2, b3, b4) where b2 and b3 are paired, b1 left of b2, b4 right of b3
T Model<T>::interior_mismatch(Base b1, Base b2, Base b3, Base b4) const {
    return dG(::nupack::interior_mismatch, b1, b2, b3, b4);
}

/******************************************************************************************/

template <class T> /// Terminal mismatch energy for (b1, b2, b3, b4) where b2 and b3 are paired, b1 left of b2, b4 right of b3
T Model<T>::terminal_mismatch(Base b1, Base b2, Base b3, Base b4) const {
    return dG(::nupack::terminal_mismatch, b1, b2, b3, b4);
}

/******************************************************************************************/

/// Coaxial stack energy for (b1, b2, b3, b4) where b1 and b2 are paired, b3 and b4 are paired,
template <class T> /// and b2 and b3 are on the same strand
T Model<T>::coaxial_stack_energy(Base b1, Base b2, Base b3, Base b4) const {
    return dG(coaxial_stack, b2, b3, b4, b1);
}

/******************************************************************************************/

// Input 5' ---seq1--> 3'
//       3' <--seq2--- 5'
template <class T> template <class S1, class S2>
T Model<T>::interior_energy(S1 const &seq1, S2 const &seq2) const {
    T en = T();
    int const n1 = len(seq1) - 2, n2 = len(seq2) - 2;
    if (n1 == 0 && n2 == 0) { // Stack loop
        return en + dG(stack, seq1[0], seq1[1], seq2[0], seq2[1]);
    } else if (n1 == 0 || n2 == 0) { // Bulge loop
        auto sz = std::max(len(seq1), len(seq2)) - 2;

        if (sz <= 30) en += dG(bulge_size, sz - 1);
        else en += dG(bulge_size, bulge_size.back()) + log(sz / 30.0) * dG(log_loop_penalty);

        if (sz == 1) {
            // add stacking term for single-base bulges. No terminal penalty here
            return en + dG(stack, front(seq1), back(seq1), front(seq2), back(seq2)) - parameters.info.loop_bias;
        } else {
            // Terminal penalty applies otherwise
            if (has_terminal_penalty) en += terminal_penalty(front(seq1), back(seq2))
                                          + terminal_penalty(front(seq2), back(seq1));
            return en;
        }
    } else if (n1 == 1 && n2 == 1) { // Interior 1x1
        return en + dG(interior_1_1, seq1[0], seq1[1], seq1[2], seq2[0], seq2[1], seq2[2]);
    } else if (n1 == 1 && n2 == 2) { // Interior 1x2
        return en + dG(interior_1_2, seq1[0], seq1[1], seq1[2], seq2[0], seq2[1], seq2[2], seq2[3]);
    } else if (n1 == 2 && n2 == 1) { // Interior 2x1
        return en + dG(interior_1_2, seq2[0], seq2[1], seq2[2], seq1[0], seq1[1], seq1[2], seq1[3]);
    } else if (n1 == 2 && n2 == 2) { // Interior 2x2
        return en + dG(interior_2_2, seq1[0], seq1[1], seq1[2], seq1[3], seq2[0], seq2[1], seq2[2], seq2[3]);
    } else { // Big interior loop
        en += interior_size_energy(n1 + n2);
    }
    // interior loops n1 > 4 x n2 > 4 are size + asymmetry + mismatch
    en += interior_asymmetry(n1, n2); // after n1 > 4 && n2 > 4 this just depends on |n1-n2|

    if ((n1 == 1 && n2 > 2) || (n2 == 1 && n1 > 2)) {
        en += interior_mismatch(Base('A'), back(seq2), front(seq1), Base('A'));
        en += interior_mismatch(Base('A'), back(seq1), front(seq2), Base('A'));
    } else {
        en += interior_mismatch(back_index(seq2, 1), back(seq2), front(seq1), front(seq1, 1));
        en += interior_mismatch(back_index(seq1, 1), back(seq1), front(seq2), front(seq2, 1));
    }
    return en;
}

/******************************************************************************************/

template <class T> template <class S>
T Model<T>::hairpin_energy(S const &seq) const {
    T en = T();
    if (len(seq) <= 32) en += dG(hairpin_size, len(seq) - 3);
    else en += dG(hairpin_size, hairpin_size.back()) + log((len(seq) - 2) / 30.0) * dG(log_loop_penalty);

    if (len(seq) == 5) { // triloop
        // terminal penalty is unintuitive but correct
        if (has_terminal_penalty) en += terminal_penalty(back(seq), front(seq));
        return en + dG(hairpin_tri, seq[0], seq[1], seq[2], seq[3], seq[4]);
    } else if (len(seq) == 6) { // tetraloop
        en += dG(hairpin_tetra, seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]);
    }
    return en + dG(hairpin_mismatch, back_index(seq, 1), back(seq), front(seq), front(seq, 1));
}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::linear_multi_energy(V const &v) const {
    auto const n_unpaired = sum(v, len) - 2 * len(v);
    return len(v) * multi_pair() + multi_init() + multi_base() * n_unpaired;
}

/******************************************************************************************/

// Total of terminal penalty and impossible closing base pair energies for multi and exterior loop
template <class T> template <class V>
T Model<T>::terminal_penalty_sum(V const &v) const {
    if (!pairable.wobble_closing) {
        if (front(front(v)) != Base('_') && front(front(v)) + back(back(v)) == 5) return *inf;
        for (auto it : iterators(v).offset(0, -1))
            if (front(it[1]) != Base('_') && front(it[1]) + back(it[0]) == 5) return *inf;
    }

    T t{};
    if (has_terminal_penalty) {
        for (auto it = begin_of(v); it != end_of(v) - 1; ++it)
            if (front(it[1]) != Base('_')) t += terminal_penalty(front(it[1]), back(it[0]));
        if (front(front(v)) != Base('_')) t += terminal_penalty(front(front(v)), back(back(v)));
    }
    return t;
}

/******************************************************************************************/

template <class T> template <class V>
T Model<T>::multi_energy(V const &v) const {
    return terminal_penalty_sum(v)
         + linear_multi_energy(v)
         + fork(ensemble_type(), [&](auto d) -> T {return stacking_energy(d, *this, v, -1);});
}

/******************************************************************************************/

// A list of sequences. nick indicates which one the strand break is before
// For example, the starting loop always has a nick of 0
template <class T> template <class V>
T Model<T>::exterior_energy(V const &v, int const nick) const {
    return terminal_penalty_sum(v)
         + fork(ensemble_type(), [&](auto d) -> T {return stacking_energy(d, *this, v, nick);});
}

/******************************************************************************************/

int find_loop_structure_nick(Complex const &, PairList const &);

/******************************************************************************************/

template <class Model>
std::map<string, double> loop_stacking_energies(Model const &m, Complex const &c, int nick=Ether) {
    std::map<string, double> out;
    auto const v = complex_to_loop(c, nick);
    if (nick != Ether || c.n_strands() > 2) {
        auto const t = m.terminal_penalty_sum(v) + (nick != Ether ? 0 : m.linear_multi_energy(v));
        enumerate_stacking_state_energies(v, nick, m, [&](auto const &p, auto e) {
            out.emplace(loop_stack_sequence_string(p), t + e);
        });
    } else {
        out.emplace(string(c.n_strands(), 'n'), m.loop_energy(v, nick));
    }
    return out;
}

/******************************************************************************************/

}
