/**
 * @brief Wrapper for energy model that holds exponentiated Boltzmann factors, etc.
 *
 * @file CachedModel.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Algebras.h"

#include "../types/Sequence.h"
#include "../types/Matrix.h"
#include "../model/ModelVariants.h"
#include "../model/ParameterSet.h"

#include <atomic>

namespace nupack::thermo {

/******************************************************************************************/

template <class T>
struct ParameterCache {
    multi_array<T, 4, 4> terminal, mismatch_b;
    multi_array<T, 4, 4, 4, 4> mismatch;

    // These members scale with length of maximum sequence so are kept mutable
    mutable Tensor<T, 2> alpha;               /// Multiloop factors
    mutable Tensor<T, 2> gamma;               /// Interior loop factors
    mutable Tensor<T, 1> asymmetry;
    T multi1, multi2, multi12, multi22, multi122;
    mutable iseq n=0;                   /// Number of bases that can be incorporated

    void clear() {n = 0; alpha.clear(); gamma.clear(); asymmetry.clear();}

    template <class I> auto int_scale(I i) const {return mantissa_at(gamma(10, span(0, n)), i);}
    template <class I> auto int_asym(I i) const {return mantissa_at(gamma(11, span(0, n)), i-4);}
    template <class I> auto int_asym(int m, I i) const {return asymmetry(i + (n - m));}

    auto bulge(span s) const {return gamma(8, s);}
    auto rbulge(span s) const {return gamma(9, s.reversed(n));}

    template <class I> auto int_size(I i) const {return gamma(12, i);}
    template <class I> auto int_size(I i, iseq j) const {return gamma(j, i);}
    auto int_rsize(span s, iseq j) const {return gamma(j+4, s.reversed(n));}

    template <class I> auto multi3(I i) const {return alpha(0, i);}
    auto multi3r(span s) const {return alpha(1, s.reversed(n));}

    NUPACK_REFLECT(ParameterCache, multi1, multi2, multi12, multi22, multi122, asymmetry, gamma, alpha, terminal, mismatch, mismatch_b, n);
};

NUPACK_DEFINE_TEMPLATE(is_parameter_cache, ParameterCache, class);


/******************************************************************************************/

/**
 * @brief Model containing cached parameter values from an underlying energy Model
 * Only threadsafe if the capacity does not have to be increased
 * @tparam Rig Algebraic rig (usually PF() or MFE())
 * @tparam Model An energy model like nupack::Model()
 */
template <class Rig, class Model>
class CachedModel : public ParameterCache<typename Model::value_type> {
    using T = typename Model::value_type;
public:
    using base_type = ParameterCache<T>;
    using rig_type = Rig;
    using model_type = Model;
    using value_type = T;
    using is_cached_model = True;
    model_type energy_model;
    iseq int_max = std::numeric_limits<std::make_signed_t<iseq>>::max() / 2; // this number should be at least 8

    NUPACK_EXTEND_REFLECT(CachedModel, base_type, energy_model, int_max);

    /**************************************************************************************/

    constexpr auto rig() const {return Rig();}
    constexpr T zero() const {return Rig::zero();}
    constexpr T one() const {return Rig::one();}

    /**************************************************************************************/

    template <class ...Ts>
    bool can_pair(Ts &&...ts) const {return energy_model.pairable(fw<Ts>(ts)...);}
    bool can_close(Base b, Base c) const {return energy_model.pairable.can_close(b, c);}

    /// Boltzmann factor for partition function contribution
    template <bool B=true, NUPACK_IF(B && !Rig::logarithmic::value)>
    T boltz(T e) const {return energy_model.boltz(e);}
    /// Boltzmann factor for MFE is just Identity
    template <bool B=true, NUPACK_IF(B && Rig::logarithmic::value)>
    constexpr T boltz(T e) const {return e;}

    template <class E, bool B=true, NUPACK_IF(B && !Rig::logarithmic::value)>
    auto as_log(E e) const {return log(mantissa(e)) + exponent(e) * LogOf2;}

    template <class E, bool B=true, NUPACK_IF(B && Rig::logarithmic::value)>
    auto as_log(E e) const {return e;}

    template <class E, bool B=true, NUPACK_IF(B && !Rig::logarithmic::value)>
    auto free_energy(E e) const {return -as_log(e) / energy_model.beta;}
    /// Boltzmann factor for MFE is just Identity
    template <class E, bool B=true, NUPACK_IF(B && Rig::logarithmic::value)>
    constexpr auto free_energy(E e) const {return e;}

    T terminal(Base i, Base j) const {return base_type::terminal[i][j];}
    T mismatch(Base i, Base d, Base e, Base j) const {return base_type::mismatch[i][d][e][j];}
    T mismatch(Base d, Base e) const {return base_type::mismatch_b[d][e];}

    template <class Seq>
    T hairpin(Seq const &s) const {return boltz(energy_model.hairpin_energy(s));}
    template <class Seq>
    T interior(Seq const &s, Seq const &t) const {return boltz(energy_model.interior_energy(s, t));}

    /**
     * @brief Get multistranded complex partition function by applying join penalty and rotational symmetry
     * @param t Raw partition function
     * @param v List of sequences
     */
    template <class V> T complex_result(T const t, V const &v) const {
        T join = (len(v) - 1) * energy_model.join_penalty();
        return t + (Rig::logarithmic::value ? join : -energy_model.beta * join);
    }

    /**************************************************************************************/

    // Dangle used for non-coaxial algorithm - assumes the paired bases are Watson Crick complement
    template <class Seq>
    T dangle(iseq i, iseq j, Seq const &s) const {
        if (energy_model.ensemble == Ensemble::stacking) return one();
        if (len(s.nicks()) > 1) return zero();
        if (len(s.nicks()) == 1 && (i == 0 || j == len(s) - 1)) return zero();

        return energy_model.dangle_switch([&](auto const &dangle) {
            T d3 = 0, d5 = 0;
            if (i != 0) d5 = dangle.energy5(complement(s[i-1]), s[i-1], s[i]);
            if (j != len(s) - 1) d3 = dangle.energy3(s[j], s[j+1], complement(s[j+1]));
            if (i == 0) return boltz(d3);
            if (j == len(s) - 1) return boltz(d5);
            return boltz(dangle.reduce(d3, d5, j - i + 3)); // last argument to match definition of size in ModelVariants
        });
    }

    // Full dangle function used for coaxial stacking algorithm
    template <class Seq>
    T dangle(iseq d3, iseq b3, iseq b5, iseq d5, Seq const &s) const {
        if (!can_pair(s[b3], s[b5])) return zero();
        T out;
        if (d5 != b5 && d3 != b3) out = energy_model.terminal_mismatch(s[b3-1], s[b3], s[b5], s[b5+1]);
        else if (d5 != b5) out = energy_model.dG(dangle5, s[b3], s[b5], s[b5+1]);
        else if (d3 != b3) out = energy_model.dG(dangle3, s[b3-1], s[b3], s[b5]);
        else out = 0;
        return boltz(out);
    }

    // i j paired, k l paired, k = j + 1
    T coaxial(Base i, Base j, Base k, Base l) const {
        if (!can_pair(i, j) || !can_pair(k, l)) return zero();
        return boltz(energy_model.coaxial_stack_energy(i, j, k, l));
    }

    CachedModel() = default;

    CachedModel(Rig, Model mod) : energy_model(std::move(mod)) {
        NUPACK_ASSERT(energy_model.valid(), "Empty parameters");
        using C = base_type;

        for (auto i : CanonicalBases) for (auto j : CanonicalBases) {
            C::terminal[i][j] = boltz(energy_model.terminal_penalty(i, j));
            C::mismatch_b[i][j] = boltz(energy_model.interior_mismatch(Base('A'), i, j, Base('A')));
            for (auto d : CanonicalBases) for (auto e : CanonicalBases)
                C::mismatch[i][d][e][j] = boltz(energy_model.interior_mismatch(i, d, e, j));
        }

        C::multi1 = boltz(energy_model.multi_init());
        C::multi2 = boltz(energy_model.multi_pair());
        C::multi12 = boltz(energy_model.multi_pair() + energy_model.multi_init());
        C::multi22 = boltz(energy_model.multi_pair() * 2);
        C::multi122 = boltz(2 * energy_model.multi_pair() + energy_model.multi_init());
    }

    explicit CachedModel(Model mod) : CachedModel(Rig(), std::move(mod)) {}

    /// Change the temperature, which means the cache must be cleared
    void set_beta(T f) {energy_model.beta = f; *this = CachedModel(std::move(energy_model));}

    /// Interconversion between CachedModel() of a different type
    template <class M>
    explicit CachedModel(CachedModel<Rig, M> const &mod) : energy_model(mod.energy_model) {}

    /// Calculate cached elements for calculation of sequence up to length n
    void force_reserve(iseq m) const {
        using C = base_type;
        auto const all = span(0, m);
        NUPACK_ASSERT(energy_model.valid(), "Empty model");
        // Newly inserted check. If parameters are perturbed to be very negative, an infinite Boltzmann factor
        // could be encountered. For instance, a 1000 length multiloop with a2 = -10.0.
        // It's assumed that this can be safely reassigned to 0. This was only
        // encountered during some parameter optimization trials.
        auto Q = [&](auto dG) {
            auto out = boltz(dG);
            return (Rig::logarithmic::value || std::isfinite(out)) ? out : 0;
        };

        C::alpha.resize(2, m);
        C::alpha(0, all).map([&](auto i) {return Q(i * energy_model.multi_base());});
        C::alpha(1, all).assign(reversed(C::alpha(0, all)));

        C::gamma.resize(13, m);
        for (auto min : range(4)) for (auto i : range(min ? 0 : 1, m))
            *C::gamma(4+min, m-i-1) = *C::gamma(min, i) = Q(std::min(i * energy_model.dG(ninio, min-1), energy_model.dG(ninio, ninio.back())) + energy_model.interior_size_energy(i+min+min));

        for (auto s : range(1, min(m, bulge_size.size())))
            *C::gamma(8, s) = Q(energy_model.dG(bulge_size, s - 1));
        for (auto s : lrange(bulge_size.size(), m))
            *C::gamma(8, s) = Q(energy_model.dG(bulge_size, bulge_size.back()) + std::log(s / 30.0) * energy_model.dG(log_loop_penalty));
        C::gamma(9, all).assign(reversed(C::gamma(8, all)));

        C::gamma(10, all).map([&](auto i) {
            return i == 0 ? 0 : Q(energy_model.interior_size_energy(i+2) - energy_model.interior_size_energy(i));
        });
        C::gamma(11, all).map([&](auto i) {
            return Q(energy_model.interior_asymmetry(i, 4) + energy_model.interior_size_energy(i+4));
        });
        C::gamma(12, all).map([&](auto i) {
            return i == 0 ? 0 : Q(energy_model.interior_size_energy(i));
        });

        C::asymmetry.resize(2 * m);
        C::asymmetry(span(0, 2*m)).map([&](auto i) {
            return Q(energy_model.interior_asymmetry(i > m ? 4+i-m : 4+m-i, 4));
        });
        C::n = m;
    };

    iseq capacity() const {return this->n;}
    bool reserve(iseq m) const {return m > capacity() ? (force_reserve(m), true) : false;}
};

// template <class Rig, class Model>
// CachedModel(Model mod) -> CachedModel<Rig, Model>;

NUPACK_DEFINE_TEMPLATE(is_cached_model, CachedModel, class, class);

/******************************************************************************************/

}



