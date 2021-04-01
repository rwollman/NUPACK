#include <nupack/design/ThermoWrapper.h>
#include <nupack/design/Models.h>
#include <nupack/design/DesignComponents.h>
#include <nupack/model/Model.h>
#include <nupack/thermo/Engine.h>

namespace nupack { namespace newdesign {

EngineObserver NullEngineObserver {};

/**
 * @brief An adapter for thermo::dynamic_program()
 * @details Used to avoid multiple recomputation of dynamic programming
 *     algorithm code.
 *
 * @param seqs The sequence for which to compute the partition function
 * @param models References to two models with the same conditions but with
 *     single- and double-precision datatypes
 *
 * @return The logarithm of the partition function
 */
real partition_function(Local const &env, ::nupack::Complex const &seqs,
        models_type const &models, EngineObserver &engobs) {
    if (!engobs.slowdown)
        return thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, models);

    decltype(thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, models)) ret;
    auto time = time_it(engobs.slowdown, [&](){
        ret = thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, models);
    });
    engobs.log("thermo", "partition function", len(seqs), time, false);
    return ret;
}


/**
 * @brief An adapter for thermo::pair_probability()
 * @details Used to avoid multiple recomputation of dynamic programming
 *     algorithm code.
 *
 * @param seqs The sequence for which to compute the pair probabilities matrix
 * @param models References to two models with the same conditions but with
 *     single- and double-precision datatypes
 *
 * @return A pair with the pair probabilities matrix (Tensor) and the
 *     logarithm of the partition function.
 */
std::pair<thermo::Tensor<real, 2>, real> pair_probability(Local const &env,
        ::nupack::Complex const &seqs, models_type const &models, EngineObserver &engobs) {
    if (!engobs.slowdown)
        return thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, models);

    decltype(thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, models)) ret;
    auto time = time_it(engobs.slowdown, [&]{
        ret = thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, models);
    });
    engobs.log("thermo", "pair probability", len(seqs), time, false);
    return ret;
}

real partition_function(Local const &env, ::nupack::Complex const &seqs, ThermoEnviron &t_env,
        EngineObserver &engobs) {
    return fork(std::get<0>(t_env.models).energy_model.ensemble_type(), [&](auto x) {
        auto &underlying_cache = std::get<DesignCache<decltype(x)>>(t_env.cache);
        decltype(thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache)) ret;
        if (!engobs.slowdown)
            return thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache);

        auto time = time_it(engobs.slowdown, [&]{
            ret = thermo::dynamic_program<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache);
        });
        engobs.log("thermo", "partition function", len(seqs), time, true);
        return ret;
    });
}


std::pair<thermo::Tensor<real, 2>, real> pair_probability(
        Local const &env, ::nupack::Complex const &seqs,
        ThermoEnviron &t_env, EngineObserver &engobs) {
    return fork(std::get<0>(t_env.models).energy_model.ensemble_type(), [&](auto x) {
        auto &underlying_cache = std::get<DesignCache<decltype(x)>>(t_env.cache);

        real pfunc;
        auto obs = [&](auto const &m) {if (m.sequences == seqs.views()) pfunc = m.raw_result;};

        decltype(thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache, obs)) ret;
        try {
            if (!engobs.slowdown) {
                ret = thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache, obs);
            } else {
                auto time = time_it(engobs.slowdown, [&]{
                    ret = thermo::pair_probability<3, 0, 0, 1, 1>(env, seqs, t_env.doubled(), underlying_cache, obs);
                });
                engobs.log("thermo", "pair probability", len(seqs), time, true);
            }
        } catch (Error const &e) {
            auto it = underlying_cache.find(seqs);
            if (it != end_of(underlying_cache)) {
                BEEP(it->first);
                BEEP(fork(it->second, [](auto const &x) {return std::get<1>(x.contents).size();}));
            }

            BEEP(seqs);
            throw e;
        }

        ret.second = pfunc;
        return ret;
    });
}


/**
 * @brief Adapts thermo::pair_probability() for a sequence where fixed_pairs
 *     are forced to to pair by adding a bonus energy.
 * @details Essential function for computing conditional pair probabilities
 *     and partition functions for nodes in the decomposition tree. This
 *     function removes the extra bonus contributions from the partition
 *     function and pair probabilities before returning.
 *
 * @param seqs The sequence for which to compute the clamped pair
 *     probabilities matrix
 * @param models References to two CachedModels with the same conditions but
 *     with single- and double-precision datatypes
 * @param fixed_pairs pairs which are enforced by multiplying their \f$Q^b\f$
 *     element by bonus
 * @param bonus the energy whose boltzmann factor is multiplied into the
 *     \f$Q^b\f$ elements of the enforced pairs to enforce them
 * @return A pair with the pair probabilities matrix (Tensor) and the
 *     logarithm of the partition function. Both have bonuses removed already.
 */
std::pair<thermo::Tensor<real, 2>, real> pair_probability(
        Local const &env,
        ::nupack::Complex const &seqs,
        models_type const &models,
        vec<SplitPoint> const &fixed_pairs,
        real bonus,
        EngineObserver &engobs
        ) {
    auto const &cachedmod = std::get<1>(models);
    auto const &mod = cachedmod.energy_model;
    real exp_bonus = mod.boltz(bonus);
    std::tuple<float, real, overflow<real32>, overflow<real64>> bonuses;
    std::get<0>(bonuses) = exp_bonus;
    std::get<1>(bonuses) = exp_bonus;
    std::get<2>(bonuses) = simd::ifrexp(exp_bonus);
    std::get<3>(bonuses) = simd::ifrexp(exp_bonus);

    auto pairing = [&, n=len(seqs)](auto i, auto j, bool can_pair, auto const & A,
                auto const &block, auto const & s, auto const &model, auto && recursion) {
        auto orig_i = i + s.offset, orig_j = j + s.offset;
        auto distance = std::abs(int(orig_i) - int(orig_j));

        /* modulo and swap necessary so that bonuses are added to both Q^b(i,j)
        and Q^b(j, i+n) */
        i = orig_i % n;
        j = orig_j % n;
        if (j < i) std::swap(i, j);
        // normal behavior; early exit
        bool fixed = contains(fixed_pairs, SplitPoint{i, j});
        bool adjacent = distance == 1;
        bool normal = !fixed && can_pair;
        using V = value_type_of<decltype(block)>;

        return A.sum(
            (normal)             ? A.maybe() & recursion() : A.zero(),
            (fixed && adjacent)  ? A.maybe() & std::get<V>(bonuses) : A.zero(),
            (fixed && !adjacent) ? A.maybe() & A.product(recursion(), std::get<V>(bonuses)) : A.zero()
        );

    };

    real pfunc;
    bool use_B = contains(fixed_pairs, SplitPoint{0, len(seqs) - 1});
    auto obs = [&](auto const &m) {if (m.sequences == seqs.views()) pfunc = m.raw_result;};

    decltype(thermo::bonus_pair_probability<3, 0, 0, 1, 1>(env, seqs, models, False(), obs, pairing, use_B)) ret;

    if (!engobs.slowdown) {
        ret = thermo::bonus_pair_probability<3, 0, 0, 1, 1>(env, seqs, models, False(), obs, pairing, use_B);
    } else {
        auto time = time_it(engobs.slowdown, [&]{
            ret = thermo::bonus_pair_probability<3, 0, 0, 1, 1>(env, seqs, models, False(), obs, pairing, use_B);
        });
        engobs.log("thermo", "bonused pair probability", len(seqs), time, true);
    }

    // remove extra bonuses from fixed base pairs
    auto & pair_probs = ret.first;
    for (auto const &sp : fixed_pairs) {
        uint i, j;
        std::tie(i, j) = sp;
        *pair_probs(i, j) = *pair_probs(j, i) = *pair_probs(i, j) / exp_bonus;
    }

    /* fix the diagonal unpaired probabilities */
    for (auto i : indices(pair_probs)) {
        real acc = 0;
        for (auto j : indices(pair_probs)) acc += (i == j) ? 0 : *pair_probs(i, j);
        *pair_probs(i, i) = 1.0 - acc;
    }

     real terminal_correction = 1.0;
     if (mod.has_terminal_penalty && use_B) {
        terminal_correction *= cachedmod.terminal(at(seqs.catenated, 0), at(seqs.catenated, len(seqs) - 1));
     }

    /* assumes unique pairs */
    ret.second = pfunc - len(fixed_pairs) * std::log(exp_bonus) - std::log(terminal_correction);

    if (std::isnan(ret.second)) NUPACK_ERROR("bonused DPA generated NaN", seqs, pfunc, len(fixed_pairs), std::log(exp_bonus));
    return ret;
}

}}
