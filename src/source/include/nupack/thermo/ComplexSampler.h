/**
 * @brief Sampling structures out of an ensemble of different complexes
 *
 * @file ComplexSampler.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../Forward.h"
#include "../common/Random.h"
#include "../common/Error.h"
#include "../iteration/Patterns.h"
#include "../standard/Vec.h"
#include "../types/Sequence.h"
#include "../types/PairList.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

// V: a vector of sequences that does not need to be sorted
// PFs: a sorted vector of sequence sets to PFs
// scale: a multiplicative penalty factor for complex joins that should usually be < 1
template <class V, class Q>
auto disconnected_weights(V const &v, Q const &logq, real const scale) {
    vec<std::pair<vec<Complex>, real>> ret;
    for_partitions(false, indices(v), [&](auto const &i) {
        real total_logq = (len(v) - len(i)) * log(scale);
        vec<Complex> complexes;
        for (auto const &x : i) {
            complexes.emplace_back(indexed_view(x, v));
            complexes.back().rotate_lowest();
            auto it = binary_search(logq, complexes.back(), first_of);
            NUPACK_REQUIRE(it, !=, end_of(logq), v);
            total_logq += it->second;
        }
        ret.emplace_back(std::move(complexes), total_logq);
    });
    auto max = maximum(ret, second_of).second;
    if (is_finite(exp(max))) max = 0; // biggest one is OK magnitude
    for (auto &i : ret) i.second = exp(i.second - max);
    return ret;
}

/******************************************************************************************/

struct ComplexSampler {

    template <class V, class Q>
    ComplexSampler(V &&seqs, Q &&logq, real const scale) : strands(fw<V>(seqs)),
        complex_logq(fw<Q>(logq)), weights(disconnected_weights(strands, complex_logq, scale)), distribution(item_view(weights)),
        strand_starts(prefixes<vec<iseq>>(true, indirect_view(strands, len))) {}

    StrandList strands; // # strands
    vec<std::pair<Complex, real>> complex_logq; // # complexes
    vec<std::pair<vec<Complex>, real>> weights; // # sets of complexes
    std::discrete_distribution<> distribution;
    vec<iseq> strand_starts; // # strands -> start of sequence

    NUPACK_REFLECT(ComplexSampler, strands, complex_logq, weights, distribution, strand_starts);

    template <class Env, class Ms, class RNG = decltype(StaticRNG) &>
    auto operator()(Env &&env, Ms const &models, usize n, bool shuffle=true, RNG &&gen=StaticRNG)
    {
        auto samples = vmap(complex_logq, [](auto const &p) {
            return std::make_tuple(p.first, vec<PairList>(), usize());
        });

        auto const picks = weighted_samples(distribution, n, gen);

        zip(picks, weights, [&](auto n, auto const &p) {
            for (auto const &x : p.first) third_of(*binary_search(samples, x, first_of)) += n;
        });

        for (auto &i : samples)
            if (third_of(i)) second_of(i) = first_of(sample(env, third_of(i), 1, first_of(i), models));

        vec<PairList> ret(n, PairList(sum(strands, len)));
        auto r = begin_of(ret), s = r;
        zip(picks, weights, [&](auto n, auto const &p) {
            small_vec<bool> mask(len(strands), true);
            for (auto const &x : p.first) {
                auto &sample = *binary_search(samples, x, first_of);
                auto map = reserved<vec<iseq>>(len(x));
                for (auto s : x.views()) {
                    auto i = find_with_mask(strands, s, mask).second - begin_of(mask);
                    mask[i] = false;
                    extend(map, indices(s).shift(strand_starts[i]));
                }
                s = r;
                for (auto rep : range(n)) {
                    izip(second_of(sample).back(), [&](auto i, auto j) {(*s)[map[i]] = map[j];});
                    second_of(sample).pop_back();
                    ++s;
                }
            }
            r = s;
        });
        if (shuffle) random_shuffle(ret, gen);
        return ret;
    }
};

/******************************************************************************************/

}

}
