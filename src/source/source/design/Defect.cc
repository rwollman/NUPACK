#include <nupack/design/Defect.h>

namespace nupack { namespace newdesign {

Defect::Defect(vec<real> const &defs, real normalization) {
    defect_vec temp;
    izip(defs, [&](auto i, auto d) {if (d > 0) temp.emplace_back(i, d / normalization);});
    contributions = std::move(temp);
}

// Defect evaluate_defect(Design & design) {
//     return Defect();
// }

/**
 * @brief combine defects from repeated underlying nucleotide variables
 *
 * @return Defect contains the collapsed vector (must have same total)
 */
Defect Defect::reduced() const {
    std::map<uint, real> temp;

    for (auto const &c : contributions) {
        auto it = temp.find(c.first);
        if (it == end_of(temp)) it = temp.emplace(c.first, 0).first;
        it->second += c.second;
    }

    return {defect_vec(view(temp))};
}



Defect Defect::weighted(vec<real> const &weights) const {
    NUPACK_REQUIRE(len(weights), ==, len(contributions), "can only apply weights equally");

    Defect temp(*this);
    zip(weights, temp.contributions, [](auto const &w, auto &c) {
        c.second *= w;
    });

    return temp;
}


Defect Defect::scaled(real weight) const {
    Defect temp(*this);
    for (auto &x : temp.contributions) x.second *= weight;
    return temp;
}


/**
 * @brief actually computes the per-nucleotide complex ensemble defect
 * @details
 *
 * @param pp computed pair probabilities to compare against the target structure
 * @param s the target structure
 *
 * @return vector of defects where the index of the vector is the nucleotide
 *     index in the structure/pair probabilities matrix
 */
vec<real> nucleotide_defects(ProbabilityMatrix const &pp, Structure const &s) {
    auto ret = vec<real>(len(s), 0.0);
    for (auto i : indices(s)) ret[i] = 1.0 - pp(i, s[i]);
    return ret;
}


/**
 * @brief actually computes the per-nucleotide complex ensemble defect
 * @details
 *
 * @param pp computed pair probabilities to compare against the target structure
 * @param s the target structure
 *
 * @return vector of defects where the index of the vector is the nucleotide
 *     index in the structure/pair probabilities matrix
 */
vec<real> nucleotide_defects(Tensor<real, 2> const &pp, Structure const &s) {
    auto ret = vec<real>(len(s), 0.0);
    for (auto i : indices(s)) ret[i] = 1.0 - *pp(i, s[i]);
    return ret;
}


/**
 * @brief sample num positions proportional to their contribution to the defect
 * @details [long description]
 *
 * @param num the number of positions to sample. The actual number sampled
 *     will not be higher than the total number.
 * @return the nucleotides sample according to their contributions
 */
vec<uint> Defect::sample_nucleotides(uint num) const {
    auto get_defect = [](auto const &d) {return d.second;};
    auto get_nuc = [](auto const &d) {return d.first;};
    /* edge case: trying to mutate more variables than exist */
    if (num >= len(contributions)) return vmap(contributions, get_nuc);

    vec<uint> sampled_nucs;
    auto distribution = contributions;
    while (num > 0) {
        // recompute sum to avoid floating point error accumulation. can optimize later if necessary.
        auto stop = random_float() * sum(indirect_view(distribution, get_defect));
        auto it = find_cumulative(distribution, stop, get_defect).first;

        sampled_nucs.push_back(get_nuc(*it));

        /* remove sampled element from distribution (sampling w/o replacement) */
        distribution.erase(it);

        --num;
    }
    return sampled_nucs;
}

/**
 * ostensibly used when there is a way of reweighting the defects, for instance by mutation
 * correlation (which is why this was added initially)
 *
 * @return vec<uint> the nucleotides sample according to their reweighted contributions
 */
// vec<uint> Defect::sample_nucleotides(real_mat correlation) const {
//     // correlation = normalise(correlation, 1, 1);

//     uint n = correlation.n_rows;
//     real_col unweighted(n, arma::fill::zeros);
//     for (auto const &c: contributions) unweighted(c.first) += c.second;
//     unweighted = unweighted - (accu(unweighted) / n);
//     unweighted.transform([](auto i) -> real {return (i > 0) ? i : 0;});
//     // real_col reweighted = correlation * unweighted;
//     real_col reweighted = unweighted;

//     real stop = random_float() * accu(reweighted);
//     auto it = find_cumulative(reweighted, stop);
//     uint mut_var =  it.first - reweighted.begin();
//     return {mut_var};
// }


}}
