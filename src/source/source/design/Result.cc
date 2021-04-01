#include <nupack/design/Result.h>
#include <nupack/design/Designer.h>


namespace nupack { namespace newdesign {

Defect Result::defect(uint i) const {
    return at(defects, i);
}

Defect Result::weighted_defect(uint i) const {
    if (len(weights) != len(defects)) NUPACK_BUG("weights must be same length as defects", len(defects), len(weights), defects, weights);
    return at(defects, i).scaled(at(weights, i));
}

vec<Defect> Result::weighted_defects() const {
    return vmap(range(len(defects)), [&](auto i) { return weighted_defect(i); });
}

real Result::total() const {
    return sum(totals());
}


real Result::weighted_total() const {
    return sum(weighted_totals());
}


vec<real> Result::totals() const {
    return vmap(defects, [](auto const &d) {return d.total();});
}


vec<real> Result::weighted_totals() const {
    if (len(weights) != len(defects)) NUPACK_BUG("weights must be same length as defects", len(defects), len(weights), defects, weights);
    vec<real> ret;
    zip(defects, weights, [&](auto const &d, auto w) {
        ret.emplace_back(d.total() * w);
    });
    return ret;
}


SingleResult const & Result::full_evaluation(Designer const &designer) const {
    /* evaluate if not done yet */
    if (len(evaluated.domains) == 0)
        evaluated = SingleResult(designer, *this);
    return evaluated;
}


bool Result::operator<(Result const &other) const {
    if (*this == other) return false;
    bool less_than {true};
    zip(this->totals(), other.totals(), [&](auto a, auto b) {
        if (a > b) less_than = false;
    });
    return less_than;
}


bool Result::operator>(Result const &other) const {
    return other < *this;
}


bool Result::operator==(Result const &other) const {
    return this->totals() == other.totals();
}


bool Result::operator!=(Result const &other) const {
    return !(*this == other);
}


/************************************************************************************************************/
/* sampling */

/**
 * @brief sample from the defect contributions of whatever the first defect is
 *
 * @param res a multiobjective set of defects
 * @param num number of nucleotide indices to mutate
 * @return vec<uint> the nucleotides to mutate; len(ret) = num
 */
vec<uint> first_defect_sample(Result const &res, uint num) {
    return res.defect().sample_nucleotides(num);
}


/**
 * @brief first sample a defect (i.e. an objective) to minimize based on the relative defect total;
 * then sample nucleotides according to the contributions to that defect
 *
 * @param res a multiobjective set of defects
 * @param num number of nucleotide indices to mutate
 * @return vec<uint> the nucleotides to mutate; len(ret) = num
 */
vec<uint> stochastic_hierarchical_sample(Result const &res, uint num) {
    auto stop = random_float() * sum(res.weighted_totals());
    auto it = find_cumulative(res.defects, stop, [] (auto const &d) {return d.total();}).first;

    return it->sample_nucleotides(num);
}


/**
 * @brief sum all the individual defects together and normalized by number of
 * defects (currently no weights)
 *
 * @param res a multiobjective set of defects
 * @param num number of nucleotide indices to mutate
 * @return vec<uint> the nucleotides to mutate; len(ret) = num
 */
vec<uint> scalarized_sample(Result const &res, uint num) {
    real n = len(res.defects);
    vec<real> combined_defects(len(res.sequence), 0.0);
    for (auto const &def : res.weighted_defects()) {
        for (auto const &p : def.contributions)
            at(combined_defects, p.first) += p.second;
    }

    Defect defect(combined_defects, n);

    return defect.sample_nucleotides(num);
}


/**
 * @brief sample nucleotides without considering defects at all (discrete
 * uniform distribution over underlying variables w/o replacement)
 *
 * @param res a multiobjective set of defects
 * @param num number of nucleotide indices to mutate
 * @return vec<uint> the nucleotides to mutate; len(ret) = num
 */
vec<uint> uniform_sample(Result const &res, uint num) {
    vec<uint> choices {indices(res.sequence)};
    vec<uint> ret;

    while (num > 0) {
        auto it = random_choice(choices);
        ret.emplace_back(*it);
        choices.erase(it);

        num--;
        if (len(choices) == 0) num = 0;
    }

    return ret;
}


}}
