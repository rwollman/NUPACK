#include <nupack/design/Designer.h>
#include <nupack/design/Archive.h>


namespace nupack { namespace newdesign {


void Archive::reevaluate(Local const &env, Designer &designer, uint depth, EnsemblePartition const &part) {
    for (auto &res : results) {
        res = designer.reevaluate_objectives(env, res, depth, part);
        /* save full evaluation while grabbing elements from cache is cheap */
        if (depth == 0 && part.all_active()) res.full_evaluation(designer);
    }
}

/**
 * @brief remove any dominated results.
 * relevant after reevaluating defects at a more accurate estimate.
 *
 * @return true some results removed
 * @return false no results removed
 */
uint Archive::remove_dominated() {
    auto reference = results;
    erase_if(results, [&](auto const &res) {
        return any_of(reference, [&](auto const &res2) {return res2 < res;});
    });
    return len(reference) - len(results);
}


/**
 * @brief remove any dominated results.
 * relevant after reevaluating defects at a more accurate estimate.
 *
 * @return true some results removed
 * @return false no results removed
 */
uint Archive::remove_dominated_by(Result const &res) {
    uint init_len = len(results);
    erase_if(results, [&](auto const &res1) { return res1 > res; });
    return init_len - len(results);
}


/**
 * @brief
 *
 * @param env
 * @param designer
 * @param depth
 * @param part
 * @return true
 * @return false
 */
uint Archive::update_estimates(Local const &env, Designer &designer, uint depth, EnsemblePartition const &part) {
    reevaluate(env, designer, depth, part);
    return remove_dominated();
}


/**
 * @brief resolve current results with incoming result,
 * rejecting result if it is dominated by anything in the results,
 * removing anything dominated by the result and adding automatically in that case,
 * and adding if mutually non-dominating with results and there is space
 *
 * @param res result to try to add to the results
 * @return true res was added
 * @return false res wasn't added
 */
std::pair<uint, uint> Archive::attempt_add(Result const &res) {
    // BEEP(len(results), max_size);

    if (any_of(results, [&](auto const &res1) { return res1 < res || res1 == res; })) {
        // BEEP("new dominated");
        return {0, 0};
    }

    uint new_dominates = remove_dominated_by(res);
    if (new_dominates) {
        results.emplace_back(res);
        // BEEP("new dominates");
        return {1, new_dominates};
    }

    if (!full()) {
        results.emplace_back(res);
        // BEEP("results not full");
        return {1, 0};
    }

    /* condition for diversity promotion */
    auto ds = densities();
    auto min_el = min_element(ds);
    if (min_el == end_of(ds)) return {0, 0};
    if (density(res) > *min_el) {
        uint pos = min_el - begin_of(ds);
        results.erase(begin_of(results) + pos);
        results.emplace_back(res);

        // BEEP("improves diversity");
        return {1, 1};
    }

    // BEEP("nothing else");
    return {0, 0};
}



std::pair<uint, uint> Archive::merge(Archive const &other) {
    auto temp = vmap(other.results, [&](auto const &res) { return attempt_add(res); });
    uint added = sum(key_view(temp));
    uint removed = sum(item_view(temp));
    return {added, removed};
}


/**
 * @brief compute nearest neighbor densities for results members
 *
 * @return vec<real> densities
 */
vec<real> Archive::densities() const {
    return vmap(results, [&](auto const &r) { return density(r); });
}


/**
 * @brief compute density of a particular
 *
 * @param res
 * @return real
 */
real Archive::density(Result const &res) const {
    auto distances = vmap(results, [&](auto const &r) { return distance(res, r); });
    erase_if(distances, [](auto x) { return x <= 0; });
    return *min_element(distances);
}


/**
 * @brief compute distance between defect totals vectors.
 * Average L1 distance normalized by
 *
 * @param res1
 * @param res2
 * @return real
 */
real Archive::distance(Result const &res1, Result const &res2) const {
    NUPACK_REQUIRE(len(res1.defects), ==, len(res2.defects), "incomparable vectors");

    real total_distance = 0;
    zip(res1.totals(), res2.totals(), [&](auto a, auto b) {
        total_distance += std::abs(a - b);
    });
    return total_distance / len(res1.totals());
}


}}
