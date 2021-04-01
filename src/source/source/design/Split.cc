#include <nupack/design/Split.h>

namespace nupack { namespace newdesign {



/**
 * @brief threshold out small pair probabilities and return sparse matrix of the
 * rest
 *
 * @param in
 * @return ProbabilityMatrix
 */
ProbabilityMatrix sparsify(Tensor<real, 2> const &in, real f_sparse) {
    vec<arma::uword> Is;
    vec<arma::uword> Js;
    vec<real> vec_values;

    for (auto i : range(len(in))) for (auto j : range(len(in))) {
        auto value = *in(i, j);
        if (value >= f_sparse) {
            Is.emplace_back(i);
            Js.emplace_back(j);
            vec_values.emplace_back(value);
        }
    }

    real_col values(vec_values);

    uint N = len(Is);
    auto IJ = catted<vec<arma::uword>>(Is, Js);
    arma::umat locations(IJ);
    locations.reshape(N, 2);

    return ProbabilityMatrix(locations.t(), values, len(in), len(in), true);
}


/**
 * @brief check if the children implied by a split point are both larger than
 *     the minimum size.
 * @details [long description]
 *
 * @param sp the split point used to divide the complex
 * @param n the number of nucleotides in the complex
 * @param min_size the minimum number of nucleotides that must be in a child
 *     resulting from the split
 *
 * @return whether both children are larger than min_size
 */
bool is_large_enough(SplitPoint sp, int n, int min_size) {
    int i = sp.first, j = sp.second, left = (i + 1) + (n - j), right = j - i + 1;
    return (left >= min_size && right >= min_size);
}


/**
 * @brief check if a nucleotide is far enough away from the ends of its
 *     strand.
 * @details Ensure that there are at least min_helix nucleotides 5' and 3' of
 *     the given index on the same strand.
 *
 * @param i a nucleotide position
 * @param bounds the index of the first nucleotide of each strand, including
 *     one past-the-end index (same as a nicks array with 0 prepended)
 * @param min_helix the number of nucleotides to 5' and 3' of i required to be
 *     on the same strand.
 *
 * @return whether i - min_helix and i + min_helix are on the same strand.
 */
bool is_padded(int i, Nicks const &bounds, int min_helix) {
    return any_of(range(0, len(bounds) - 1), [&](auto r) {
        return i - min_helix >= int(bounds[r]) && i + min_helix < int(bounds[r+1]);
    });
}


/**
 * @brief produce a vector of all positions that are far enough from their
 *     strand ends.
 * @details used as a first pass to filter out nucleotide positions that are
 *     too close to the 3' or 5' end of the strand they are on.
 *
 * @param bounds the index of the first nucleotide of each strand, including
 *     one past-the-end index (same as a nicks array with 0 prepended)
 * @param min_helix the number of nucleotides to 5' and 3' of i required to be
 *     on the same strand.
 *
 * @return all indices that could be involved in a base pair that individually
 *     meet the padding requirements
 */
small_vec<uint> padded(small_vec<uint> bounds, uint min_helix) {
    small_vec<uint> ret;
    for (auto i : range(back(bounds))) if (is_padded(i, bounds, min_helix)) ret.push_back(i);
    return ret;
}


/**
 * @brief checks if a proposed split point meets child size and helix
 *     requirements to be considered valid.
 * @details Ensures that relative to a target structure, a split point is
 *     surrounded by an adequate number of pairs, both nucleotides are far
 *     enough from the ends of the strands they are on, and
 *
 * @param sp the split point to check for validity
 * @param s the target structure to check the split point against
 * @param min_size the minimum number of nucleotides that must be in both
 *     children implied by sp for sp to be valid
 * @param min_helix the minimum number of stacked base pairs on each side of
 *     the split point base pair that must be present in the target structure
 *
 * @return whether the helix and child size constraints are met for the split
 *     point sp in structure s
 */
bool is_valid(SplitPoint sp, Structure const &s, uint min_size, uint min_helix) {
    int i = sp.first, j = sp.second, n = len(s);
    Nicks bounds {0};
    cat(bounds, s.nicks);

    if (is_padded(i, bounds, min_helix)
            && is_padded(j, bounds, min_helix)
            && is_large_enough(sp, n, min_size)) {
        if (any_of(range(-int(min_helix), int(min_helix + 1)), [&](auto r) {return s[i - r] != j + r;}))
            return false;
    } else
        return false;

    return true;
}


/**
 * @brief returns the collection of base pairs in a target structure that are valid
 *     split points.
 * @details used in generating a list of all split points that meet the child
 *     size and padding requirements.
 *
 * @param s the target structure
 * @param min_size the minimum number of nucleotides that must be in both
 *     children implied by a valid split point
 * @param min_helix the minimum number of stacked base pairs on each side of
 *     a split point base pair that must be present in the target structure
 *
 * @return a list of all split points that are both base pairs in the
 *     structure s and meet the padding and resultant child size requirements
 */
vec<SplitPoint> valid_split_points(Structure const &s, uint min_size, uint min_helix) {
    vec<SplitPoint> splits;
    s.for_each_pair([&] (int i, int j) {
        if (is_valid({i, j}, s, min_size, min_helix)) splits.emplace_back(i, j);
    });
    // BEEP(s.dp(), splits);
    return splits;
}


/**
 * @brief the sum of the cost proxy for evaluating the two children implied by
 *     the split point.
 * @details The asymptotic cost of running the dynamic programs for a sequence
 *     of length \f$n\f$ is \f$O(n^3)\f$. This computes the cost of evaluating
 *     the two children resulting from the split point. The cost is \f$n_l^3 +
 *     n_r^3\f$ for children with lengths \f$n_l\f$ and \f$n_r\f$ (these
 *     lengths are computed from the split point and total length).
 *
 * @param p the split point
 * @param n the length of the complex being split
 *
 * @return the sum of the costs of evaluating
 */
real children_cost(SplitPoint p, uint n) {
    int i = p.first, j = p.second, left = (i + 1) + (n - j), right = j - i + 1;
    return std::pow(left, 3) + std::pow(right, 3);
}


/**
 * @brief sorts the split points from lowest to highest cost.
 * @details used in structure-guided decomposition to produce the split point
 *     whose implied children are cheapest to evaluate.
 *
 * @param splits the split points to sort by child cost
 * @param n the size of the complex associated with the split points
 *
 * @return the same split points as in splits, but sorted from least to most
 *     expensive child evaluation cost
 */
vec<SplitPoint> ascending_cost_splits(vec<SplitPoint> splits, uint n) {
    auto compare_cost = [&](auto const &l, auto const &r) {return children_cost(l, n) < children_cost(r, n);};
    sort(splits, compare_cost);
    // BEEP(n, splits);
    return splits;
}



/**
 * @brief finds all splits consistent with a given structure (if valid) and
 *     meeting the minimum probability threshold based on the
 *     pair_probabilities
 * @details used primarily as a subsidiary function to minimal_splits(), this
 *     function produces the lists of splits points that are used in the
 *     branch and bound procedure in that function. As minimal_splits() is
 *     flexible in that it can produce split points for pure probability-
 *     guided decomposition or structure- and probability-guided
 *     decomposition, this function produces both structure-related and
 *     probability-related lists of split points. The split points coming from
 *     the structure must result in children that are large enough and must be
 *     embedded in a large enough helical region. The split points coming from
 *     the probability matrix must also meet the child size requirement as
 *     well as having the probability of the minimum padding base pair exceed
 *     a threshold (so that pairs flanked by nucleotides that can't possibly
 *     pair are not considered).
 *
 * @param probs the pair probabilities matrix used in finding valid split
 *     points for probability-guided decomposition
 * @param min_size the minimum number of nucleotides that must be in both
 *     children implied by a valid split point
 * @param min_helix the minimum number of stacked base pairs on each side of
 *     a split point base pair that must be present in the target structure
 * @param s the structure to generate structure-guided split points from. If
 *     not a valid structure, it must at least contain a valid nicks array as
 *     this is used in the bounds checking.
 * @return a tuple containing the valid structure-based split points first and
 *     the valid probability-based split points second
 */
auto possible_splits(ProbabilityMatrix const &probs, uint min_size, uint min_helix, Structure const &s)
        -> std::tuple<vec<ProbabilitySplit>, vec<ProbabilitySplit>> {
    /* TODO: can move this into design parameters, possibly, but fine here for now.
        Could just make this 0 to eliminate truly impossible pairs only?
    */
    constexpr real threshold = 0.001;

    auto n = probs.n_rows;
    auto min_prob = [&](int i, int j) {
        return minimum(indirect_view(range(-int(min_helix), int(min_helix + 1)), [&](auto r){return probs(i-r, j+r);}));
        // return *probs(i, j);
    };

    /* Prepare split points from target structure */
    vec<ProbabilitySplit> structure_splits;
    if (s.valid()) {
        structure_splits = vmap(valid_split_points(s, min_size, min_helix), [&](auto const &t) {
            return ProbabilitySplit(t.first, t.second, min_prob(t.first, t.second), children_cost(t, n));
        });
    }

    Nicks bounds {0};
    cat(bounds, s.nicks);
    auto valid_nucs = padded(bounds, min_helix);

    vec<ProbabilitySplit> probability_splits;
    for (auto i : valid_nucs) {
        for (auto j : valid_nucs) {
            bool in_structure = s.valid() && s[i] == j;
            if (!in_structure && i < j && is_large_enough({i, j}, back(bounds), min_size)) {
                probability_splits.emplace_back(i, j, min_prob(i, j), children_cost({i, j}, n));
            }
        }
    }


    auto compare_prob = [](auto const &l, auto const &r) {return l.prob > r.prob;};
    sort(structure_splits, compare_prob);
    sort(probability_splits, compare_prob);
    erase_if(probability_splits, [&] (auto const &p) {return p.prob < threshold;});
    return std::make_tuple(std::move(structure_splits), std::move(probability_splits));
}


/**
 * @brief returns the minimal cost set of exclusive split points whose
 *     collective probability exceeds f_split.
 * @details Applies a branch and bound procedure to produce a set of split
 *     points which capture at least f_split of the total structural ensemble
 *     and do so more cheaply than computing the parent. This function handles
 *     the minimal split case for both pure probability-guided decomposition
 *     and structure- and probability-guided decomposition. It is possible
 *     that there is no set of split points that are mutually exclusive and
 *     meet the f_split requirement, in this case an empty set is returned.
 *     Furthermore, it is possible that the only sets which meet the f_split
 *     requirement have costs exceeding the parent cost, e.g. any set of 4 or
 *     more split points. In this case, an empty set is also returned.
 *
 * @param probs the pair probabilities matrix used in finding valid split
 *     points for probability-guided decomposition
 * @param f_split the minimum probability that must be captured by the set of
 *     split points returned
 * @param min_size the minimum number of nucleotides that must be in both
 *     children implied by a valid split point
 * @param min_helix the minimum number of stacked base pairs on each side of
 *     a split point base pair that must be present in the target structure
 * @param s the structure to generate structure-guided split points from. If
 *     not a valid structure, it must at least contain a valid nicks array as
 *     this is used in the bounds checking.
 * @return a set of split points that meets the child size, minimum
 *     probability, padding, and cost requirements. If this is empty, the
 *     complex cannot be decomposed any further for this set of parameters.
 */
vec<SplitPoint> minimal_splits(ProbabilityMatrix const &probs, real f_split, uint min_size, uint min_helix, Structure const &s) {
    /* Starting splits */
    vec<ProbabilitySplit> structure_splits;
    vec<ProbabilitySplit> probability_splits;
    std::tie(structure_splits, probability_splits) = possible_splits(probs, min_size, min_helix, s);

    // auto cost_sort = [](auto const &a, auto const &b) {return a.cost < b.cost; };

    // BEEP(probs.n_rows);
    // BEEP(structure_splits);
    // BEEP(sorted(structure_splits, cost_sort));
    // BEEP(probability_splits);
    // BEEP(sorted(probability_splits, cost_sort));



    /* Initialize with cost of parent; should NOT decompose at all if the lowest
        cost decomposition meeting the probability constraint is more
        expensive than the parent computation.
    */
    real best_cost = pow(probs.n_rows, 3);
    // real best_cost = std::numeric_limits<real>::max(); // old behavior, sans limit on number of split points
    vec<ProbabilitySplit> best_splits;

    vec<ProbabilitySplit> cur_splits;
    vec<uint> included_positions;
    real cur_prob {0};
    real cur_cost {0};

    /* add spl to the current split set (assuming it has been previously checked to not conflict) */
    auto push = [&](auto const &spl, auto pos) {
        cur_splits.emplace_back(spl);
        included_positions.emplace_back(pos);
        cur_prob += spl.prob;
        cur_cost += spl.cost;
    };

    auto pop = [&] {
        auto ret = back(included_positions);
        cur_splits.pop_back();
        included_positions.pop_back();
        /* recalculate instead of repeated addition and subtraction accumulating error */
        cur_prob = sum(cur_splits, [](auto const &c) {return c.prob;});
        cur_cost = sum(cur_splits, [](auto const &c) {return c.cost;});
        return ret;
    };

    uint struc_pos = 0;
    do {
        uint cur_pos = 0;

        /* Only interact with structure splits if there is a valid structure. However,
            because the first split point must come from the structure splits,
            if there are no valid split points in the structure, this node
            cannot be decomposed through structure- and probability-guided
            decomposition. */
        if (s.valid()) {
            while (!cur_splits.empty()) pop();

            /* Add next structure split if there are more remaining */
            if (struc_pos < len(structure_splits)) {
                push(structure_splits[struc_pos], -1);
                ++struc_pos;
            } else {
                /* If there are no more structure splits, we can end immediately as the entire
                    search tree has been checked. */
                break;
            }

            if (cur_prob >= f_split) {
                if (cur_cost < best_cost) {
                    best_splits = cur_splits;
                    best_cost = cur_cost;
                }
                // cur_pos = pop();
            }
        }

        bool add_probability = (cur_prob < f_split);

        /* sentinel to make sure that collective probability wasn't satisfied with
            single split point from structure. */
        if (add_probability && len(probability_splits) > 0) {
            do {
                // try to add current current position to existing splits
                auto const &spl = probability_splits.at(cur_pos);
                if (all_of(cur_splits, [&](auto const &c) {return crosses(c, spl);})) push(spl, cur_pos);

                /* bound and trim if too expensive */
                if (cur_cost > best_cost) cur_pos = pop();
                /* stop exploring  */
                if (cur_prob >= f_split) {
                    if (cur_cost < best_cost) {
                        best_splits = cur_splits;
                        best_cost = cur_cost;
                    }
                    cur_pos = pop();
                }

                ++cur_pos;
                while (cur_pos == len(probability_splits) && len(cur_splits) > int(s.valid())) cur_pos = pop() + 1;
            } while (cur_pos < len(probability_splits) || len(cur_splits) > int(s.valid()));
        }
    } while (len(cur_splits) > 0);

    // BEEP(sum(best_splits, [](auto const &c) {return c.prob;}));

    // BEEP(best_splits);
    return vmap<vec<SplitPoint>>(best_splits, [](auto const& b) {return SplitPoint{b.first, b.second};});
}


/**
 * @brief create the two child structures implied by the parent structure and the split point
 * @details The pair of resulting child structures are split at the specified
 *     SplitPoint such that both include the nucleotides of the SplitPoint.
 *     Base pairs in the parent structure are sorted into the left, the right
 *     child, both children (the SplitPoint only), or neither (if the pair
 *     crosses the SplitPoint). This function also works for degenerate
 *     structures containing only a list of nicks, and produces two child
 *     degenerate structures with the appropriate nicks.
 *
 * @param sp the SplitPoint at which to divide the structure
 * @param s the parent structure
 *
 * @return a pair of the two child structures. The first is the "left" child,
 *     the second is the "right" child.
 */
std::pair<Structure, Structure> split(SplitPoint const &sp, Structure const &s) {
    Structure left, right;
    uint i, j; std::tie(i, j) = sp;
    auto on_left = [&](auto x) {return x <= i || x >= j;};
    auto on_right = [&](auto x) {return x >= i && x <= j;};
    auto to_left = [&](auto x) {return x <= i ? x : x - j + 1 + i;};
    auto to_right = [&](auto x) {return x - i;};
    /* move pairs into child structures */
    if (s.valid()) {
        auto n = len(s);

        /* create correctly sized and blanked underlying data arrays for the children */
        using data_type = decltype(left.values);
        left.values = linspace<data_type>(i + 1 + n - j);
        right.values = linspace<data_type>(j - i + 1);

        s.for_each_pair([&](auto d, auto e) {
            /* skip pairs made incompatible by the split point */
            if (crosses({d, e}, sp)) return;
            if (on_left(d) && on_left(e)) left.toggle_pair(to_left(d), to_left(e));
            if (on_right(d) && on_right(e)) right.toggle_pair(to_right(d), to_right(e));
        });
    }
    /* create new nicks arrays */
    for_each(s.nicks, [&](auto n) {
        if (on_left(n)) {
            left.nicks.emplace_back(to_left(n));
        }
        if (on_right(n)) {
            right.nicks.emplace_back(to_right(n));
        }
    });
    if (!contains(left.nicks, i+1)) left.nicks.emplace_back(i+1);
    sort(left.nicks);
    if (!contains(right.nicks, j-i+1)) right.nicks.emplace_back(j-i+1);
    sort(right.nicks);

    return {std::move(left), std::move(right)};
}


/**
 * @brief divide the sequence for the parent structure into the two child
 *     sequences implied by the split point
 * @details Both children include the split point nucleotides.
 *
 * @param sp the SplitPoint at which to divide the complex sequence
 * @param seq the parent sequence
 *
 * @return a pair of the two child sequences. The first is the "left" child,
 *     the second is the "right" child.
 */
std::pair<vec<StrandView>, vec<StrandView>> split(SplitPoint const &sp, vec<StrandView> const &seq) {
    uint i, j; std::tie(i, j) = sp;
    uint begin = 0, end = -1, total = 0;
    vec<StrandView> left, right;
    auto on_strand = [&](auto x) {return x >= begin && x <= end;};

    for (auto strand : seq) {
        begin = end + 1; /* begin is inclusive */
        total += len(strand); /* end is exclusive */
        end = total - 1;

        auto si = [&](auto x) {return x - begin;}; // index in current strand

        /* strand entirely in left */
        if (begin > j || end < i) left.emplace_back(strand);
        /* strand entirely in right */
        else if (begin > i && end < j) right.emplace_back(strand);
        /* split is within one strand */
        else if (on_strand(i) && on_strand(j)) {
            left.emplace_back(strand.slice(si(begin), si(i)));
            right.emplace_back(strand.slice(si(i), si(j)));
            left.emplace_back(strand.slice(si(j), si(end)));
        }
        /* split is on separate strands */
        else if (on_strand(i)) {
            left.emplace_back(strand.slice(si(begin), si(i)));
            right.emplace_back(strand.slice(si(i), si(end)));
        }
        else if (on_strand(j)) {
            right.emplace_back(strand.slice(si(begin), si(j)));
            left.emplace_back(strand.slice(si(j), si(end)));
        }
    }


    return {std::move(left), std::move(right)};
}


/**
 * @brief divides a set of enforced pairs in a parent nodes indexing scheme into
 * the appropriate children using their indexing scheme (including the split
 * point).
 *
 * @param sp the split point dividing the parent
 * @param pairs the set of enforced pairs in the parent
 * @return a pair of vectors of enforced pairs. First is for the "left" child,
 * second is for the "right" child.
 */
std::pair<vec<SplitPoint>, vec<SplitPoint>> split(SplitPoint const &sp, vec<SplitPoint> const &pairs) {
    vec<SplitPoint> left, right;
    uint i, j; std::tie(i, j) = sp;
    auto on_left = [&](auto x) {return x <= i || x >= j;};
    auto on_right = [&](auto x) {return x >= i && x <= j;};
    auto to_left = [&](auto x) {return x <= i ? x : x - j + 1 + i;};
    auto to_right = [&](auto x) {return x - i;};

    for (auto const &p : pairs) {
        NUPACK_REQUIRE(crosses(p, sp), ==, false);
        uint d, e; std::tie(d, e) = p;
        if (on_left(d) && on_left(e)) left.emplace_back(to_left(d), to_left(e));
        else right.emplace_back(to_right(d), to_right(e));
    }
    left.emplace_back(to_left(i), to_left(j));
    right.emplace_back(to_right(i), to_right(j));

    return {std::move(left), std::move(right)};
}

}}
