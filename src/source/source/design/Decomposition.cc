#include <nupack/design/Decomposition.h>
#include <nupack/design/DesignParameters.h>

namespace nupack { namespace newdesign {




/**
 * @brief combine the information from two children into the implied parent estimate
 * @details From the perspective of the parent node, this takes the data from
 *     the two children implied by the SplitPoint and maps the pair
 *     probabilities data into the parent pair probabilities matrix and
 *     multiplies the two independent conditional child partition functions
 *     together.
 *
 * @param sp the base pair in the parent which separates the two children
 * @param left the child whose indices are before and after the split point indices
 * @param right the child whose indices are between the split point indices
 * @return pair probabilities and partition function estimate of the parent by the children
 */
ThermoData join_children(SplitPoint sp, ThermoData const &left, ThermoData const &right) {
    auto pfunc = left.second + right.second;
    auto const n = left.first.n_rows + right.first.n_rows - 2;

    ProbabilityMatrix joined_probs(n, n);

    uint i, j; std::tie(i, j) = sp;
    auto from_left  = [&](auto x) {return x <= i ? x : x - i + j - 1;};
    auto from_right = [&](auto x) {return x + i;};

    auto map_probs = [&](auto const &probs, auto f) {
        for (auto it = probs.begin(); it != probs.end(); ++it)
            joined_probs(f(it.row()), f(it.col())) = *it;
        // for (auto i : indices(probs)) for (auto j : indices(probs))
        //     *joined_probs(f(i), f(j)) = *probs(i, j);
    };
    map_probs(left.first, from_left);
    map_probs(right.first, from_right);

    return std::make_pair(std::move(joined_probs), pfunc);
}


/**
 * @brief Sums up the data from alternative children resulting from multiple
 *     exclusive split points
 * @details Sums up the partition functions from each pair of children and
 *     then combines the alternate pair probabilities matrices by weighting
 *     the elements according to the fraction of the total partition function
 *     estimate they contribute
 *
 * @param results the thermodynamic data from each pair of children
 * @return The merged pair probabilities and partition function
 */
ThermoData merge_alternatives(vec<ThermoData> const &results, real f_sparse) {
    /* single decomposition */
    if (len(results) == 1) return results[0];

    real total_pf = accumulate(item_view(results), [](auto &a, auto b) {a = log_sum_exp(a, b);});

    auto const n = results[0].first.n_rows;
    ProbabilityMatrix total_pprobs(n, n);

    for (auto const &res : results) {
        real fraction = std::exp(res.second - total_pf);

        total_pprobs = total_pprobs + fraction * res.first;
        // zip(total_pprobs.storage, res.first.storage, [&](auto &a, auto const &b) {
        //     a += fraction * b;
        // });
    }

    total_pprobs.transform([&](auto val) -> real {
        return (val < f_sparse) ? 0.0 : val;
    });

    return {total_pprobs, total_pf};
}


/**
 * @brief adds the pair of children implied by the SplitPoint to the children
 *     vector
 * @details Used in both structure_decompose() and probability_decompose(),
 *     this function splits the structure, seqeunce, and list of enforced base
 *     pairs appropriately based on sp. The two resulting children are then added to the vec
 *
 * @param sp the SplitPoint implying the two child nodes
 */
void ComplexNode::add_child(SplitPoint sp) {
    auto child_strucs = split(sp, structure);
    auto child_seqs = split(sp, sequence);
    auto child_enforced = split(sp, enforced_pairs);

    children.emplace_back(sp,
            std::make_pair(ComplexNode(child_seqs.first, child_strucs.first, child_enforced.first),
                        ComplexNode(child_seqs.second, child_strucs.second, child_enforced.second)));
}


/**
 * @brief Recursively divide the current node into the lowest cost pair of children
 *     consistent with the minimum size and helix padding parameters
 *
 * @param min_size the minimum number of nucleotides that must be in a child
 *     node for the decomposition to be valid
 * @param min_helix the minimum number of flanking base pairs that must exist
 *     on either side of a SplitPoint for it to be valid
 */
void ComplexNode::structure_decompose(uint min_size, uint min_helix) {
    auto splits = valid_split_points(structure, min_size, min_helix);
    if (splits.empty()) return;
    auto best = ascending_cost_splits(splits, len(*this))[0];

    children.clear();
    add_child(best);

    child_op([&](auto &c) {c.structure_decompose(min_size, min_helix);});
}


/**
 * @brief Probability-guided recursive decomposition with multiple exclusive
 *     split points
 * @details For each node, compute the pair probability matrix given the
 *     Sequence s and use this and the target structure (if available) to
 *     decompose the current node via multiple exclusive split points subject to
 *     the minimum size, minimum helix padding, and minimum flanking probability
 *     constraints. Then, repeat the procedure recursively on the newly
 *     generated child nodes. If a node cannot be decomposed while meeting the
 *     constraints, decomposition stops on this branch.
 *
 * @param min_size the minimum number of nucleotides that must be in a child
 *     node for the decomposition to be valid
 * @param min_helix the minimum number of flanking base pairs that must exist on
 *     either side of a SplitPoint for it to be valid
 * @param f_split the minimum probability that must be captured by the exclusive
 *     split points.
 * @param s the sequence to use when computing pair probabilities
 * @param mods the set of models consistent with the parent complex model
 *     specification to use when computing thermodynamic information
 */
bool ComplexNode::probability_decompose(DecompositionParameters const &params,
        Sequence const &s, ThermoEnviron &t_env, int depth, LevelSpecification const &indiv, EngineObserver &obs) {
    bool revoke_cache = false;
    depth = indiv.get_depth(index, depth);

    if (depth <= 0 || (children.empty() && !indiv)) {
        auto probs = dynamic_program(Local(), t_env, s, 0, params, {}, obs).first;
        auto optimal_splits = minimal_splits(probs, params.f_split, params.N_split, params.H_split, structure);
        // BEEP(key_view(children), optimal_splits);


        auto init_len = len(children);
        erase_if(children, [&](auto const &c) {return !contains(optimal_splits, c.first);});
        revoke_cache = init_len > len(children);

        for (auto const &o : optimal_splits) {
            if (!contains(key_view(children), o)) {
                add_child(o);
                revoke_cache = true;
            }
        }
    }

    child_op([&](auto &c) {
        revoke_cache |= c.probability_decompose(params, s, t_env, depth-1, indiv, obs);
    });
    if (revoke_cache) cache.revoke_non_root();
    return revoke_cache;
}


/**
 * @brief Compute the pair probabilities and partition function of a given
 *     node, or merge the results from its children
 * @details For nodes at depth == 0 (or nodes above this depth with
 *     no subtree beneath them), thermodynamic information is computed. It is
 *     then merged recursively up the tree: all nodes at greater depth join
 *     the results from their left and right children and merge the
 *     alternatives in the case of multiple exclusive split points.
 *
 * @param mods the set of models consistent with the parent complex model
 *     specification to use when computing thermodynamic information
 * @param s the sequence to use when computing thermodynamic information
 * @param depth the number of levels lower the function must merge recursively
 *     before actually computing thermodynamic information
 * @return either the computed value for this node, or the merged results of
 *     its children
 */
ThermoData ComplexNode::dynamic_program(Local env, ThermoEnviron &t_env, Sequence const &s, uint depth,
        DecompositionParameters const &params, LevelSpecification const &indiv, EngineObserver &obs) const {
    // if (!obs.slowdown) throw; /* enable this to debug cases where obs is not getting passed down */
    auto seq = to_nick_sequence(sequence, s);
    auto mods = t_env.doubled();
    bool use_higher_cache = !(indiv); // ignore if there are any exception nodes
    depth = indiv.get_depth(index, depth);

    auto can_pair = [&](auto pair) {
        uint i, j; std::tie(i, j) = pair;
        return std::get<0>(mods).can_pair(seq[i], seq[j]);
    };

    /* if all children are predicated on pairs that can't form, this node should be evaluated */
    if (depth == 0 || children.empty() || none_of(key_view(children), can_pair)) {
    // if (depth == 0 || children.empty()) {
        if (cache.match(seq, 0)) {
            // BEEP("cache hit!");
            return cache.get(0);
        }


        decltype(newdesign::pair_probability(threshold(env, *this), seq, t_env, obs)) ret;
        if (len(enforced_pairs) == 0) {
            ret = newdesign::pair_probability(threshold(env, *this), seq, t_env, obs);
        } else {
            ret = newdesign::pair_probability(threshold(env, *this), seq, mods, enforced_pairs, params.dG_clamp, obs);
        }

        /* convert raw N^2 probability matrix to sparse matrix */
        ThermoData sparse_ret {sparsify(ret.first, params.f_sparse), ret.second};
        cache.add(seq, sparse_ret, 0);
        return sparse_ret;
    } else {
        if (use_higher_cache && cache.match(seq, depth)) {
            // BEEP("cache hit!");
            return cache.get(depth);
        }

        // new

        /* flat vector of all children to evaluate */
        vec<ComplexNode const *> to_evaluate;
        vec<SplitPoint> to_join_on;
        for (auto const &c : children) {
            if (can_pair(c.first)) {
                to_join_on.emplace_back(c.first);
                to_evaluate.emplace_back(&(c.second.first));
                to_evaluate.emplace_back(&(c.second.second));
            }
        }

        /* evaluate all children, potentially in parallel */
        auto eval = [&](auto const &env, auto i) {
            auto const &c = *at(to_evaluate, i);
            return c.dynamic_program(env, t_env, s, depth-1, params, indiv, obs);
        };
        vec<ThermoData> child_results = env.map(len(to_evaluate), 1, eval);

        /* join pairs of children, potentially in parallel */
        auto join = [&](auto const &env, auto i) {
            auto pair = at(to_join_on, i);
            auto const &left = at(child_results, 2*i);
            auto const &right = at(child_results, 2*i+1);
            return join_children(pair, left, right);
        };
        vec<ThermoData> joined_results = env.map(len(to_join_on), 1, join);

        auto ret = merge_alternatives(joined_results, params.f_sparse);

        // end new

        if (use_higher_cache) cache.add(seq, ret, depth);
        // BEEP("cache miss!");
        return ret;
    }
}

/**
 * @brief compute the depth of the subtree beneath (and including) this node
 *
 * @return the subtree depth
 */
uint ComplexNode::depth() const {
    if (children.empty()) return 0;
    uint max_depth = 0;

    child_op([&](auto &c) {
        auto cur_depth = c.depth();
        if (max_depth < cur_depth) max_depth = cur_depth;
    });

    return max_depth + 1;
}

/**
 * @brief recursive function that adds the index of all nodes which are beneath
 * the current node by depth levels
 *
 * @param registered a vector of all indices at depth "depth" beneath this node
 * @param depth the current number of levels above the desired level
 */
void ComplexNode::register_indices(vec<int> &registered, uint depth, bool include_leaves) const {
    if (depth == 0 && (include_leaves || !children.empty())) registered.emplace_back(index);

    if (depth > 0) child_op([&](auto const &c) {
        c.register_indices(registered, depth - 1, include_leaves);
    });
}

}}
