#include <nupack/design/Complex.h>
#include <nupack/design/ThermoWrapper.h>

namespace nupack { namespace newdesign {

/**
 * @brief Compute the logarithm of the partition function for the complex
 * @details compute the logarithm of the partition function taking into account symmetry and join penalties
 *
 * @param env the compute resources to pass to the thermo code.
 * @param map [description]
 * @param s [description]
 *
 * @return [description]
 */
real Complex::log_pfunc(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs) const {
    env = threshold(env, *this);
    return newdesign::partition_function(env, to_nick_sequence(strands, s), target.environment(map), obs)
         - symmetry_correction()
         - target.model.beta * bonus;
};

#if 1

real Complex::log_pf_single_strands(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs) const {
    auto seqs = to_nick_sequence(strands, s);
    real log_pf = 0;
    auto &t_env = target.environment(map);

    for (auto const &s : seqs.strands()) {
        ::nupack::Complex seq {s};
        auto val = t_env.get_pfunc(seq);
        if (std::holds_alternative<real>(val)) {
            log_pf += std::get<real>(val);
        }
        else {
            auto temp = newdesign::partition_function(env, seq, t_env, obs);
            t_env.add_pfunc(seq, temp);
            log_pf += temp;
        }
    }

    return log_pf - symmetry_correction() + join_penalty(map);
}

#endif


/**
 * @brief compute the pair probabilities matrix for a complex
 * @details compute the pure pair probabilities matrix (no enforced base
 *     pairs) for a target complex
 *
 * @param env the compute resources to pass to the thermo code.
 * @param map the model cache to draw the model to evaluate with from
 * @param s the current sequence from which to pull the sequence for the
 *     complex to evaluate
 *
 * @return The pair probabilities matrix
 */
Tensor<real, 2> Complex::pair_probabilities(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs) const {
    /* skip entirely if no target structure as it won't be used */
    if (!is_on_target()) return {};
    env = threshold(env, *this);
    return newdesign::pair_probability(env, to_nick_sequence(strands, s), target.environment(map), obs).first;
};


/**
 * @brief computes the complex ensemble defect for the complex
 * @details wrapper function taking a complex, evaluating and then forwarding
 *     the pair probabilities and target structure to the function for actual
 *     calculation.
 *
 * @param env the compute resources to pass to the thermo code.
 * @param complex the complex to evaluate
 * @param s the sequence to map onto the complex
 *
 * @return A defect object
 */
Defect Complex::defect(Local env, ModelMap const &map, Sequence const &s, EngineObserver &obs) const {
    /* skip if no target structure */
    if (!target.has_structure()) return Defect();

    auto defects = nucleotide_defects(pair_probabilities(env, map, s, obs), target.structure);

    vec<real> mapped_defects(len(s), 0.0);
    zip(to_indices(), defects, [&] (auto i, auto d) {
        mapped_defects[i] += d;
    });

    defect_vec defs;
    izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d);});

    return {defs};
}


/**
 * @brief the partition function estimate at a given depth
 * @details computes the thermodynamic data at a given depth and then
 *     recursively merges this information up the decomposition tree. The
 *     symmetry correction and join penalty (which are left off to get the
 *     correct fragment answers during node evaluation) are added on the
 *     result of the tree computation.
 *
 * @param env the compute resources to pass to the thermo code.
 * @param map the cache of models to draw the evaluation model from
 * @param s the sequence to use when computing thermodynamic information
 * @param depth the level in the tree
 * @return the symmetry- and join penalty-corrected partition function estimate
 */
real Complex::log_pfunc(Local env, ModelMap const &map, Sequence const &s,
        uint depth, LevelSpecification const &indiv, EngineObserver &obs) const {
    auto ret = decomposition.dynamic_program(env, target.environment(map), s, depth, params, indiv, obs).second;
    if (std::isnan(ret)) {
        auto info = const_cast<Complex *>(this)->hierarchical_pfunc(map, s, depth);
        NUPACK_ERROR("partition function estimate is NaN", to_nick_sequence(strands, s), depth, info);
    }
    // auto extra = (depth == 0 && !(indiv)) ? real(0) : join_penalty(map);
    // return ret - symmetry_correction() + extra;
    return ret - symmetry_correction() + join_penalty(map) - target.model.beta * bonus;
}


/**
 * @brief the pair probabilities estimate at a given depth
 * @details computes the thermodynamic data at a given depth and then
 *     recursively merges this information up the decomposition tree. The
 *     symmetry correction and join penalty (which are left off to get the
 *     correct fragment answers during node evaluation) are added on the
 *     result of the tree computation.
 *
 * @param env the compute resources to pass to the thermo code.
 * @param map the cache of models to draw the evaluation model from
 * @param s the sequence to use when computing thermodynamic information
 * @param depth the level in the tree
 * @return the symmetry- and join penalty-corrected pair probabilities matrix estimate
 */
ProbabilityMatrix Complex::pair_probabilities(Local env, ModelMap const &map, Sequence const &s,
        uint depth, LevelSpecification const &indiv, EngineObserver &obs) const {
    return decomposition.dynamic_program(env, target.environment(map), s, depth, params, indiv, obs).first;
}


/**
 * @brief compute an approximation to the defect at the given depth
 * @details computes the estimate of the pair probabilities matrix at the given
 *      depth and then forwards this and the target structure to the function
 *      computing per nucleotide defects
 *
 * @param env the compute resources to pass to the thermo code.
 * @param map the cache of models to draw the evaluation model from
 * @param s the sequence to use when computing thermodynamic information
 * @param depth the level in the tree
 * @return the defect estimate per nucleotide
 */
Defect Complex::defect(Local env, ModelMap const &map, Sequence const &s,
        uint depth, LevelSpecification const &indiv, EngineObserver &obs) const {
    if (!target.has_structure()) return Defect();

    auto defects = nucleotide_defects(pair_probabilities(env, map, s, depth, indiv, obs), target.structure);

    // vec<real> mapped_defects(len(s), 0.0);
    // zip(to_indices(), defects, [&] (auto i, auto d) {mapped_defects[i] += d;});

    // defect_vec defs;
    // izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d);});
    defect_vec defs;
    zip(to_indices(), defects, [&] (auto i, auto d) {defs.emplace_back(i, d);});

    return {defs};
}


/**
 * @brief decompose the complex according to its target structure
 * @details recursively decompose the structure from the root downwards,
 *     stopping only when a node cannot be split into children that meet the
 *     minimum size and minimum helix padding requirements.
 *
 * @param min_size the minimum number of nucleotides that must be in a child
 *     node for the decomposition to be valid
 * @param min_helix the minimum number of flanking base pairs that must exist
 *     on either side of a SplitPoint for it to be valid
 */
void Complex::structure_decompose() {
    decomposition.structure_decompose(params.N_split, params.H_split);
}


/**
 * @brief decompose the complex according to its pair probabilities matrix
 *     (and optionally, its target structure)
 * @details recursively decompose the complex from the root downwards using
 *     multiple exclusive split points, stopping only when a node cannot be
 *     split into children that meet the minimum size, minimum helix padding,
 *     and minimum captured partition function fraction requirements.
 *
 * @param min_size the minimum number of nucleotides that must be in a child
 *     node for the decomposition to be valid
 * @param min_helix the minimum number of flanking base pairs that must exist
 *     on either side of a SplitPoint for it to be valid
 * @param f_split the minimum partition function fraction that must be
   captured by the exclusive split points
 * @param s the sequence to use when computing pair probabilities
 * @param map the cache of models to draw the evaluation model from
 */
bool Complex::probability_decompose(Sequence const &s, ModelMap const &map,
        uint depth, LevelSpecification const &indiv, EngineObserver &obs) {
    return decomposition.probability_decompose(params, s, target.environment(map), depth, indiv, obs);
}

/**
 * @brief label each node in the complex's decomposition tree with a unique index.
 *
 */
void Complex::index_nodes() {
    int i = 0;

    auto func = [&](auto &c) {c.index = i++;};
    decomposition.apply_recursive(func);
}

/**
 * @brief return a vector of the indices of all the nodes at the given depth
 *
 * @param depth the depth in the decomposition tree at which to find nodes
 * @return the indices of all the nodes at the given depth
 */
vec<int> Complex::get_node_indices(uint depth, bool include_leaves) const {
    vec<int> ret;
    decomposition.register_indices(ret, depth, include_leaves);
    return ret;
}


/**
 * @brief For debugging. Compute the log_pfunc at a given level of the
 * decomposition tree and produce a string representation of the value of every
 * node hit during the calculation.
 *
 * @param map the model cache to draw the model to evaluate with from
 * @param s the current sequence from which to pull the sequence for the complex
 * to evaluate
 * @param depth the depth at which to ultimately evaluate at
 * @return string information about every node encountered during evaluation of
 * the estimate at a given depth
 */
string Complex::hierarchical_pfunc(ModelMap const &map, Sequence const &s, uint depth, EngineObserver &obs) {
    index_nodes();
    auto &ms = target.environment(map);
    std::ostringstream ss;
    Local env;

    for (auto i : range(depth + 1)) {
        auto nodes = get_node_indices(i);
        auto cur_depth = depth - i;
        print_os(ss, i, nodes, cur_depth);
        auto f = [&] (auto &x) {
            if (contains(nodes, x.index)) {
                auto nickseq = to_nick_sequence(x.sequence, s);
                print_os(ss, x.index, x.enforced_pairs, x.sequence, x.structure, nickseq, x.dynamic_program(env, ms, s, cur_depth, params, {}, obs).second, newdesign::partition_function(env, nickseq, ms, obs));
            }
        };
        decomposition.apply_recursive(f);
    }

    return ss.str();
}


string Complex::decomposition_connectivity() const {
    std::ostringstream ss;
    decomposition.apply_recursive([&] (auto const &node) {
        ss << node.index << " -> (";
        node.child_op([&] (auto const &child) {ss << child.index << ", ";});
        ss << ")" << std::endl;
    });
    return ss.str();
}


string Complex::json_decomposition() const {
    std::stringstream ss;
    ss << json(decomposition) << std::flush;
    return ss.str();
}

}}
