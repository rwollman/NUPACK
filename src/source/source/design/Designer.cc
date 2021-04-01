#include <nupack/design/Designer.h>
#include <nupack/design/OutputResult.h>
#include <nupack/design/Specification.h>
#include <nupack/reflect/Serialize.h>

namespace nupack { namespace newdesign {

void DesignParameters::init_rng() const {
    auto seed = rng_seed;
    std::random_device rd;
    while (seed == 0) seed = rd();
    StaticRNG.seed(seed);
}

/**
 * @brief encapsulates an expression in a lambda so that it only gets evaluated
 * if the value is actually needed (by logging).
 *
 */
#define LAZY(x) [&] {return x;}


/**
 * @brief initializes random sequence consistent with constraint and
 * conditionally decomposes active structures. Creates models necessary for
 * evaluating properties of each complex so that further access is threadsafe.
 *
 * @param decompose whether or not to decompose complexes (to allow root-only
 * design)
 */
void Designer::initialize(bool decompose) {
    parameters.init_rng();

    timer = Timer().start();

    design.initialize_sequence();
    // disable constraint solver switching if deterministic (see comment in Constraints.cc)
    if (parameters.rng_seed != 0) design.sequences.constraints.msec_cutoff = 0;

    /* burn in for constraint co-variation */
    // auto temp_seq = design.sequence();
    // for (auto i : range(design.sequences.real_variables)) design.mutate_sequence({i});
    // design.set_sequence(temp_seq);
    /* end burn in */

    if (decompose) design.initialize_decomposition(Psi);

    // {
    //     for (auto &c : design.complexes) {
    //         c.index_nodes();
    //         for (auto depth : range(c.depth())) {
    //             BEEP(depth);
    //             auto inds = c.get_node_indices(depth);
    //             c.decomposition.apply_recursive([&](auto const &node) {
    //                 if (contains(inds, node.index)) {
    //                     BEEP(len(node), node.structure);
    //                 }
    //             });
    //         }
    //     }
    // }

    /* ensure models are already created before parallel access */
    for (auto const &c : design.complexes) {
        for_each(c.target.cached_models(design.models), [&](auto const& m) {m.reserve(2 * len(c));});
    }

    /* reserve cache for models */
    design.models.create_caches(parameters.cache_bytes_of_RAM);

    for (auto &o : objectives) o.initialize(design);

    max_depth = design.max_depth();

    /* no need to check weights during design if  */
    if (bool(weights)) weights.resolve_weights(design);
}


/**
 * @brief redecompose all active complexes at the given depth
 *
 * @param depth the depth at which redecomposition starts in each complex
 */
void Designer::redecompose_active(Local const &env, uint depth) {
    design.redecompose_active(env, depth, Psi);
    max_depth = design.max_depth();
}


/**
 * @brief decompose parent nodes at a given depth in descending order of the
 * underestimate in the defect resulting from replacing the given node with its
 * children. Does this until either the difference between child and parent
 * defects falls beneath a threshhold (success) or until all parents have
 * attempted to be redecomposed (failure).
 *
 * @param depth the parent depth
 * @param parent the defect estimate computed at the parent depth
 * @param init_child the defect estimate computed one level deeper
 * @return true the child defect after redecomposition is within threshold distance from parent
 * @return false the child defect still underestimates the parent by too much.
 */
bool Designer::redecompose(uint depth, Sequence const &sequence) {
    Local env;

    auto saved_seq = design.sequence();
    design.set_sequence(sequence);

    /* recompute (or pull from cache) normalized defects for */
    auto parent = design.normalized_defect(env, depth, Psi, {}, weights, obs);
    // auto init_child = design.normalized_defect(env, depth+1, Psi, {}, weights, obs);

    auto multitube_position = value_of(find_multitube(objectives));
    auto init_child = at(best.forest, depth+1).defect(multitube_position);

    /* evaluate underestimates caused by replacing each non-leaf node at depth
    "depth" with its children */
    vec<std::pair<std::pair<uint, int>, real>> child_replaced;

    // BEEP(parent.total(), init_child.total());
    // auto yada = [&] {
    //     for (auto i : Psi.actives()) {
    //         auto const &c = at(design.complexes, i);
    //         auto parent_pf = c.log_pfunc(env, design.models, design.sequence(), depth);
    //         auto child_pf = c.log_pfunc(env, design.models, design.sequence(), depth+1);
    //         auto parent_defect = c.defect(env, design.models, design.sequence(), depth).total();
    //         auto child_defect = c.defect(env, design.models, design.sequence(), depth+1).total();
    //         auto fraction_captured = std::exp(child_pf - parent_pf);
    //         auto len_children = len(c.get_node_indices(depth, false));
    //         BEEP(c.name, depth, len_children, parent_pf, child_pf, fraction_captured, parent_defect, child_defect);
    //         if (len_children == 0 && fraction_captured < 1) {
    //             BEEP(c.log_pfunc(env, design.models, design.sequence(), depth));
    //             BEEP(c.decomposition);
    //             BEEP(c.log_pfunc(env, design.models, design.sequence(), depth+1));
    //             BEEP(c.decomposition);
    //         }
    //     }
    //     for (auto const &t : design.tubes) {
    //         auto parent_concs = t.concentrations(env, design.models, design.complexes, design.sequence(), depth, Psi);

    //         vec<string> targets;
    //         for (auto const &tar : t.targets) {
    //             targets.emplace_back(at(design.complexes, tar.complex_index).name);
    //         }

    //         auto child_concs = t.concentrations(env, design.models, design.complexes, design.sequence(), depth+1, Psi);

    //         auto parent_defect = t.normalized_defect(env, design.models, design.complexes, design.sequence(), depth, Psi).total();
    //         auto child_defect = t.normalized_defect(env, design.models, design.complexes, design.sequence(), depth+1, Psi).total();
    //         BEEP(t.name, parent_defect, child_defect);
    //         zip(targets, parent_concs, child_concs, [](auto target, auto parent, auto child) {
    //              BEEP(target, parent, child);
    //         });
    //     }
    // };
    // yada();

    for (auto index: Psi.actives()) {
        auto &c = at(design.complexes, index);
        c.index_nodes();
        auto nodes = c.get_node_indices(depth, false);

        for (auto node: nodes) {
            LevelSpecification indiv_spec;
            indiv_spec.add_exception(node, 1);

            EnsembleLevelSpecification ens_spec;
            ens_spec.add_level_spec(index, indiv_spec);

            real underestimate = parent.total() - design.normalized_defect(env, depth, Psi, ens_spec, weights, obs).total();
            child_replaced.emplace_back(std::pair<uint, int>{index, node}, underestimate);

            // { // testing block
            //     BEEP(index, node, c.name);
            //     c.decomposition.apply_recursive([&](auto const &n) {
            //         if (n.index == node) BEEP(len(n.children));
            //     });
            // }
        }
    }

    if (len(child_replaced) == 0) {
        return true;
    }

    /* sort nodes in descending order of underestimate */
    sort(child_replaced, [](auto const &a, auto const &b) {return a.second > b.second;});

    auto cur = begin_of(child_replaced);

    auto child_defect = design.normalized_defect(env, depth+1, Psi, {}, weights, obs);
    real cutoff = parameters.f_redecomp * (parent.total() - init_child.total() / parameters.f_stringent);
    // BEEP(child_replaced);
    auto condition = [&](auto const &x) {
        // BEEP(parent.total(), init_child.total(), x.total());
        return (parent.total() - x.total() / parameters.f_stringent) > cutoff;
    };

    bool any_changed = false;
    std::set<uint> changed_complex_inds;
    while (cur != end_of(child_replaced) && condition(child_defect)) {
        // decompose node implied by cur
        uint comp_index;
        int node_index;
        std::tie(comp_index, node_index) = cur->first;

        auto underestimate = cur->second;
        auto &c = at(design.complexes, comp_index);

        LevelSpecification spec;
        spec.add_exception(node_index, 0);

        bool changed = c.probability_decompose(design.sequence(), design.models, c.depth()+1, spec, obs);
        any_changed |= changed;

        if (changed) changed_complex_inds.emplace(comp_index);

        child_defect = design.normalized_defect(env, depth+1, Psi, {}, weights, obs);

        logs.log("basic", time_elapsed(), "redecomposed", depth+1, Psi.num_active(), Psi.num_inactive(),
                design.sequences.json_domains(), LAZY(child_defect.total()));
        ++cur;
    }

    max_depth = design.max_depth();

    /* restore original sequence */
    design.set_sequence(saved_seq);

    /* log changed decompositions */

    for (auto i : changed_complex_inds) {
        auto const &c = at(design.complexes, i);
        logs.log("decomposition", i, c.name, LAZY(c.json_decomposition()));
    }

    // BEEP(condition(child_defect), parent.total(), child_defect.total());
    // yada();
    return !condition(child_defect);
}


/**
 * @brief decompose a subset of all complexes at a given depth
 *
 * @param subset indices into complexes vector of complexes that should be decomposed
 * @param depth the depth at which decomposition starts in each complex
 */
void Designer::subset_decompose(vec<uint> subset, uint depth) {
    for (auto c : subset) {
        at(design.complexes, c).probability_decompose(design.sequence(), design.models, depth, {}, obs);
    }

    for (auto i : subset) {
        auto const &c = at(design.complexes, i);
        logs.log("decomposition", i, c.name, LAZY(c.json_decomposition()));
    }

    max_depth = design.max_depth();
}


/**
 * @brief add off-targets to active set until difference between full ensemble
 * defect and focused estimate is small enough
 *
 * @param env compute resources allowing for potential parallel execution
 * @param full the defect calculated for the full ensemble
 * @param init_estimate the defect calculated for the focused estimate at the
 * root level before adding more off-targets during this refocus operation
 */
void Designer::refocus(Local const &env, Sequence const &sequence) {
    auto saved_seq = design.sequence();
    design.set_sequence(sequence);

    auto full = design.normalized_defect(env, 0, {}, {}, weights, obs);
    auto init_estimate = design.normalized_defect(env, 0, Psi, {}, weights, obs);
    if (Psi.all_active()) NUPACK_ERROR("can't refocus if all complexes are already active", full, init_estimate);

    /* determine the order to add off-targets in Psi_passive into Psi_active
      based on fractional contribution to concentration defect */
    vec<real> fractions(len(design.complexes), 0.0);
    auto log_pfuncs = design.log_pfuncs(env, 0, {}, {}, obs);
    for (auto const &tube : design.tubes) {
        zip(tube.targets, tube.fractions(log_pfuncs),
        [&] (auto const &c, auto frac) {
            auto i = c.complex_index;
            if (!Psi.active(i)) at(fractions, i) += frac;
        });
    }

    vec<std::pair<uint, real>> passive;
    izip(Psi.mask, [&](auto i, auto b) {if (!b) passive.emplace_back(i, at(fractions, i));});
    sort(passive, [](auto const &a, auto const &b) {return a.second > b.second;});

    auto order = key_view(passive);
    auto cur = begin_of(order);
    auto part = Psi;

    if (cur == end_of(order)) NUPACK_ERROR("first passive complex to add out of range",
            part, len(part), cur, passive, len(passive),
            full.total(), init_estimate.total());
    at(part.mask, *cur) = true;

    auto estimate = design.normalized_defect(env, 0, part, {}, weights, obs);
    // {design.sequence(), {design.normalized_defect(env, 0, part)}};
    logs.log("basic", time_elapsed(), "refocused", 0, part.num_active(), part.num_inactive(), design.sequences.json_domains(), LAZY(estimate.total()));


    real cutoff = parameters.f_refocus * (full.total() - init_estimate.total());
    while (full.total() - estimate.total() > cutoff) {
        // if (cur >= end_of(order)) break;
        at(part.mask, *(++cur)) = true;
        estimate = design.normalized_defect(env, 0, part, {}, weights, obs);
        // {design.sequence(), {design.normalized_defect(env, 0, part)}};
        logs.log("basic", time_elapsed(), "refocused", 0, part.num_active(), part.num_inactive(), design.sequences.json_domains(), LAZY(estimate.total()));
    }

    vec<uint> changed;
    izip(part.mask, Psi.mask, [&](auto i, auto n, auto o) {if (n && !o) changed.emplace_back(i);});
    subset_decompose(Psi.actives());
    stats.offtargets_added_per_refocus.emplace_back(len(changed));

    Psi = part;
    known_bads.clear();


    /* restore original sequence */
    design.set_sequence(saved_seq);
}


/**
 * @brief top-level entry into design algorithm. the main loop at this level
 * checks whether the root-level, full-ensemble multitube ensemble defect is
 * better than f_stop or better than the focused, root-level estimate.
 *
 * @param env compute resources allowing for potential parallel execution
 * @return the best discovered root-level and full-ensemble defect
 */
Result Designer::optimize_tubes(Local const &env) {
    /* print headers for CSV log files */
    logs.log("basic", "time", "type", "depth", "psi_active", "psi_passive", "sequence", "defect");
    // logs.log("thermo", "type", "length", "time", "cache possible");
    obs.log("thermo", "type", "length", "time", "cache possible");
    logs.log("decomposition", "index", "name", "decomposition");

    /* initial logging of active decompositions */
    for (auto i : Psi.actives()) {
        auto const &c = at(design.complexes, i);
        logs.log("decomposition", i, c.name, LAZY(c.json_decomposition()));
    }

    // return alternate_optimize_tubes(env);
    return optimize_tubes_impl(env);
}


Result Designer::optimize_tubes_impl(Local const &env) {
    /******************************************************************/
    max_depth = design.max_depth();

    auto estimate = optimize_forest(env, design.sequence());
    design.set_sequence(estimate.sequence);
    Result full = reevaluate_objectives(env, estimate, 0, {}, weights);
    // {design.sequence(), {design.normalized_defect(env)}};
    if (full.weighted_total() < best.full.weighted_total()) best.full = full;
    logs.log("basic", time_elapsed(), "root accepted", 0, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(full.sequence), LAZY(full.weighted_total()));


    while (full.weighted_total() > max(parameters.f_stop, estimate.weighted_total())) {
        checkpoint(*this, false);
        refocus(env, full.sequence);
        estimate = optimize_forest(env, full.sequence);
        design.set_sequence(estimate.sequence);
        full = reevaluate_objectives(env, estimate, 0, {}, weights);
        // {design.sequence(), {design.normalized_defect(env)}};
        if (full.weighted_total() < best.full.weighted_total()) {
            logs.log("basic", time_elapsed(), "root accepted", 0, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(full.sequence), LAZY(full.weighted_total()));
            best.full = full;
        } else {
            logs.log("basic", time_elapsed(), "root rejected", 0, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(full.sequence), LAZY(full.weighted_total()));
        }
    }
    checkpoint(*this, true);

    /******************************************************************/
    stats.design_time += timer.stop(); // if checkpointed and restarted, the += will make the output stats reflect the total design time instead of just the most recent segment
    stats.final_Psi = Psi;

    if (parameters.time_analysis) time_analysis(env);

    // std::cout << json(design.sequences.correlation_matrix) << std::endl;

    return best.full;
}




/**
 * @brief measure time taken to compute partition function and pair probabilities for all complexes
 *
 * @param env compute resources allowing for potential parallel execution
 */
void Designer::time_analysis(Local const &env) {
    for (auto &c : design.complexes)
        c.decomposition.apply_recursive([](auto &node) {
            node.cache = ComplexNode::Cache();
        });

    /* to force actually recomputing everything */
    design.models.clear_caches();

    auto t = Timer().start();
    (void) design.normalized_defect(env, 0, {}, {}, weights, obs);
    stats.analysis_time = t.stop();

    // auto eval = [&]() {
    //     for (auto &c: design.complexes) (void) c.pair_probabilities(env, design.models, design.sequence());
    //     // (void) design.normalized_defect(env);
    // };
    // stats.analysis_time = time_it(30, eval);
}


/**
 * @brief manages merging decomposed estimates of the defect after
 * leaf-optimization finishes, either accepting merges or calling for
 * redecomposition to improve the estimates at lower levels. Exits when
 * root-level estimate meets stop condition or doesn't appreciably add defect to
 * level beneath
 *
 * @param env compute resources allowing for potential parallel execution
 * @param seq the starting sequence to pass on to leaf optimization
 * @return the best discovered root-level, focused ensemble defect estimate
 */
Result Designer::optimize_forest(Local const &env, Sequence seq) {
    // archive.reset(archive.forest); archive.resize_forest(max_depth + 1);

    best.reset(best.forest); best.resize_forest(max_depth + 1);
    at(best.forest, max_depth).sequence = seq;

    bool merge_successful = false;

    while (!merge_successful) {
        auto & leaf_best = at(best.forest, max_depth);
        leaf_best = optimize_leaves(env, leaf_best.sequence);
        design.set_sequence(leaf_best.sequence);

        // auto blah = at(archive.forest, max_depth).merge(archive.leaf_opt);
        // BEEP(max_depth, blah);

        auto depth = max_depth - 1;
        merge_successful = true;
        while (depth >= 0 && merge_successful) {
            Result cur_result = reevaluate_objectives(env, at(best.forest, depth + 1), depth, Psi, weights);
            // {design.sequence(), {design.normalized_defect(env, depth, Psi)}};
            // auto temp_archive = at(archive.forest, depth+1);
            // temp_archive.update_estimates(env, *this, depth, Psi);
            // auto blah = at(archive.forest, depth).merge(temp_archive);
            // BEEP(depth, blah);

            if (cur_result.weighted_total() < at(best.forest, depth).weighted_total()) {
                at(best.forest, depth) = cur_result;
                logs.log("basic", time_elapsed(), "best merge", depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(cur_result.sequence), LAZY(cur_result.weighted_total()));
            }

            auto f_d_stop = parameters.f_stop * pow(parameters.f_stringent, depth);
            auto child_defect = at(best.forest, depth + 1).weighted_total();

            if (cur_result.weighted_total() > std::max(f_d_stop, child_defect / parameters.f_stringent)) {
                checkpoint(*this, false);
                logs.log("basic", time_elapsed(), "merge unsuccessful", depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(cur_result.sequence), LAZY(cur_result.weighted_total()));

                merge_successful = false;
                redecompose(depth, at(best.forest, depth + 1).sequence);

                /* log decomposition failure at current level */
                if (len(stats.num_redecompositions) <= max_depth) stats.num_redecompositions.resize(max_depth + 1, 0);
                ++at(stats.num_redecompositions, depth);

                best.resize_forest(max_depth + 1);
                for_each(view(best.forest, depth + 1, max_depth + 1), [&](auto &b) {b = inf_result;});

                // archive.resize_forest(max_depth + 1);
                // for_each(view(archive.forest, depth + 1, max_depth + 1), [&](auto &b) {b = archive.dfault;});

                at(best.forest, max_depth).sequence = design.sequence();
                known_bads.emplace(design.sequence());
            } else {
                logs.log("basic", time_elapsed(), "merge successful", depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(cur_result.sequence), LAZY(cur_result.weighted_total()));
            }

            --depth;
        }
    }

    return at(best.forest, 0);
}


/**
 * @brief manages reseeding in the case that leaf mutation fails to reach leaf
 * stop condition. Exits with best encountered leaf-level defect sequence.
 *
 * @param env compute resources allowing for potential parallel execution
 * @param seq the starting sequence to pass on to leaf mutation
 * @return the best leaf-level defect estimate and sequence
 */
Result Designer::optimize_leaves(Local const &env, Sequence seq) {
    // archive.reset(archive.leaf_opt);
    best.leaf_opt = mutate_leaves(env, seq);
    // auto blah = archive.leaf_opt.merge(archive.leaf_mut);
    // BEEP(blah);

    uint m_reopt = 0;
    auto f_D_stop = parameters.f_stop * pow(parameters.f_stringent, max_depth);
    while (best.leaf_opt.weighted_total() > f_D_stop && m_reopt < parameters.M_reopt) {
        checkpoint(*this, false);
        /* reseed from best sequence */
        design.set_sequence(best.leaf_opt.sequence);
        auto sampled_nucs = scalarized_sample(best.leaf_opt, parameters.M_reseed);
        bool mutation_succeeded = design.mutate_sequence(sampled_nucs);

        if (!mutation_succeeded) {
            ++m_reopt;
            continue;
        }

        Result temp = evaluate_objectives(env, max_depth, Psi, weights);
        logs.log("basic", time_elapsed(), "reseeded", max_depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(), LAZY(temp.weighted_total()));
        ++stats.num_reseeds;

        auto cur_result = mutate_leaves(env, design.sequence());
        // auto blah = archive.leaf_opt.merge(archive.leaf_mut);
        // BEEP(blah);
        if (cur_result.weighted_total() < best.leaf_opt.weighted_total()) {
            best.leaf_opt = cur_result;
            m_reopt = 0;
        } else {
            ++m_reopt;
        }
    }
    return best.leaf_opt;
}


/**
 * @brief starting from the provided initial sequence, attempts to find a
 * sequence with a defect estimate less than the leaf stop condition through
 * directed single nucleotide variable mutation. Exits when such a sequence is
 * found or when enough failed sequences are evaluated without finding an
 * improving sequence.
 *
 * @param env compute resources allowing for potential parallel execution
 * @param seq the starting sequence to pass to begin mutation from
 * @return the sequence with the best encountered leaf-level defect estimate
 */
Result Designer::mutate_leaves(Local const &env, Sequence seq) {
    // archive.reset(archive.leaf_mut);
    /**
     * \gamma_{bad} in pseudocode. Here we initialize with sequences which
     * would otherwise cause cycling following redecomposition.
     */
    std::set<Sequence> bad_seqs = known_bads;
    // vec<Sequence> bad_seqs;

    design.set_sequence(seq);
    best.leaf_mut = evaluate_objectives(env, max_depth, Psi, weights);
    // archive.leaf_mut.attempt_add(best.leaf_mut);
    // {design.sequence(), {design.normalized_defect(env, max_depth, Psi)}};
    ++stats.num_leaf_evaluations;


    uint m_bad = 0;
    while (nupack::contains(bad_seqs, best.leaf_mut.sequence) && m_bad < parameters.M_bad) {
        auto sampled_nucs = scalarized_sample(best.leaf_mut);
        design.mutate_sequence(sampled_nucs);
        best.leaf_mut = evaluate_objectives(env, max_depth, Psi, weights);
        // archive.leaf_mut.attempt_add(best.leaf_mut);
        // {design.sequence(), {design.normalized_defect(env, max_depth, Psi)}};

        ++stats.num_leaf_evaluations;
        ++m_bad;
    }

    logs.log("basic", time_elapsed(), "mutation accepted", max_depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(best.leaf_mut.sequence), LAZY(best.leaf_mut.weighted_total()));
    m_bad = 0;
    auto f_D_stop = parameters.f_stop * pow(parameters.f_stringent, max_depth);

    uint num_muts = 0;
    vec<uint> muts = {0};
    vec<real> defects = {best.leaf_mut.weighted_total()};

    while (best.leaf_mut.weighted_total() > f_D_stop && m_bad < parameters.M_bad) { // && !improvement_slowing(muts, defects)) {
        checkpoint(*this, false);
        /* mutate away from best encountered sequence */
        design.set_sequence(best.leaf_mut.sequence);
        auto sampled_nucs = scalarized_sample(best.leaf_mut);
        bool mutation_succeeded = design.mutate_sequence(sampled_nucs);

        if (nupack::contains(bad_seqs, design.sequence()) || !mutation_succeeded) {
            ++m_bad;
        } else {
            // auto cur_defect = design.normalized_defect(env, max_depth, Psi);
            auto cur_result = evaluate_objectives(env, max_depth, Psi, weights);
            // archive.leaf_mut.attempt_add(cur_result);
            ++stats.num_leaf_evaluations;

            ++num_muts;
            if (cur_result.weighted_total() < best.leaf_mut.weighted_total()) {
                best.leaf_mut = cur_result;
                logs.log("basic", time_elapsed(), "mutation accepted", max_depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(best.leaf_mut.sequence), LAZY(cur_result.weighted_total()));

                muts.emplace_back(num_muts);
                defects.emplace_back(cur_result.weighted_total());

                /* reset counter and tabu */
                bad_seqs = known_bads;
                bad_seqs.clear();
                m_bad = 0;

            } else {
                logs.log("basic", time_elapsed(), "mutation rejected", max_depth, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(), LAZY(cur_result.weighted_total()));
                bad_seqs.emplace(design.sequence());
                ++m_bad;
            }
        }
    }

    return best.leaf_mut;
}

Sequence Designer::best_sequence(Local const &env) {
    /* save initial sequence state */
    auto temp = design.sequence();

    Sequence cur_best = best.leaf_mut.sequence;
    auto compete = [&](auto &other, uint depth, EnsemblePartition const &part={}) {
        design.set_sequence(cur_best);
        auto cur_result = evaluate_objectives(env, depth, part, weights);
        if (cur_result.weighted_total() < other.weighted_total()) {
            other = cur_result; // update with better checkpointing information
        } else {
            cur_best = other.sequence;
        }
    };

    compete(best.leaf_opt, max_depth, Psi);
    for (auto depth : ~indices(best.forest)) compete(at(best.forest, depth), depth, Psi);
    compete(best.full, 0);

    /* return to initial sequence state */
    design.set_sequence(temp);

    return cur_best;
}

bool Designer::improvement_slowing(vec<uint> const &x, vec<real> const &y) {
    real threshold_slope = -0.0001;
    uint max_allowed = 1;

    if (len(x) > max_allowed) {
        auto n = len(x);
        /* too few mutations to say it was truly poor improvement */
        if ((x[n-1] - x[n-2]) < 20) return false;

        auto count = 0;
        for (auto i : range(len(x) - max_allowed, len(x)))
            if ((y[i] - y[i-1]) / (x[i] - x[i-1]) > threshold_slope)
                ++count;
        return !(count < max_allowed);
    }
    return false;
}


Result Designer::evaluate_objectives(Local const &env, uint depth, EnsemblePartition const &part, Weights const &weights) {
    auto seq = design.sequence();
    vec<Defect> defects;
    zip(objectives, [&](auto const &o) {
        defects.emplace_back(o.evaluate(env, design, depth, part, weights, obs));
    });
    return {seq, defects, weights.objective_weights};
}


Result Designer::reevaluate_objectives(Local const &env, Result const &res, uint depth, EnsemblePartition const &part, Weights const &weights) {
    /* store sequence state */
    auto seq = design.sequence();

    design.set_sequence(res.sequence);
    vec<Defect> defects;
    zip(objectives, res.defects, [&](auto const &o, auto const &orig_defect) {
        auto reeval_defect = o.reevaluate(env, design, depth, part, weights, obs);
        auto defect = bool(reeval_defect) ? value_of(reeval_defect) : orig_defect;
        defects.emplace_back(defect);
    });
    Result ret{res.sequence, std::move(defects), weights.objective_weights};
    if (depth == 0 && part.all_active()) ret.full_evaluation(*this);

    /* restore sequence state */
    design.set_sequence(seq);

    return ret;
}


}

}
