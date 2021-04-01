#include <nupack/design/Design.h>

namespace nupack { namespace newdesign {

/**
 * @brief construct a new complex in the design
 * @details If the complex is an on-target in some tube (it has a valid
 *     structure), this structure will be added as its target. Otherwise, a
 *     target without a pairlist vector will be augmented to included a nicks
 *     vector with positions from the sequences.
 *
 * @param strands names of strands in the linearized circular order they appear
 * in the complex
 * @param key description of model for this complex
 * @param name optional name
 * @param struc optional target structure
 */
void Design::add_complex(vec<string> const &strands, Model<real> key,
        string const &name, Structure struc, DecompositionParameters params, real bonus) {
    auto seq = vmap<vec<StrandView>>(strands, [&](auto const &s) {return sequences.get_strand(s);});
    if (!struc.valid()) struc.nicks = prefixes<Nicks>(false, indirect_view(seq, len));
    complexes.emplace_back(std::move(seq), Target(std::move(key), struc), name, std::move(params), bonus);
}

/**
 * @brief Add a tube defined by the complexes implied by indices with the target
 * concentration in concs
 * @details [long description]
 *
 * @param indices the positions in the complexes member of the included
 * complexes in the tube
 * @param concs target concentrations for each of the complexes
 * @param name a name for the tube
 */
void Design::add_tube(vec<uint> const &indices, vec<real> const &concs, string const &name) {
    vec<TubeTarget> targets;
    zip(indices, concs, [&](auto complex, auto target_conc) {targets.emplace_back(complex, target_conc);});
    tubes.emplace_back(targets, name, complexes);
}

/**
 * @brief For every target structure, add a complementarity constraint for every
 * base-pair in the structure.
 *
 */
void Design::add_structure_complementarity() {
    for_each(complexes, [&](auto const &c) {
        if (c.is_on_target()) {
            auto indices = c.to_indices();
            c.target.structure.for_each_pair([&](auto i, auto j) {
                auto var_i = at(indices, i), var_j = at(indices, j);
                // bool duplicate = any_of(sequences.constraints.get_constraints(), [&](auto const &c) {
                //     auto ptr = maybe_get<CompConstraint>(c);
                //     if (ptr == nullptr) return false;
                //     auto vars = ptr->get_constrained_vars();
                //     if (contains(vars, var_i) && contains(vars, var_j)) return true;
                //     return false;
                // });
                // if (!duplicate) sequences.constraints.add_constraint(CompConstraint(var_i, var_j, NUPACK_CS_STRONG));
                sequences.constraints.complementarity_constraint(var_i, var_j, sequences.wobble_mutations);
            });
        }
    });
}


vec<real> Design::log_pfuncs(Local const &env, uint depth,
        EnsemblePartition const &part, EnsembleLevelSpecification const &indiv,
        EngineObserver &obs) const {
    auto compute_log_pfuncs = [&](auto const &env, auto j) -> real {
        if (len(part) == 0 || part.active(j)) {
            auto const &c = at(complexes, j);
            return c.log_pfunc(env, models, sequence(), depth, indiv.get_level_spec(j), obs);
        } else {
            return nan("");
        }
    };

    return env.map(len(complexes), 1, compute_log_pfuncs, AffinitySplit());
}


vec<Defect> Design::complex_defects(Local const &env, uint depth,
        EnsemblePartition const &part, EnsembleLevelSpecification const &indiv,
        EngineObserver &obs) const {
    auto compute_defects = [&](auto const &env, auto j) -> Defect {
        if (len(part) == 0 || part.active(j)) {
            auto const &c = at(complexes, j);
            return c.defect(env, models, sequence(), depth, indiv.get_level_spec(j), obs);
        } else {
            return {};
        }
    };

    return env.map(len(complexes), 1, compute_defects, AffinitySplit());
}

/**
 * @brief compute the nucleotide contributions to the normalized multitube ensemble defect
 *
 * @param depth what depth of the forest to evaluate
 * @param part partition of complexes into \f$\Psi^{\text{active}}\f$ and \f$\Psi^{\text{passive}}\f$
 * @return the normalized multitube ensemble defect
 */
Defect Design::normalized_defect(Local const &env, uint depth,
        EnsemblePartition const &part, EnsembleLevelSpecification const &indiv,
        Weights const &weights, EngineObserver &obs) const {
    // auto compute_active = [&](auto const &env, auto j, auto) {
    //     if (len(part) == 0 || part.active(j)) {
    //         auto const &c = at(complexes, j);
    //         (void) c.pair_probabilities(env, models, sequence(), depth, indiv.get_level_spec(j), obs);
    //     }
    // };
    // env.spread(vec<uint>(indices(complexes)), 1, compute_active, AffinitySplit());

    vec<real> lpfs = log_pfuncs(env, depth, part, indiv, obs);
    vec<Defect> cdefs = complex_defects(env, depth, part, indiv, obs);

    auto f = [&](auto const &env, auto i) {
        auto const &t = at(tubes, i);
        auto tube_weights = bool(weights) ? weights.per_tube.at(i) : ComplexWeights();
        vec<real> defects(len(sequence()), 0.0);
        auto tube_defect = t.normalized_defect(lpfs, cdefs, part, tube_weights);
        for (auto const &d: tube_defect.contributions) defects[d.first] += d.second;
        return defects;
    };

    vec<real> mapped_defects = env.map_reduce(len(tubes), 1, f, sum_vec);

    /* scale contributions by number of tubes */
    auto num_tubes = len(tubes);
    for (auto &m: mapped_defects) m /= num_tubes;

    /* repackage vector of defects into vector of pairs from nucleotides to defects */
    defect_vec defs;
    izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d);});
    return Defect(defs);
}


/**
 * @brief decompose complexes in the active set: structure-based (single split
 * points) if complex has target structure, probability-based (potential
 * multiple exclusive split points) if not.
 *
 * @param min_size the minimum number of nucleotides that must be in a child
 *     node for the decomposition to be valid
 * @param min_helix the minimum number of flanking base pairs that must exist on
 *     either side of a SplitPoint for it to be valid
 * @param f_split the minimum partition function fraction that must be captured
 *      by the exclusive split points
 * @param part partition of complexes into \f$\Psi^{\text{active}}\f$ and
 *      \f$\Psi^{\text{passive}}\f$
 */
void Design::initialize_decomposition(EnsemblePartition const &part) {
    izip(complexes, [&](auto i, auto &c) {
        if (len(part) > 0 && !part.active(i)) return;
        if (c.is_on_target()) c.structure_decompose();
        else c.probability_decompose(sequence(), models);
    });
}


/**
 * @brief redecompose complexes in the active set: probability- and
 * structure-based if complex has target structure, probability-based if not.
 * Both with potential multiple exclusive split points.
 *
 * @param min_size the minimum number of nucleotides that must be in a child
 *     node for the decomposition to be valid
 * @param min_helix the minimum number of flanking base pairs that must exist on
 *     either side of a SplitPoint for it to be valid
 * @param f_split the minimum partition function fraction that must be captured
 *      by the exclusive split points
 * @param part partition of complexes into \f$\Psi^{\text{active}}\f$ and
 *      \f$\Psi^{\text{passive}}\f$
 */
void Design::redecompose_active(Local const &env, uint depth, EnsemblePartition const &part) {
    auto decomp = [&](auto, auto i, auto) {
        auto &c = at(complexes, i);
        if (len(part) == 0 || part.active(i))
            c.probability_decompose(sequence(), models, depth);
    };

    env.spread(vec<uint>(indices(complexes)), 1, decomp, AffinitySplit());
}


Variant<DomainView, StrandView> find_sequence_element(Design const &design, string const &name) {
    Variant<DomainView, StrandView> ret;
    try {
        ret = design.sequences.get_domain(name);
    } catch (std::out_of_range const &e) {
        try {
            ret = design.sequences.get_strand(name);
        } catch (std::out_of_range const &e) {
            NUPACK_ERROR(name + " is not a strand or domain");
        }
    }
    return ret;
}


template <class F>
uint find_index(string name, Design const &design, F &&f) {
    auto const &collection = f(design);
    auto it = find_if(collection, [&](auto const &x) {return x.name == name;});
    if (it == end_of(collection)) NUPACK_ERROR(name, "Element not found.");
    return it - begin_of(collection);
}


uint find_tube(string name, Design const &design) {
    try {
        return find_index(name, design, [](auto const &x) -> auto const & {return x.tubes;});
    } catch (Error const &e) {
        NUPACK_ERROR("tube not found", name);
    }
}

uint find_complex(string name, Design const &design) {
    try {
        return find_index(name, design, [](auto const &x) -> auto const & {return x.complexes;});
    } catch (Error const &e) {
        NUPACK_ERROR("complex not found", name);
    }
}


}}
