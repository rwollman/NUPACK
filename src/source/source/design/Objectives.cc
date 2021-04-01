#include <nupack/design/Objectives.h>
#include <nupack/design/Weights.h>
#include <nupack/types/Matrix.h>
#include <nupack/state/State.h>

namespace nupack { namespace newdesign {


void TubeObjective::initialize(Design const &design) {
    tube_id = find_tube(tube_name, design);
}


void ComplexObjective::initialize(Design const &design) {
    complex_id = find_complex(complex_name, design);
}


void PatternObjective::initialize(Design const &design) {
    /* use all strands if elements not specified */
    if (len(component_names) == 0) component_names = vec<string>(key_view(design.sequences.strands));

    elements = vmap(component_names, [&](auto const &n) { return find_sequence_element(design, n); });
    NUPACK_REQUIRE(len(component_names), ==, len(elements));

    for (auto const &p : patterns) {
        uint N = len(p);
        grouped_patterns.try_emplace(N, vec<Sequence>());
        at(grouped_patterns, N).emplace_back(p);
    }

    for (auto p: key_view(grouped_patterns))
    for (auto const &el : elements) {
        auto n = fork(el, [](auto const &x) {return len(x);});
        if (n < p) continue;
        normalization += n - p + 1;
    }
}


void SimilarityObjective::initialize(Design const &design) {
    NUPACK_REQUIRE(len(component_names), ==, len(ref_seqs));
    NUPACK_REQUIRE(len(component_names), ==, len(limits));

    /* 0.0 < lower limit < upper limit < 1.0 */
    for (auto const & lim : limits) {
        NUPACK_REQUIRE(lim.first, <, lim.second);
        NUPACK_REQUIRE(0.0, <, lim.first);
        NUPACK_REQUIRE(lim.second, <, 1.0);
    }

    elements = vmap(component_names, [&](auto const &n) { return find_sequence_element(design, n); });
    zip(elements, ref_seqs, [](auto const &el, auto const &rs) {
        auto el_len = fork(el, [](auto const &x) {return len(x);});
        NUPACK_REQUIRE(el_len, ==, len(rs), "reference sequence and design element for SimilarityObjective are different lengths");
    });
    NUPACK_REQUIRE(len(component_names), ==, len(elements));
}


void EnergyEqualizationObjective::initialize(Design const &design) {
    for (auto const &name : domain_names) {
        try {
            domains.emplace_back(std::get<DomainView>(find_sequence_element(design, name)));
        } catch (Error const &e) {
            NUPACK_ERROR("Element is not a domain", name);
        } catch (std::bad_variant_access const &e) {
            NUPACK_ERROR("Element is not a domain", name);
        }
    }
    model = at(design.complexes, 0).target.model;
}


void SSMObjective::initialize(Design const &design) {
    complex_ids = vmap<vec<uint>>(complex_names, [&](auto const &name) {
        return find_complex(name, design);
    });

    add_identicals(design.sequences);
    add_complements(design.sequences);
    process_words(design);
    process_structures(design);
    // BEEP(*this);
}

using custom_csp::IdentConstraint;
using custom_csp::CompConstraint;

template <class T>
void add_binary_relations(DesignSequence const &seqs, NucleotideRelationMap &container) {
    for (auto const &c : seqs.constraints.handler.get_constraints()) {
        auto ptr = maybe_get<T>(c);
        if (ptr == nullptr) continue;
        auto vars = ptr->get_constrained_vars();
        NUPACK_REQUIRE(len(vars), ==, 2, "binary constraint should have two variables");

        auto i = at(vars, 0), j = at(vars, 1);
        container.at(i).emplace(j);
        container.at(j).emplace(i);
    }
}

void SSMObjective::add_identicals(DesignSequence const &seqs) {
    /* add in reflexivity for consistent interface */
    for (auto i : range(len(seqs.nucleotides))) {
        uint ui = i;
        identicals.emplace(ui, std::set<uint>{ui});
    }

    add_binary_relations<IdentConstraint>(seqs, identicals);
}


void SSMObjective::add_complements(DesignSequence const &seqs) {
    for (auto i : range(len(seqs.nucleotides))) { complements.emplace(uint(i), std::set<uint>{}); }

    add_binary_relations<CompConstraint>(seqs, complements);
}


/**
 * @brief creates Ranges that represent windows in the indexing of an individual complex.
 * Converted downstream into indexing in the sequence variable indexing.
 *
 * @param nicks the positions after 3' strand ends, including that of the last strand
 * @return vec<Range<uint>>
 */
vec<Range<uint>> SSMObjective::ranges(Nicks const &nicks) const {
    uint i = 0;
    uint N = back(nicks);
    auto it = begin_of(nicks);
    vec<decltype(range(i, i))> rs;
    while (i < (N - word_size + 1)) {
        if (i + word_size > *it) {
            i = *it;
            ++it;
            continue;
        }

        rs.emplace_back(range(i, i+word_size));
        ++i;
    }
    // BEEP(rs);
    return rs;
}

/**
 * @brief determine if two words are identical based on explicit
 * constraints in this->identicals
 *
 * @param a first word
 * @param b second word
 * @return true the are identical
 * @return false they are not identical
 */
bool SSMObjective::identical(vec<uint> const &a, vec<uint> const &b) const {
    bool matching = true;
    zip(a, b, [&](auto i, auto j) {
        if (!contains(identicals.at(i), j)) matching = false;
    });
    return matching;
}


/**
 * @brief determine if two words are reverse complements of each other based on
 * constraints in this->complements
 *
 * @param a first word
 * @param b second word
 * @return true the are reverse complements
 * @return false they are not reverse complements
 */
bool SSMObjective::complementary(vec<uint> const &a, vec<uint> const &b) const {
    bool inverted = true;
    zip(a, reversed(b), [&](auto i, auto j) {
        if (!contains(complements.at(i), j)) inverted = false;
    });
    return inverted;
}


void SSMObjective::process_words(Design const &design) {
    for (auto const &c : indirect_view(complex_ids,
        [&](auto i) -> Complex const & { return at(design.complexes, i); })) {
        auto indices = c.to_indices();
        auto rs = ranges(c.target.structure.nicks);
        for (auto r : rs)
            words.emplace_back(vmap(r, [&](auto i) {return at(indices, i);}));
    }

    std::set<vec<uint>> counter(begin_of(words), end_of(words));
    normalization = len(counter);
}


/**
 * @brief fills in complement_restricted by looking for windows whose complement
 * is not a perfect duplex in the target structure (unpaired nucleotides or
 * nicks in the paired sequence).
 *
 * @param design
 */
void SSMObjective::process_structures(Design const &design) {
    auto condition = [&] (auto const &positions, auto const &struc) {
        NUPACK_REQUIRE(len(positions), ==, word_size);

        /* reverse order paired positions */
        vec<uint> paired {reversed(subview(struc, positions))};

        auto it = begin_of(paired);
        auto prev = *it;
        ++it;

        bool contiguous = true;
        while (it != end_of(paired)) {
            if (*it != prev + 1) {
                contiguous = false;
                break;
            }
            prev = *it;
            ++it;
        }

        /* AND no nicks */
        bool no_nicks = none_of(view(begin_of(paired) + 1, end_of(paired)), [&](auto i) {
            return contains(struc.nicks, i);
        });

        // BEEP(positions, paired, contiguous, no_nicks);
        return contiguous && no_nicks;
    };

    for (auto const &c : indirect_view(complex_ids,
            [&](auto i) -> Complex const & { return at(design.complexes, i); })) {
        // BEEP(c.name);
        auto indices = c.to_indices();
        // BEEP(indices);
        auto struc = c.target.structure;
        // BEEP(struc);

        auto lengths = vmap(c.strands, len);
        auto num_critons = sum(lengths, [&](int l) {
            return max(0, l - int(word_size) + 1);
        });
        auto rs = ranges(struc.nicks);
        NUPACK_REQUIRE(len(rs), ==, num_critons);

        for (auto positions : rs) {
            if (!condition(positions, struc)) {
                complement_restricted.emplace(vec<uint>(subview(indices, positions)));
            }
        }
    }
}




Defect MultitubeObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const {
    return design.normalized_defect(env, depth, part, {}, weights, obs);
}


Defect TubeObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &, EngineObserver &obs) const {
    auto const &tube = at(design.tubes, tube_id);
    auto log_pfuncs = design.log_pfuncs(env, depth, part, {}, obs);
    auto complex_defects = design.complex_defects(env, depth, part, {}, obs);
    // auto tube_weights = bool(weights) ? weights.per_tube.at(tube_id) : ComplexWeights();
    return tube.normalized_defect(log_pfuncs, complex_defects, part);
}


Defect ComplexObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &, Weights const &, EngineObserver &obs) const {
    auto const &comp = at(design.complexes, complex_id);
    auto defect = comp.defect(env, design.models, design.sequence(), depth, {}, obs);
    uint N = len(comp);
    for (auto &i : defect.contributions) i.second /= N;
    return defect;
}


Defect PatternObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &) const {
    vec<real> defs(len(design.sequence()), 0.0);

    for (auto const &el : elements) {
        auto seq = fork(el, [&](auto const &x) { return x.to_sequence(design.sequence()); });
        auto indices = fork(el, [&](auto const &x) { return x.to_indices(); });

        for (auto const &[n, ps] : grouped_patterns) {
            real per_nuc = 1.0 / real(n);

            for (auto start : range(len(seq) - n + 1)) {
                auto sp = span(start, start+n);
                auto sub_seq = subview(seq, sp);

                bool matched = any_of(ps, [&](auto const &p) {
                    bool all_matched = true;
                    zip(sub_seq, p, [&](auto a, auto b) {
                        if (!is_base_specialization(b, a)) all_matched = false;
                    });
                    return all_matched;
                });

                if (matched) for (auto i : subview(indices, sp)) at(defs, i) += per_nuc;
            }
        }
    }

    return Defect(defs, normalization);
}


/**
 * @brief determine whether each real sequence is above or below the matching
 * limits and penalize accordingly
 *
 * @param env compute environment (single or multithreaded)
 * @param design the active design
 * @param depth ignored
 * @param part ignored
 * @return Defect average fraction of nucleotides not in similarity ranges out
 * of maximum number of nucleotides that could be incorrect.
 */
Defect SimilarityObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &) const {
    real normalization = 0;
    zip(ref_seqs, limits, [&](auto const &ref, auto const &lim) {
        /* maximum number of nucleotides that can be incorrectly matched per sequence */
        normalization += len(ref) * std::max(lim.first, 1.0 - lim.second);
    });

    vec<real> mapped_defects(len(design.sequence()), 0.0);
    zip(elements, ref_seqs, limits, [&](auto const &el, auto const &ref, auto const &lim) {
        auto seq = fork(el, [&](auto const &x) { return x.to_sequence(design.sequence()); });
        vec<uint> matches(len(seq), 0);
        zip(seq, ref, matches, [](auto const &s, auto const &r, auto &m) {
            m = is_base_specialization(r, s);
        });

        real num_matches = sum(matches);
        real N = len(seq);
        real frac = num_matches / N;
        auto indices = fork(el, [&](auto const &x) { return x.to_indices(); });
        if (frac < lim.first) { // need more matches
            auto per_nuc = (lim.first - frac) / frac;
            zip(indices, matches, [&](auto i, auto m) { if (!m) at(mapped_defects, i) += per_nuc; });
        } else if (frac > lim.second) { // need fewer matches
            auto per_nuc = (frac - lim.second) / frac;
            zip(indices, matches, [&](auto i, auto m) { if (m) at(mapped_defects, i) += per_nuc; });
        }
    });

    defect_vec defs;
    izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d / normalization);});

    return {defs};
}


/**
 * @brief the sum of the absolute distances from either the median energy or the
 * reference energy dividied by (this sum + |reference energy/median energy|)
 *
 * @param env compute environment (single or multithreaded)
 * @param design the active design
 * @param depth ignored
 * @param part ignored
 * @return Defect measure of how close to reference energy / median energy
 */
Defect EnergyEqualizationObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &) const {
    /* convert to duplex of length l */
    auto pl = [](uint l) -> PairList {
        std::ostringstream ss;
        ss << "(" << l << "+)" << l;
        return Structure(ss.str());
    };

    /* compute energies for each of the domains + their complements using
    structure_energy */
    vec<real> v_energies = vmap<vec<real>>(domains, [&](auto const &domain) {
        auto seq = domain.to_sequence(design.sequence());
        vec<Sequence> seqs = {seq, reverse_complement(seq)};
        return structure_energy(seqs, pl(len(seq)), model);
    });
    real_col energies(v_energies);

    /* use arma to calculate median and absolute diffs from it */
    real goal = bool(ref_energy) ? *ref_energy : median(energies);
    real_col diffs = energies - goal;

    /* TODO: include in design parameters somehow */
    real scale = 10; // kcal/mol
    diffs = 1 - arma::exp((-1 / scale) * arma::abs(diffs));
    real denom = len(domains);

    /* used x / (x + a) function to find defects per domain */
    real_col col_per_domain = diffs / denom;
    vec<real> per_domain = {col_per_domain.begin(), col_per_domain.end()};

    /* equally split defect per domain into component nucleotides */
    vec<real> mapped_defects(len(design.sequence()), 0.0);
    zip(domains, per_domain, [&] (auto const &d, auto p) {
        p = p / len(d);
        for (auto i : d.to_indices()) at(mapped_defects, i) += p;
    });
    defect_vec defs;
    izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d);});

    return {defs};
}


/**
 * @brief for each sequence of length word_size actually appearing as a
 * contiguous substrand in the set of complexes, compute the number of unrelated
 * sequences that use this word and penalize each accordingly
 *
 * @param env ignored
 * @param design the active design
 * @param depth ignored
 * @param part ignored
 * @return Defect number of spurious appearances of the same sequence
 */
Defect SSMObjective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &) const {
    std::map<Sequence, Index_Map> critons;

    auto create_or_expand = [&](auto const &seq, auto const &inds) {
        critons.try_emplace(seq, Index_Map());
        critons.at(seq).add(inds);
    };

    auto sequences = design.sequence();
    for (auto const &word : words) {
        Sequence seq = Sequence(vmap(word, [&](auto i) { return at(sequences, i); }));
        create_or_expand(seq, word);
        if (complement_restricted.count(word) && !is_palindromic(seq)) {
            create_or_expand(reverse_complement(seq), word);
        }
    }

    real accum = 0;
    vec<real> mapped_defects(len(design.sequence()), 0.0);
    // BEEP(design.sequences.json_domains());
    for (auto &c : critons) {
        auto const &seq = c.first;
        auto &ind_map = c.second;

        /* find number of distinct conflicting underlying variables with same sequence */
        ind_map.resolve_groups(*this);

        /* find number of distinct conflicting underlying variables with same sequence */
        if (is_palindromic(seq)) ind_map.num_violations += 1;

        real increment = ind_map.assign_blame(mapped_defects);
        // if (increment > 0) BEEP(seq, ind_map);
        accum += increment;
    }
    accum /= normalization;

    defect_vec defs;
    izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d / normalization);});
    Defect ret {defs};

    NUPACK_REQUIRE(ret.total(), ==, about(accum));
    return ret;
}


void Index_Map::resolve_groups(SSMObjective const &obj) {
    uint group_num = 0;

    for (auto it = begin_of(used); it != end_of(used); ++it) {
        if (it->assigned()) continue;

        it->group = group_num;

        for (auto it2 = it + 1; it2 != end_of(used); ++it2) {
            if (it2->assigned()) continue;

            auto const & a = it->indices;
            auto const & b = it2->indices;

            if (obj.identical(a, b) || obj.complementary(a, b))
                it2->group = group_num;
        }

        ++group_num;
    }

    num_violations = group_num - 1;
}


real Index_Map::assign_blame(vec<real> &defects) const {
    real per_nuc = real(num_violations) / real(sum(used, [](auto const &u) { return len(u.indices);}));
    for (auto const &u : used) {
        for (auto i : u.indices) {
            at(defects, i) += per_nuc;
        }
    }

    return real(num_violations);
}


Optional<Defect> MultitubeObjective::reevaluate(Local const &env,
        Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const {
    return evaluate(env, design, depth, part, weights, obs);
}


Optional<Defect> TubeObjective::reevaluate(Local const &env,
        Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const {
    return evaluate(env, design, depth, part, weights, obs);
}


Optional<Defect> ComplexObjective::reevaluate(Local const &env,
        Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const {
    return evaluate(env, design, depth, part, weights, obs);
}







void Objective::initialize(Design const &design) {
    fork(variant, [&](auto &x) { x.initialize(design); });
}

Defect Objective::evaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const {
    return fork(variant, [&](auto const &x) {
        return x.evaluate(env, design, depth, part, weights, obs);
    });
}

Optional<Defect> Objective::reevaluate(Local const &env, Design const &design, uint depth, EnsemblePartition const &part, Weights const &weights, EngineObserver &obs) const {
    return fork(variant, [&](auto const &x) {
        return x.reevaluate(env, design, depth, part, weights, obs);
    });
}




}}
