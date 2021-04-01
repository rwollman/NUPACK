#include <nupack/design/Specification.h>
#include <nupack/design/OutputResult.h>
// #include <spdlog/spdlog.h>
// #include <spdlog/sinks/basic_file_sink.h>

namespace nupack { namespace newdesign {

// using custom_csp::CompConstraint;
// using custom_csp::CompConstraint;
// using custom_csp::IdentConstraint;
// using custom_csp::PatternConstraint;
// using custom_csp::WordConstraint;
// using custom_csp::MatchConstraint;
// using custom_csp::NUPACK_CS_STRONG;
// using custom_csp::trinary;


Specification::operator Designer() const {
    DesignSequence seqs;
    seqs.wobble_mutations = wobble_mutations;

    /** Sequence-level operations */
    /* add domains */
    for_each(domains, [&](auto const &x) {seqs.add_domain(x);});
    /* add strands */
    for_each(strands, [&](auto const &x) {seqs.add_strand(x);});

    seqs.make_sequence();

    /** constraint-level operations */
    /** match constraints */
    for_each(constraints.match, [&](auto const &c) {
        auto vars = c.get_variables(seqs);
        zip(vars.first, vars.second, [&](int i, int j) {
            // seqs.constraints.add_constraint(IdentConstraint(i, j));
            seqs.constraints.match_constraint(i, j);
        });
    });
    /** complementarity constraints */
    for_each(constraints.complementarity, [&](auto const &c) {
        auto vars = c.get_variables(seqs);
        zip(vars.first, reversed(vars.second), [&](int i, int j) {
            // seqs.constraints.add_constraint(CompConstraint(i, j, NUPACK_CS_STRONG));
            seqs.constraints.complementarity_constraint(i, j, wobble_mutations);
        });
    });
    /** pattern constraints */
    for_each(constraints.pattern,     [&](auto const &c) {c.add_constraint(seqs);});
    /** diversity constraints */
    for_each(constraints.diversity,   [&](auto const &c) {c.add_constraint(seqs);});
    /** word (library and window) constraints */
    for_each(constraints.word,        [&](auto const &c) {c.add_constraint(seqs);});
    /** similarity constraints */
    for_each(constraints.similarity,  [&](auto const &c) {c.add_constraint(seqs);});

    /** Design level operations */
    Design design(std::move(seqs));
    /* add complexes */
    auto comp_name = [](auto const &c) {
        if (c.name != "") return c.name;
        std::stringstream s; s << at(c.strands, 0);
        for (auto i : range(1, len(c.strands))) s << "-" << c.strands[i];
        return s.str();
    };

    DecompositionParameters params {parameters.H_split, parameters.N_split, parameters.f_split, parameters.f_sparse, parameters.dG_clamp};

    for_each(complexes, [&](auto const &x) {design.add_complex(x.strands, model, comp_name(x), x.structure, params);});

    design.add_structure_complementarity();

    auto comp = [&](auto const &x) {return complex_index(x);};

    /* add tubes */
    for_each(tubes, [&](auto const &x) {
        auto indices = vec<uint>(indirect_view(key_view(x.targets), comp));
        auto concs = vec<real>(item_view(x.targets));
        design.add_tube(indices, concs, x.name);
    });

    Designer ret(std::move(design), std::move(objectives), std::move(weights), std::move(parameters));
    return ret;
}


/**
 * @brief returns a pair of the variables for the left and right concatenated sequences
 *
 * @param seqs the DesignSequence to pull the variable mapping from
 * @return a pair with the left indices first and the right indices second
 */
std::pair<vec<int>, vec<int>> DualListSpec::get_variables(DesignSequence const &seqs) const {
    return {extract_variables(left, seqs), extract_variables(right, seqs)};
}


/**
 * @brief converts the specification into a PatternConstraint and adds it to the
 *     list of constraints in seqs
 * @details if name is the empty string, pattern is prevented in every strand
 *     in seqs. Otherwise, name is looked up as either a strand or a domain
 *     and the pattern is prevented in that one element.
 *
 * @param seqs the DesignSequence to which to add the constraint implied by
 *     the specification
 */
void PatternSpec::add_constraint(DesignSequence &seqs) const {
    /* handle global case */
    // auto poss = seqs.constraints.get_possible_nucleotides();
    if (domains.empty()) {
        auto vars = vmap(seqs.strands, [&](auto const &s) {
            return vmap(s.second.to_indices(), [](auto i) {return int(i);});
        });
        for_each(vars, [&](auto const &v) {
            // seqs.constraints.add_constraint(PatternConstraint(v, pattern, poss));
            seqs.constraints.pattern_constraint(v, Sequence(pattern));
        });
    } else {
        // seqs.constraints.add_constraint(PatternConstraint(extract_element(name, seqs), pattern, poss));
        seqs.constraints.pattern_constraint(extract_variables(domains, seqs), Sequence(pattern));
    }
}


void DiversitySpec::add_constraint(DesignSequence &seqs) const {
    if (domains.empty()) {
        auto vars = vmap(seqs.strands, [&](auto const &s) {
            return vmap(s.second.to_indices(), [](auto i) {return int(i);});
        });
        for_each(vars, [&](auto const &v) {
            seqs.constraints.diversity_constraint(v, word_length, min_nucleotide_types);
        });
    } else {
        seqs.constraints.diversity_constraint(extract_variables(domains, seqs), word_length, min_nucleotide_types);
    }
}

/**
 * @brief converts the specification into a WordConstraint and adds it to the
 *     list of constraints in seqs
 * @details Adding a WordConstraint must be preceded by adding a variable to
 *     the ConstraintHandler to represent which of the enumerated words is
 *     still accessible during propagation. Both adding the variable and then
 *     adding the constraint are handled in this function.
 *
 * @param seqs the DesignSequence to which to add the constraint implied by
 *     the specification
 */
void WordSpec::add_constraint(DesignSequence &seqs) const {
    auto vars = extract_variables(domains, seqs);
    int i = 0;
    for_each(comparisons, [&](auto const &c) {
        auto length = len(at(c, 0));
        // auto supp_var = seqs.constraints.add_variable(vec<trinary>(len(c), true));
        auto cur = vec<int>(subview(vars, span{i, i+length}));

        // seqs.constraints.add_constraint(WordConstraint(cur, c, supp_var));
        auto ref_seqs = vmap<vec<Sequence>>(c, [](auto const &x) {return Sequence(x);});
        seqs.constraints.word_constraint(cur, ref_seqs);
        i += length;
    });
}


/**
 * @brief converts the specification into a MatchConstraint and adds it to the
 *     list of constraints in seqs
 *
 * @param seqs the DesignSequence to which to add the constraint implied by
 *     the specification
 */
void SimilaritySpec::add_constraint(DesignSequence &seqs) const {
    // seqs.constraints.add_constraint(MatchConstraint(extract_element(name, seqs), reference, {range.first}, {range.second}));
    seqs.constraints.similarity_constraint(extract_variables(domains, seqs), Sequence(reference), range);
}


/**
 * @brief return a concatenation of the indices associated with the elements
 *     in names if they are valid elements of seqs
 *
 * @param names the names of strands or domains
 * @param seqs the DesignSequence from which to pull the elements
 *
 * @return the concatenation of the individual indices from each of the
 *     elements named in names
 */
vec<int> extract_variables(vec<string> const &names, DesignSequence const &seqs) {
    auto temp = vmap(names, [&](auto const &name) {
        return extract_element(name, seqs);
    });
    return join(temp);
}


/**
 * @brief return the variable indices corresponding to the named element
 * @details attempt to find either a domain or strand with the given name,
 *     convert the DomainSpec or StrandSpec into its vector of indices and
 *     return. If no element is found with the given name, an exception is
 *     thrown.
 *
 * @param name the name of a strand or domain
 * @param seqs the DesignSequence from which to pull the DomainSpec or
 *     StrandSpec
 * @throws nupack::Error if name is neither a strand nor a domain in seqs
 *
 * @return a vector of the indices corresponding to the named strand or domain
 */
vec<int> extract_element(string name, DesignSequence const &seqs) {
    vec<uint> vars;
    try {
        auto domain = seqs.get_domain(name);
        vars = domain.to_indices();
    } catch (std::out_of_range const &e) {
        try {
            auto strand = seqs.get_strand(name);
            vars = strand.to_indices();
        } catch (std::out_of_range const &e) {
            NUPACK_ERROR(name + " is not a strand or domain");
        }
    }
    return vmap(vars, [](auto i) {return int(i);});
}




vec<uint> Specification::ensure_compatibility(Specification const &spec, SingleResult const &res) {
    /* domains */
    NUPACK_REQUIRE(len(spec.domains), ==, len(res.domains), "mismatched number of domains");
    for (auto const &d : spec.domains) {
        auto it = res.domains.find(d.name);
        NUPACK_ASSERT(it != end_of(res.domains), "result missing domain", d);
        NUPACK_REQUIRE(len(d.allowed_bases), ==, len(it->second), "different domain lengths", d, *it);
    }

    /* strands */
    NUPACK_REQUIRE(len(spec.strands), ==, len(res.strands), "mismatched number of strands");
    for (auto const &s : spec.strands) {
        auto it = res.strands.find(s.name);
        NUPACK_ASSERT(it != end_of(res.strands), "result missing strand", s);
        auto sslength = sum(s.domain_names, [&](auto const &name) {
            return len(find_if(spec.domains, [&](auto const &d) {return d.name == name;})->allowed_bases);
        });
        NUPACK_REQUIRE(sslength, ==, len(it->second), "different strand lengths", s, *it);
    }

    auto get_strand_seq = [&](auto s) {return res.strands.at(s);};

    /* complexes */
    auto mapping = vmap<vec<uint>>(spec.complexes, [&](auto const &c) {
        auto it = find_if(res.complexes, [&] (auto const &c2) {
            return lowest_rotation(vmap<StrandList>(c.strands, get_strand_seq))
                    == lowest_rotation(c2.sequence.strands());
        });
        NUPACK_REQUIRE(it, !=, end_of(res.complexes),
                "strands in checkpoint complex do not match expected sequence given strand names in specification", c);
        return it - begin_of(res.complexes);
    });


    // compare TubeComplex:name -(find)-> Complexes
    /* tubes */
    NUPACK_REQUIRE(len(spec.tubes), ==, len(res.tubes), "mismatched number of tubes");

    auto comp_name = [](auto const &a, auto const &b) {return a.name < b.name;};

    zip(sorted(spec.tubes, comp_name), sorted(res.tubes, comp_name), [&] (auto const &s, auto const &r) {
        auto temp = vmap(r.complexes, [&](auto const &c) {
            auto const &named_complex = *find_if(res.complexes, [&](auto const &rc) {return c.name == rc.name;});
            return lowest_rotation(named_complex.sequence.strands());
        });

        std::set<StrandList> res_seqs(begin_of(temp), end_of(temp));

        for (auto const &c : s.targets) {
            auto tubename = s.name;
            auto strand_names = c.first;
            auto const &comp = at(spec.complexes, spec.complex_index(strand_names));
            auto strand_seqs = vmap<StrandList>(comp.strands, get_strand_seq);
            NUPACK_ASSERT(res_seqs.find(lowest_rotation(strand_seqs)) != end_of(res_seqs),
                    "complex in given tube in specification not found in matching tube in result",
                    tubename, comp, strand_seqs);
        }
    });

    return mapping;
}

}}
