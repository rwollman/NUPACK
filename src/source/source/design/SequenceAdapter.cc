#include <nupack/design/SequenceAdapter.h>
#include <nupack/reflect/Serialize.h>

namespace nupack { namespace newdesign {




/**
 * @brief return a view of the substrand between position begin and end, INCLUSIVE
 */
StrandView StrandView::slice(uint begin, uint end) const {
    auto length = size();
    if (begin > end) throw std::runtime_error("invalid slice. begin must be <= end.");
    if (begin >= length) throw std::out_of_range("begin is not in range.");
    if (end >= length) throw std::out_of_range("end is not in range.");

    auto lengths = prefixes<vec<uint>>(false, indirect_view(domains, len));
    auto prev = prefixes<vec<uint>>(true, indirect_view(domains, len));

    auto first = upper_bound(lengths, begin) - begin_of(lengths);
    auto last = upper_bound(lengths, end) - begin_of(lengths);
    auto b = begin - prev[first]; // begin relative to first domain
    auto e = end - prev[last] + 1; // end relative to last domain


    // single domain edge case
    if (first == last) return decltype(domains) {DomainView(domains[first].start() + b, domains[first].start() + e)};

    decltype(domains) new_domains;
    new_domains.emplace_back(domains[first].start() + b, domains[first].stop());
    for (auto i : range(first+1, last)) new_domains.emplace_back(domains[i]);
    new_domains.emplace_back(domains[last].start(), domains[last].start() + e);
    return new_domains;
}


/**
 * @brief add any missing domain complements. this function is idempotent once
 *     non-complement domains are constant.
 */
void DesignSequence::add_domain_complements() {
    auto cur_domains = domain_specs;
    for (auto const &dom : cur_domains) {
        if (back(dom.name) == '*') continue;

        auto complement_name = dom.name + "*";
        if (none_of(cur_domains, [&](auto const &c) {return complement_name == c.name;})) {
            // can eventually add funcitonality to restrict domains for
            // possible performance improvement, but it shouldn't matter for
            // correctness of applying constraints.
            add_domain(complement_name, string(len(dom), 'N'));
        }
    }
}


/**
 * @brief checks if any nucleotide in the domains is free to vary. Used for
 * determining if a design can be evaluated uniquely, i.e. prerequisite for
 * multitubedefect behavior
 *
 * @return std::pair<bool, string> the bool is true if all nucleotides are
 * fixed; if the bool is false, the string contains the name of one domain which
 * is not fixed, otherwise it is the empty string
 */
std::pair<bool, string> DesignSequence::all_nucleotides_fixed() {
    add_domain_complements();

    for (auto const &dom: domain_specs) {
        string name = dom.name;
        if (back(name) == '*') // complement domain
            continue;
        string complement_name = name + "*";

        Sequence domain(dom.allowed_bases);
        auto it = find_if(domain_specs, [&](auto const &d) {return d.name == complement_name;});
        if (it == end_of(domain_specs)) NUPACK_ERROR("complement was not added correctly", complement_name);
        Sequence comp_domain(it->allowed_bases);

        if (!all_of(domain, is_canonical) && !all_of(comp_domain, is_canonical))
            return {false, name};
    }

    return {true, ""};
}


/**
 * @brief convert current set of domain and strand specs into underlying
 *     sequence and views on the sequence for the domains and strands.
 */
void DesignSequence::make_sequence() {
    domains.clear(); strands.clear();
    // constraints = ConstraintHandler();

    add_domain_complements();

    nucleotides = join(indirect_view(domain_specs, [](auto const & el) {
        return Sequence(el.allowed_bases);
    }));

    int last = 0;
    for (auto const & d: domain_specs) {
        auto beg = last; last += len(d); auto end = last;
        domains.emplace(d.name, DomainView(beg, end));
    }

    for (auto const & s: strand_specs) {
        strands.emplace(s.name, vmap<vec<DomainView>>(s.domain_names, [&](auto const &dn) {
            return get_domain(dn);
        }));
    }

    // add variables to ConstraintHandler
    // for_each(nucleotides, [&](auto const &n) {constraints.add_nucleotide_variable(n);});
    constraints = Constraints(nucleotides);
    add_complementarity_constraints();
}


/**
 * @brief initialize nucleotides by initializing constraints and
 *     converting variables to nucleotides
 * @details nucleotides are resized and blanked out as the addition of any
 *     WordConstraints increases the number of variables in the
 *     ConstraintHandler. This resize only needs to be done once at
 *     initialization assuming no further constraints are added after this
 *     point.
 */
void DesignSequence::initialize_sequence() {
    // auto vars = constraints.init_random();
    // nucleotides = Sequence(len(vars), 'N');
    // transform(vars, nucleotides, [](auto i) {return i < 16 ? Base::from_index(i) : Base('_');});
    auto result = constraints.initial_sequence();
    if (result) {
        nucleotides = result.value();
    } else {
        NUPACK_ERROR("unable to find sequence satisfying all constraints");
    }
}


/**
 * @brief Add complementarity constraints between domain x and x* for all domains x
 */
void DesignSequence::add_complementarity_constraints() {
    for (auto const &dom : domains) {
        auto name = dom.first;
        if (back(name) != '*') {
            auto comp_name = name + "*";
            auto it = domains.find(comp_name);
            if (it != domains.end())
                zip(dom.second.indices, ~(it->second.indices), [&](auto i, auto j) {
                    // constraints.add_constraint(CompConstraint(i, j, NUPACK_CS_STRONG));
                    constraints.complementarity_constraint(i, j, wobble_mutations); // only Watson Crick for now
                });
        }
    }
}

/**
 * @brief set a given domain to the sequence
 *
 * @param name the domain name
 * @param in the sequence to change the domain to
 */
void DesignSequence::set_domain(string const &name, Sequence const &in) {
    auto domain = get_domain(name);
    auto dom_spec = Sequence(find_if(domain_specs, [&](auto const &d) {return d.name == name;})->allowed_bases);

    NUPACK_REQUIRE(len(domain), ==, len(in), "input sequence does not match domain length", name);
    NUPACK_ASSERT(all_determined(in), "Cannot assign degenerate base codes to domain");
    NUPACK_ASSERT(is_sequence_specialization(dom_spec, in), "in nucleotides are not compatible with domain", dom_spec, in);

    izip(domain.indices, [&](auto i, auto j) {nucleotides[j] = in[i];});
    NUPACK_REQUIRE(domain.to_sequence(nucleotides), ==, in);
}


/**
 * @brief mutate sequence at the given positions but maintain constraint satisfaction.
 *
 * @param vars positions to change
 */
bool DesignSequence::mutate_sequence(vec<uint> const &vars) {
    vec<int> muts{view(vars)};
    // vec<int> current_variables{view(nucleotides)};
    auto old_sequence = nucleotides;
    // transform(constraints.make_mutation(muts, current_variables), nucleotides, [](auto i) {
    //     /** handling WordConstraint non-nucleotide variables with values exceeding Base
    //       * representation
    //       */
    //     return i < len(Base::names) ? Base::from_index(i) : Base('_');
    // });

    if (auto result = constraints.make_mutation(nucleotides, muts)) {
        nucleotides = std::move(*result);
        return true;
    } else {
        return false;
    }
}


string DesignSequence::json_domains(Sequence s) const {
    if (len(s) == 0) s = nucleotides;

    std::map<string, string> temp;
    for (auto const &domain : domains) temp.emplace(domain.first, domain.second.to_sequence(s));

    std::ostringstream ss;
    ss << json(temp);
    return ss.str();
}


}}
