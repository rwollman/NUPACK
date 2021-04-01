#include <nupack/design/Weights.h>
#include <nupack/design/Design.h>

namespace nupack { namespace newdesign {

ReversedComplex::ReversedComplex(Design const &design, uint index) {
    reverse_map(design, index);
}


void ReversedComplex::reverse_map(Design const &design, uint index) {
    auto const &complex = at(design.complexes, index);
    _N = len(complex);

    auto reverser = [] (auto const &x) { return std::make_pair(x.second, x.first); };

    std::map<StrandView, string> strand_map(indirect_view(design.sequences.strands, reverser));
    std::map<DomainView, string> domain_map(indirect_view(design.sequences.domains, reverser));

    /* map the range of nucleotides corresponding to a strand to that strand name */
    uint i = 0;
    for (auto const &s : complex.strands) {
        auto it = strand_map.find(s);
        if (it == strand_map.end()) NUPACK_BUG("discovered new strand, and that shouldn't be possible", s);
        uint cur_length = len(it->first);
        _strands.emplace(std::make_pair(i, i + cur_length), it->second);
        i += cur_length;
    }
    NUPACK_REQUIRE(i, ==, _N, "sum of strands does not equal complex length");

    /* convert the ordered list of strands to an ordered list of domains that
    occur in the complex */
    vec<DomainView> complex_domains;
    for (auto const &s : complex.strands) cat(complex_domains, s.domains);

    /* map the range of nucleotides corresponding to a domain to that domain name */
    i = 0;
    for (auto const &d : complex_domains) {
        auto it = domain_map.find(d);
        if (it == domain_map.end()) NUPACK_BUG("discovered new domain, and that shouldn't be possible", d);
        uint cur_length = len(it->first);
        _domains.emplace(std::make_pair(i, i + cur_length), it->second);
        i += cur_length;
    }
    NUPACK_REQUIRE(i, ==, _N, "sum of domains does not equal complex length");
}


vec<string> ReversedComplex::domains() const {
    vec<string> ret(_N);
    for (auto const &el : _domains) {
        auto r = range(el.first.first, el.first.second);
        for (auto i : r) at(ret, i) = el.second;
    }
    return ret;
}


vec<string> ReversedComplex::strands() const {
    vec<string> ret(_N);
    for (auto const &el : _strands) {
        auto r = range(el.first.first, el.first.second);
        for (auto i : r) at(ret, i) = el.second;
    }
    return ret;
}


void Weights::resolve_weights(Design const &design) {
    auto it = partition(specifications, [](auto const &spec) {
        return bool(spec.tube); // has tube component
    });
    auto tube_specific = view(begin_of(specifications), it);
    auto complex_specific = view(it, end_of(specifications));

    /* build up list of on-targets for design and */
    vec<uint> on_targets;
    izip(design.complexes, [&](auto i, auto const &c) {
        if (c.is_on_target()) {
            per_complex.emplace(i, vec<real>(len(c), 1.0));
            on_targets.emplace_back(i);
        }
    });

    /* create ReversedComplexes for each of the on_targets */
    make_reversed_complexes(design, on_targets);

    /* process complex level weights */
    for (auto const &w : complex_specific) {
        auto complexes = bool(w.complex) ? vec<uint>{find_complex(w.complex.value(), design)} : on_targets;

        for (auto i : complexes) resolve_single_complex(per_complex, i, w);
    }

    /* copy complex level weights into appropriate tube weights */
    izip(design.tubes, [&](auto i, auto const &t) {
        ComplexWeights temp;
        for (auto const &targ : t.targets) {
            if (targ.is_on_target()) {
                uint j = targ.complex_index;
                temp.emplace(j, at(per_complex, j));
            }
        }
        per_tube.emplace(i, std::move(temp));
    });

    /* process tube specific weights */
    for (auto const &w : tube_specific) {
        uint tube_index = find_tube(w.tube.value(), design);
        auto &current_tube_complexes = at(per_tube, tube_index);

        vec<uint> tube_on_targets(key_view(current_tube_complexes));
        auto complexes = bool(w.complex) ? vec<uint>{find_complex(w.complex.value(), design)} : tube_on_targets;

        /* complain if on-target isn't in tube */
        if (bool(w.complex) && !contains(tube_on_targets, at(complexes, 0))) {
            NUPACK_ERROR("Tube does not contain this on-target", w.tube.value(), w.complex.value());
        }

        for (auto i : complexes) resolve_single_complex(current_tube_complexes, i, w);
    }
}


/**
 * @brief process alternative possibilities for the weight and multiply out the
 * weight to the correct nucleotide positions
 *
 * @param cws map of complex indices to weights for that complex
 * @param index the index of the complex
 * @param w the weight to apply
 */
void Weights::resolve_single_complex(ComplexWeights &cws, uint index, Weight const &w) {
    auto it = cws.find(index);
    if (it == end_of(cws)) NUPACK_ERROR("Not an on-target");
    auto &comp_weights = it->second;
    auto strands = at(reversed_complexes, index).strands();
    auto domains = at(reversed_complexes, index).domains();

    /* only alternative conditions as complex and tube are fixed at this level */
    if (!bool(w.strand) && !bool(w.domain)) {
        for (auto &n : comp_weights) n *= w.weight;
    } else if (!bool(w.strand) && bool(w.domain)) {
        zip(comp_weights, domains, [&](auto &n, auto const &d) {
            if (d == w.domain) n *= w.weight;
        });
    } else if (bool(w.strand) && !bool(w.domain)) {
        zip(comp_weights, strands, [&](auto &n, auto const &s) {
            if (s == w.strand) n *= w.weight;
        });
    } else if (bool(w.strand) && bool(w.domain)) {
        zip(comp_weights, strands, domains, [&](auto &n, auto const &s, auto const &d) {
            if (s == w.strand && d == w.domain) n *= w.weight;
        });
    }
}


void Weights::make_reversed_complexes(Design const &design, vec<uint> const &on_targets) {
    for (auto i: on_targets) reversed_complexes.try_emplace(i, ReversedComplex(design, i));
}


Weight::Weight(Optional<string> t, Optional<string> c, Optional<string> s, Optional<string> d, real w) :
        tube(t), complex(c), strand(s), domain(d), weight(w) {
    if (!(bool(tube) || bool(complex) || bool(strand) || bool(domain))) {
        NUPACK_ERROR("weight must have a scope: tube, complex, strand, and/or domain");
    }
}




}}
