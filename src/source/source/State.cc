/** \file State.cc
 * @brief Defines any non-templated non-inline functions for State classes
 */
#include <nupack/state/State.h>
#include <nupack/model/Model.h>

namespace nupack {

StateBase::StateBase(std::shared_ptr<System const> y, PairList p) : sys(std::move(y)), pairs(std::move(p)) {
    if (!sys) NUPACK_ERROR("Empty System pointer");
    auto const n = len(sys->total_sequence);
    if (pairs.empty()) pairs = PairList(n);
    else if (len(pairs) == sys->n_bases())
        pairs = pairs.with_null_bases(small_vec<iseq>{
            indirect_view(sys->strands, [](auto const &s) {return len(s) - 2u;})});
    else if (len(pairs) != n)
        NUPACK_ERROR("pair list length doesn't match sequence length", len(pairs), n);
    // for_pseudoknots(pairs, [](auto i, auto j, auto k, auto l) { // doesn't work if strands out of order.
    //     NUPACK_ERROR("structure contains at least one pseudoknot", i, j, k, l);
    // });
}

/******************************************************************************************/

PairList StateBase::aligned_pairs(StateBase const &w) const {
    // Make a list of the {index in *this, position in *this, length of strand} for w's strands
    auto data = reserved<vec<std::array<iseq, 3>>>(len(sys->strands));
    for (auto const &s : w.sys->strands) for (auto i : indices<iseq>(sys->strands))
        if (sys->strands[i] == s && !contains(key_view(data), i))
            data.push_back({i, iseq(sum(view(sys->strands, 0, i), len)), iseq(len(s))});

    // Make a base to base map of index in w to index in *this
    auto map = reserved<vec<iseq>>(len(*sys));
    for (auto const &t : data)
        for (auto i : range(t[1], t[1] + t[2])) map.push_back(i);

    // Copy the data into the right places in a new pair vector
    auto pairs = w.pairs;
    iseq n = 0;
    for (auto const &t : data)
        for (auto i : range(t[1], t[1] + t[2])) pairs[i] = map[w.pairs[n++]];
    return pairs;
}

/******************************************************************************************/

std::size_t StateBase::symmetry() const {
    std::size_t const sym1 = rotational_symmetry(sys->strands);
    if (sym1 == 1) return 1;

    auto n = sys->nicks;
    std::adjacent_difference(n.begin(), n.end(), n.begin());
    std::size_t const sym2 = rotational_symmetry(n);

    if (sym2 == 1) return 1;

    std::size_t const sym3 = pairs.symmetry();
    if (sym3 == 1) return 1;

    return std::gcd(std::gcd(sym1, sym2), sym3);
}

/******************************************************************************************/

string StateBase::dp() const {
    string ret;
    for (auto const &x : complexes.complex_indices) {
        for (auto i : indices<iseq>(x)) {
            auto s = x[i];
            for (auto p = sys->nicks[s] + 1; p != sys->nicks[s + 1] - 1; ++p) {
                auto q = pairs[p];
                if (q == p) ret.push_back('.');
                else if (q < sys->nicks[s] || q >= sys->nicks[s + 1]) {
                    auto t = find_if(x, [&](iseq z) {return q < sys->nicks[z + 1] && q >= sys->nicks[z];}) - begin_of(x);
                    ret.push_back(t < i ? ')' : '(');
                }
                else ret.push_back(q < p ? ')' : '(');
            }
            ret.push_back('+');
        }
        ret.back() = ' ';
    }
    NUPACK_REQUIRE(count(ret, ')'), ==, count(ret, '('), ret);
    ret.pop_back(); return ret;
}

/******************************************************************************************/

/// Returns a vector of (index with null bases) -> (index without null bases)
pair_data_type make_pairs_map(StateBase const &w) {
    pair_data_type map(len(w.pairs), -1);
    auto ppos = begin_of(map);
    auto mpos = 0;
    for (auto s : indices<iseq>(w.sys->strands)) {
        ++ppos;
        std::iota(ppos, ppos + len(w.sys->strands[s]) - 2, mpos);
        ppos += len(w.sys->strands[s]) - 2; mpos += len(w.sys->strands[s]) - 2;
        ++ppos;
    };
    return map;
}

/******************************************************************************************/

string StateBase::sequence() const {
    std::ostringstream os;
    for (auto const &x : complexes.complex_indices) for (auto const &s : x)
        dump_os(os, sys->strands[s].offset(1, -1), '+');
    string ret = os.str(); ret.pop_back(); return ret;
}

/******************************************************************************************/

System::System(StrandList const &v) {
    reserve(sum(v, len) + 2 * len(v));
    for (auto const &s : v) {
        if (front(s) != Base('_')) total_sequence.push_back(Base('_'));
        extend(total_sequence, s);
        if (back(s) != Base('_')) total_sequence.push_back(Base('_'));
    }
    nicks = {0};
    for (auto const &s : v) {
        auto const n = len(s) + iseq(back(s) != Base('_')) + iseq(front(s) != Base('_'));
        extend(strand_map, n, len(nicks) - 1);
        nicks.push_back(nicks.back() + n);
    }
    make_strands();
}

/******************************************************************************************/

}
