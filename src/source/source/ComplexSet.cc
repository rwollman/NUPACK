#include <nupack/state/ComplexSet.h>
#include <nupack/iteration/Transform.h>

namespace nupack {

/******************************************************************************************/

void ComplexSet::update_join_rates(Index i, BaseMat<real> const &new_join_rates) {
    auto const x = strand_map[i].x;
    complex_rates.increment(x, new_join_rates - join_rates[i]);
    x_rates_sq.update(x, times_transpose(complex_rates[x]));
    join_rates[i] = new_join_rates;
}

/******************************************************************************************/

void ComplexSet::rotate(Index s) {
    auto const x = strand_map[s].x, pos = strand_map[s].pos; // complex and location within
    std::rotate(complex_indices[x].begin(), complex_indices[x].begin() + pos, complex_indices[x].end());
    for (auto i : indices<Index>(complex_indices[x])) strand_map[complex_indices[x][i]].pos = i;
};

/******************************************************************************************/

bool ComplexSet::check() const {
    for (auto p : indices<Index>(strand_map))
        if (p != complex_indices[strand_map[p].x][strand_map[p].pos]) return false;
    return true;
}

/******************************************************************************************/

void ComplexSet::register_join(Index s1, Index s2) {
    auto x1 = strand_map[s1].x, x2 = strand_map[s2].x;
    auto pos1 = strand_map[s1].pos, pos2 = strand_map[s2].pos;
    auto & xi1 = complex_indices[x1], & xi2 = complex_indices[x2];

    xi1 = catted<vec<Index>>(xi1.begin(), xi1.begin() + pos1, xi2.begin() + pos2, xi2.end(),
                             xi2.begin(), xi2.begin() + pos2, xi1.begin() + pos1, xi1.end());

    for (auto i : indices<Index>(xi1)) {
        strand_map[complex_indices[x1][i]].x = x1;
        strand_map[complex_indices[x1][i]].pos = i;
    }

    complex_rates.increment(x1, complex_rates[x2]);
    x_rates_sq.update(x1, times_transpose(complex_rates[x1]));
    complex_rates.swap_erase(x2); x_rates_sq.swap_erase(x2);

    if (x2 + 1 != complex_indices.size()) for (auto i : complex_indices.back()) strand_map[i].x = x2;
    swap_erase(complex_indices, x2);
}

/******************************************************************************************/

//// Given two strands indexed at "p" and "k", split the complex they belong to in 2.
/// If k.pos > p.pos, make a new complex from [p:k)
/// If p.pos > k.pos, make a new complex from [k:p)
void ComplexSet::register_split(Index p, Index k) {
    // Loops must be from same complex
    NUPACK_ASSERT(check());
    if (strand_map[p].x != strand_map[k].x) {
        print(strand_map);
        print(*this);
    }
    NUPACK_REQUIRE(strand_map[p].x, ==, strand_map[k].x, "tried to split two strands from different complexes");
    NUPACK_REQUIRE(complex_indices[strand_map[p].x][strand_map[p].pos], ==, p);
    NUPACK_REQUIRE(complex_indices[strand_map[k].x][strand_map[k].pos], ==, k);

    Index const old_ix = strand_map[p].x, new_ix = complex_indices.size();

    Index i1, i2; std::tie(i1, i2) = std::minmax(strand_map[k].pos, strand_map[p].pos);

    complex_indices.emplace_back();
    auto &oldx = complex_indices[old_ix], &newx = complex_indices[new_ix];
    NUPACK_REQUIRE(i2, <, len(oldx));

    newx.assign(oldx.begin() + i1, oldx.begin() + i2);
    oldx.erase(oldx.begin() + i1, oldx.begin() + i2);

    auto diff = BaseMat<real>(zero);
    for (auto j : newx) diff += join_rates[j];
    complex_rates.emplace_back(diff); complex_rates.increment(old_ix, -diff);

    x_rates_sq.emplace_back(times_transpose(complex_rates.back()));
    x_rates_sq.update(old_ix, times_transpose(complex_rates[old_ix]));

    for (auto j : indices<Index>(oldx)) {
        strand_map[oldx[j]].x = old_ix;
        strand_map[oldx[j]].pos = j;
    }
    for (auto j : indices<Index>(newx)) {
        strand_map[newx[j]].x = new_ix;
        strand_map[newx[j]].pos = j;
    }
}

/******************************************************************************************/

/// Search for "t" in the cumulative sum function C(f1) + s * C(f2)
template <class Matrix>
auto fenwick_find2(Fenwick<Matrix> const &f1, Fenwick<Matrix> const &f2, uint i, uint j, real t, real const &scale) {
    return fenwick_find<real>(t, 0, f1.total()(i, j) * scale - f2.total()(i, j), f1.size(), [&](auto k) {
        return f1.tree[k](i, j) * scale - f2.tree[k](i, j);
    });
}

/******************************************************************************************/

JoinMove ComplexSet::get_join_move_nondimensional(double r) const {
    JoinMove ret;
    Index x1, x2; // Complex indices of merging complexes
    double margin;

    // Choose a pair of bases from the total base by base (bxb) rate matrix
    auto const bxb = la::eval(times_transpose(complex_rates.total()) - x_rates_sq.total());
    auto const tup = find_cumulative_mat(bxb, r);
    ret.b1 = Base::from_index(first_of(tup));
    ret.b2 = Base::from_index(second_of(tup));
    r = third_of(tup);

    // Choose the first complex
    std::tie(x1, margin) = fenwick_find2(complex_rates, x_rates_sq, ret.b1, ret.b2, r, complex_rates.total()(ret.b2, ret.b1));

    // Choose the second complex, which can't equal the first
    // Scale down the running margin to work in the x2 space

    margin /= complex_rates[x1](ret.b1, ret.b2);

    auto x_rates_mod = complex_rates; x_rates_mod.zero(x1);
    std::tie(x2, margin) = x_rates_mod.find(margin, [i=ret.b2, j=ret.b1](auto const &m){return m(i, j);});
    NUPACK_REQUIRE(x1, !=, x2, *this);

    // Choose the strand in the second complex
    NUPACK_REQUIRE(margin, <= ,complex_rates[x2](ret.b2, ret.b1));
    auto const s2 = find_if(complex_indices[x2],
        [&](auto const &s) {return minus_if(margin, join_rates[s](ret.b2, ret.b1));});
    NUPACK_ASSERT(s2 != end_of(complex_indices[x2]));
    ret.o2 = strand_map[*s2].loop;

    // Scale back out
    margin *= complex_rates[x1](ret.b1, ret.b2) / join_rates[*s2](ret.b2, ret.b1);

    // Choose the strand in the first complex
    NUPACK_REQUIRE(margin, <=, complex_rates[x1](ret.b1, ret.b2));
    auto const s1 = find_if(complex_indices[x1], [&](auto const &s) {
        return minus_if(margin, join_rates[s](ret.b1, ret.b2));
    });
    NUPACK_ASSERT(s1 != end_of(complex_indices[x1]));
    ret.o1 = strand_map[*s1].loop;

    NUPACK_REQUIRE(margin, <=, join_rates[*s1](ret.b1, ret.b2));
    ret.margin = margin * join_rates[*s2](ret.b2, ret.b1);
    ret.scale = join_rates[*s2](ret.b2, ret.b1);

    NUPACK_REQUIRE(ret.o1, !=, ret.o2, *this);
    return ret;
}

/******************************************************************************************/

}
