/** \file PairList.h
 * @brief Simple vector wrapper containing the indices of which base each base is paired to
 */
#pragma once
#include "../standard/Vec.h"
#include "../types/IO.h"
#include "../iteration/Range.h"
#include "../algorithms/Operators.h"
#include "../algorithms/Utility.h"
#include "../reflect/Hash.h"

namespace nupack {

/******************************************************************************************/

using pair_data_type = vec<iseq>;
using Nicks = small_vec<iseq>;

struct PairList : Indexable<PairList>, TotallyOrdered {
    using data_type = pair_data_type;
    using size_type = iseq;
    using value_type = value_type_of<data_type>;
    using iterator = iterator_of<data_type>;

    /************************************************************************************/

    PairList() = default;
    explicit PairList(usize s) : values(s) {reset();}

    PairList(std::string_view s) : values(io::to_pairs<data_type>(s)) {}
    PairList(char const *s) : values(io::to_pairs<data_type>(s)) {}
    PairList(std::string const &s) : values(io::to_pairs<data_type>(s)) {}

    PairList(data_type t) : values(std::move(t)) {}

    template <class Iter>
    PairList(Iter b, Iter e) : values(b, e) {}

    /************************************************************************************/

    bool operator<(PairList const &p) const {return values < p.values;}
    bool operator==(PairList const &p) const {return values == p.values;}
    std::size_t symmetry() const;

    /************************************************************************************/

    data_type values;
    NUPACK_REFLECT(PairList, values);

    auto & iter() {return values;}

    /************************************************************************************/

    template <class N=data_type> string dp(N && nicks={}) const {
        return len(nicks) ? io::to_dp(values, nicks) : io::to_dp(values);
    }

    /// Print with dot-parens -- for multistranded this won't put '+' in though
    friend std::ostream & operator<<(std::ostream &os, PairList const &p) {return os << "PairList('" << p.dp() << "')";}

    /// Delete all base values
    void reset() {std::iota(begin_of(values), end_of(values), 0);}
    bool empty() const {return values.empty();}

    template <class...Ts> void resize(Ts &&...ts) {values.resize(fw<Ts>(ts)...);}

    /************************************************************************************/

    void add_pair(size_type i, size_type j) {at(values, i) = j; at(values, j) = i;}

    /// Add base pair from i to j if it doesn't exist, delete it if it does exist
    void toggle_pair(size_type i, size_type j) {
        if (at(values, i) == i) {at(values, i) = j; at(values, j) = i;}
                         else {at(values, i) = i; at(values, j) = j;}
    }

    operator data_type const & () const {return values;}
    operator data_type () && {return std::move(values);}

    /// returns if two subsequences of the pairlist are equivalent once offsets are subtracted
    bool submatch(cspan I, cspan J) const {
        if (len(I) != len(J)) return false;
        for (auto k : indices(I))
            if (at(values, I[k]) + J[k] != at(values, J[k]) + I[k]) return false;
        return true;
    }

    size_type operator^(PairList const &p) const {return hamming_distance(values, p.values);}

    /// Expects a list of strand lengths where the lengths exclude null bases
    PairList with_null_bases(small_vec<iseq> const &strand_lengths) const;

    /************************************************************************************/

    /// Call a functor with the indices of each base pair
    template <class F> void for_each_pair(F &&f) const {
        izip(values, [&](auto i, auto j) {if (i < j) f(i, j);});
    }

    /// Number of base values
    size_type n_pairs() const {
        size_type out = 0;
        for_each_pair([&](auto, auto) {++out;});
        return out;
    }

    void rotate(std::ptrdiff_t i);

    /// Append an independent PairList() to the end of this one
    void append(PairList const &p);

    void throw_if_invalid() const;

    vec<std::array<value_type, 4>> pseudoknots() const;

    bool is_connected(small_vec<uint> const &nicks) const;
};

/******************************************************************************************/

/// Call a functor with the indices of each pseudoknot in a PairList()
template <class P, class F>
void for_pseudoknots(P const &values, F &&f) {
    izip(values, [&](auto i, auto j) {
        for (auto k = i; k < j; ++k) if (values[k] > j) f(i, j, k, values[k]);
    });
}

/******************************************************************************************/

// This code looks pretty smart but I don't have much confidence in it since it isn't used... *shrug*
template <class V, class P>
std::size_t pairing_symmetry(V const &v, P const &values) {
    std::size_t out = 1, z = len(values);
    NUPACK_REQUIRE(sum(v, len), ==, len(values));
    prime_factorization(rotational_symmetry(v), [&](auto const n) {
        auto const s = sum(view(v, 0, len(v) / (out * n)), len);
        auto const b = std::begin(values), e = std::end(values);
        if (std::equal(b, e - s, b + s, e, [=](auto i, auto j) {return (i + s) % z == j;}) // values[:-shift] + shift == values[shift:])
         && std::equal(b, b + s, e - s, e, [=](auto i, auto j) {return (i + z - s) % z == j;})) // values[:shift] + (size - shift) == values[-shift:]
            out *= n;
    });
    return out;
}

template <class ...Ts>
PairList join_pairs(PairList p, Ts const &...ts) {NUPACK_UNPACK(p.append(ts)); return p;}

/******************************************************************************************/

/**
 * @brief Call a functor for a sequence of PairList() each of which only differs by one base pair from the last
 * @param a Beginning PairList() --> functor is not called on this one
 * @param b Ending PairList() --> functor is not called on this one
 * @param f functor
 * @return iseq number of base pair moves to get from a to b
 */
template <class F>
iseq for_pairlists_between(PairList a, PairList const &b, F &&f) {
    iseq n = 0;
    for (auto i : indices(a)) if (a[i] != i && a[i] != b[i]) {
        if (n++) f(a);
        a.toggle_pair(i, a[i]); // cleave conflicting values
    }
    for (auto i : indices(a)) if (a[i] != b[i]) {
        if (n++) f(a);
        a.toggle_pair(i, b[i]); // add missing values
    }
    return n;
}

/******************************************************************************************/

}

namespace std {

template <>
struct hash<nupack::PairList> {
    size_t operator()(nupack::PairList const &p) const;
};

}
