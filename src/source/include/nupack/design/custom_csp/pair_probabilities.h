#pragma once

#include "adapter.h"
#include "types.h"

#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <tuple>
#include <algorithm>

namespace nupack {
namespace custom_csp {

class StructureSpec;

struct PairProbTriple {
    PairProbTriple(int i, int j, real prob) : i(i), j(j), prob(prob) {}
    int max_index() const { return std::max(i, j); }

    int i;
    int j;
    real prob;
};

inline void swap(PairProbTriple & a, PairProbTriple & b) {
    using std::swap;
    swap(a.i, b.i);
    swap(a.j, b.j);
    swap(a.prob, b.prob);
}

inline bool operator<(const PairProbTriple & a, const PairProbTriple & b) {
    return std::make_tuple(a.i, a.j) < std::make_tuple(b.i, b.j);
}

inline bool operator==(const PairProbTriple & a, const PairProbTriple & b) {
    return a.i == b.i && a.j == b.j;
}

inline real operator-(const PairProbTriple & a, const PairProbTriple & b) {
    return a.prob - b.prob;
}

class PairProbs {
    mutable vec<PairProbTriple> probs;

public:
    PairProbs() {}
    PairProbs(const vec<int> & structure);
    PairProbs(const PairProbs & old, const vec<int> & remap);
    ~PairProbs() {}

    void clear();

    // Merge this with other, scaling other by other_scale and this by this_scale
    void merge(const PairProbs & other, real other_scale, real this_scale);

    void push_back(int i, int j, real ppair);
    int get_n() const;
    vec<real> get_mat(int n) const;
    vec<real> get_nuc_defects(vec<int> pairing) const;
    vec<real> get_nuc_defects(const PairProbs & target) const;
    vec<real> get_pair_probs(vec<int> i, vec<int> j) const;
    vec<std::pair<int, int>> get_inds() const;

    void serialize(std::ostream & out, int n) const;

    void clear_forbidden(const StructureSpec & spec);
};

}
}
