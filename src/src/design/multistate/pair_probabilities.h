#pragma once

#include <thermo.h>
#include "types.h"

#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <tuple>
#include <algorithm>

namespace nupack {
  class StructureSpec;
  
  struct PairProbTriple {
      PairProbTriple(int i, int j, DBL_TYPE prob) : i(i), j(j), prob(prob) {}
      int max_index() const { return std::max(i, j); } 

      int i;
      int j;
      DBL_TYPE prob;
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
  
  inline DBL_TYPE operator-(const PairProbTriple & a, const PairProbTriple & b) {
    return a.prob - b.prob;
  }

  class PairProbs {
    public:
      PairProbs() {}
      PairProbs(const std::vector<int> & structure);
      PairProbs(const PairProbs & old, const std::vector<int> & remap);
      ~PairProbs() {}
      
      void clear();
      
      // Merge this with other, scaling other by other_scale and this by this_scale
      void merge(const PairProbs & other, DBL_TYPE other_scale, DBL_TYPE this_scale);
      
      void push_back(int i, int j, DBL_TYPE ppair);
      int get_n() const;
      std::vector<DBL_TYPE> get_mat(int n) const;
      std::vector<DBL_TYPE> get_nuc_defects(std::vector<int> pairing) const;
      std::vector<DBL_TYPE> get_nuc_defects(const PairProbs & target) const;
      std::vector<DBL_TYPE> get_pair_probs(std::vector<int> i, std::vector<int> j) const;
      std::vector<std::pair<int, int>> get_inds() const; 
      
      void serialize(std::ostream & out, int n) const;
      
      void clear_forbidden(const StructureSpec & spec);
      
    private:
      mutable std::vector<PairProbTriple> probs;
  };
}
