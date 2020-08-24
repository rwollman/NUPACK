#pragma once

#include "types.h"

#include <vector>
#include <string>
#include <utility>

namespace nupack {
  class NupackInvariants;
  class SplitSet;
  class StructureSpec;
  class PairProbs;
  
  namespace StructureUtils {
    using structure_pair = std::pair<vec_structure, std::vector<int>>;
    
    structure_pair dpp_to_pairs(const std::string & struc);
    std::string pairs_to_dpp(const vec_structure & pairs, 
        const std::vector<int> & breaks);
    void get_strand_lengths(const vec_structure & pairs, 
        const std::vector<int> & breaks, std::vector<int> & lengths);
    void get_strand_lengths(const std::string & struc, 
        std::vector<int> & lengths);
    std::vector<SplitSet> get_consistent_splits(
        const std::vector<std::vector<int> > & strucs, 
        const NupackInvariants & invars);
    SplitSet get_minimal_splits(const StructureSpec & struc,
        const std::vector<SplitSet> & set, const PairProbs & ppairs, int n, 
        const NupackInvariants & invars);
  };
}
