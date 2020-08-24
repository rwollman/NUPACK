#include "complex_spec.h"

#include <algorithm>

namespace nupack {

// rotate strand list into a lexicographically minimal form
OrderSpec::OrderSpec(const std::vector<int> & strands) {
  int n_strs = strands.size();
  int symmetry = 1;
  std::vector<int> finalstrands = strands;
  auto strandcopy = strands;
  for (auto i_str = 1; i_str < n_strs; i_str++) {
    std::rotate(strandcopy.begin(), strandcopy.begin() + 1, strandcopy.end());

    if (strandcopy < finalstrands) finalstrands = strandcopy; 
    if (strandcopy == strands) {
      symmetry = n_strs / i_str;
      break;
    }
  }
  this->strands = finalstrands;
  this->symmetry = symmetry;
}
    
}
