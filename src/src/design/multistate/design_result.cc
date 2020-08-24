#include "design_result.h"
#include "design_debug.h"
#include "design_spec.h"

#include "algorithms.h"

#include <algorithm>

namespace nupack {

void DesignResult::init_random(const DesignSpec & spec) {
  const std::vector<int> vars = spec.constraints.init_random();
  NUPACK_CHECK(vars[0] >= 0, "Inconsistent constraints. No viable sequences found");
  this->eval.sequences.set_variables(vars, spec.eval.sequences);
}

bool DesignResult::is_tabu(const std::vector<int> & vars) {
  auto mutations = this->get_mutations(vars);
  return mutations.empty() || lin_contains(mutations, tabu);
}

void DesignResult::add_tabu(const std::vector<int> & vars) {
  std::vector<std::pair<int, int> > mutations = this->get_mutations(vars);
  this->tabu.push_back(mutations);
}

std::vector<std::pair<int, int> > DesignResult::get_mutations(const std::vector<int> & vars) {
  std::vector<std::pair<int, int> > mutations;
  const std::vector<int> & seqvars = this->eval.sequences.get_variables();
  
  NUPACK_CHECK(vars.size() == seqvars.size(),
      "variable size " + to_string(vars.size()) 
      + " does not match seqvar size " + to_string(seqvars.size()));
  
  for (auto i = 0; i < vars.size(); i++) {
    if (vars[i] != seqvars[i]) mutations.emplace_back(i, vars[i]);
  }
  return mutations;
}
}
