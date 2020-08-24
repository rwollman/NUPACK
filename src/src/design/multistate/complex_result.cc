#include "complex_result.h"

#include "sequence_state.h"

#include "design_debug.h"
#include "pathway_utils.h"

#include "algorithms.h"

namespace nupack {

void OrderResult::clear_evaluated() {
  evaluated = false;
  eval_time = 0;
  pfunc = 0;
}

void OrderResult::evaluate(const SequenceState & seqs, const OrderSpec & spec,
    const PhysicalParams & params, const NupackInvariants & invars) {
  evaluated = true;
  strands = spec.get_strands();
  
  update_sequence(seqs);
  
  DBL_TYPE start_time = get_current_time();

  pfunc = pfuncFull(sequence.data(), 3, invars.material, 
      invars.dangle_type, invars.temperature - ZERO_C_IN_KELVIN, 0,
      invars.sodium, invars.magnesium, invars.use_long_helix) / spec.get_symmetry();

  eval_time = get_current_time() - start_time;
}

void OrderResult::update_sequence(const SequenceState & seqs) {
  sequence.clear();
  
  const auto & strand_seqs = seqs.get_strands();
  
  // build sequence
  int n_nucs = 0;
  for (const auto & c_str : strands) {
    NUPACK_CHECK(c_str >= 0 && c_str < strand_seqs.size(),
        "Invalid strand id " + to_string(c_str));
    const auto & curseq = strand_seqs[c_str].get_nucs();
    append(sequence, curseq); 
    sequence.push_back(STRAND_PLUS);
    n_nucs += curseq.size();
  }
  sequence.back() = -1;
  
  this->n_nucs = n_nucs;
}  
}
