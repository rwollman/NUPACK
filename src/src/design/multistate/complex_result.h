#pragma once

#include "physical_spec.h"

#include <vector>
#include <string>
#include <memory>

namespace nupack {
    
    /**
      * This class holds the evaluation results for a single complex without
      * a target structure (or decomposition tree).
      */
  class OrderResult {
    public:
      OrderResult() : n_nucs(0), pfunc(0), eval_time(0), evaluated(false) {}
      void evaluate(const SequenceState & seqs, const OrderSpec & spec,
          const PhysicalParams & params, const NupackInvariants & invars);

      void clear_evaluated();
      bool get_evaluated() const { return this->evaluated; }
      DBL_TYPE get_pfunc() const { return this->pfunc; }
      const std::vector<int> & get_strands() const { return this->strands; }
      DBL_TYPE get_eval_time() const { return this->eval_time; }
      int size() const { return this->n_nucs; }
      
      void update_sequence(const SequenceState & seqs);

    private:
      std::vector<int> strands;
      std::vector<int> sequence;
      int n_nucs;
      DBL_TYPE pfunc;
      DBL_TYPE eval_time;
      bool evaluated;
  };

}
