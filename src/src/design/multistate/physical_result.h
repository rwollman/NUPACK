#pragma once

#include "sequence_state.h"
#include "structure_result.h"
#include "complex_result.h"
#include "tube_result.h"

#include "physical_spec.h"
#include "pair_probabilities.h"

#include "types.h"

#include <vector>
#include <string>

namespace nupack {
  
  /**
   * This holds the physical results for a set of tubes using a
   * single material and at a single temperature (so using a single
   * parameter set).
   */
  class SingleParamResult {
    public:
      void replace_node(const SingleParamResult & other, 
          const SingleParamSpec & spec,
          int i, int k, const NupackInvariants & invars);

      void evaluate(const SequenceState & seqs, 
          const SingleParamSpec & spec,
          const NupackInvariants & invars);

      void print_res_files(const SingleParamSpec & spec,
          const NupackInvariants & invars, std::string file_prefix) const;
      void serialize(const SingleParamSpec & spec,
          const NupackInvariants & invars, std::ostream & out,
          int indent, std::string prefix) const;

      std::vector<DBL_TYPE> get_tube_defects();

      /*
       * Nucleotide defect i at position j
       * nuc_defects[i][j]
       */
      const std::vector<TubeResult> & get_tubes() const { return this->tubes; }
      const std::vector<StructureResult> & get_strucs() const { return this->strucs; }
      const std::vector<OrderResult> & get_orders() const { return this->orders; }

      int get_n_strucs();
      int get_n_tubes();
      int get_max_depth() const;

      DBL_TYPE get_eval_time() const;

      friend class SingleParamSpec;

    protected:
      std::vector<StructureResult> strucs;
      std::vector<OrderResult> orders;

      std::vector<TubeResult> tubes;

    private:
      PhysicalParams params;
      std::vector<DBL_TYPE> order_pfuncs;
      std::vector<bool> order_included;

      Map nuc_defects;
  };

  /**
   * This class collects SingleParamResults.
   */
  class PhysicalResult {
    public:
      void print_res_files(const PhysicalSpec & spec,
          const NupackInvariants & invars, std::string file_prefix) const;
      void replace_node(const PhysicalResult & other, 
          const PhysicalSpec & spec, int s, int i, int k, 
          const NupackInvariants & invars);
      void evaluate(const SequenceState & seqs, const PhysicalSpec & spec,
          const NupackInvariants & invars);

      void serialize(const PhysicalSpec & spec, NupackInvariants & invars,
          std::ostream & out, int indent = 0, std::string prefix = "") const;

      const std::vector<SingleParamResult> & get_results() const { return this->results; }
      int get_max_depth() const;
      DBL_TYPE get_eval_time() const;

      friend class PhysicalSpec;

    protected:
      std::vector<SingleParamResult> results;
  };

};

