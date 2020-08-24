#pragma once

#include "physical_spec.h"
#include "pair_probabilities.h"
#include "node_result.h"

#include "types.h"

#include <vector>
#include <string>
#include <memory>

namespace nupack {
  class NodeResult;
    /**
      * This class holds the evaluation results for a single complex that
      * may have a decomposition tree.
      */
  class StructureResult {
    public:
      StructureResult() : defect(DBL_MAX) {}

      void evaluate(const SequenceState & seqs, 
          const StructureSpec & spec,
          const PhysicalParams & params, const NupackInvariants & invars);

      int get_n_leaves() const { return this->tree.get_n_leaves(); }
      const Map & get_nuc_defects(int i_target) const;

      DBL_TYPE get_pfunc(const PhysicalParams & params) const;
      DBL_TYPE get_structural_defect(int i_target) const;
      DBL_TYPE get_normalized_defect(int i_target) const;

      void serialize(const StructureSpec & spec, const PhysicalParams & params,
          const NupackInvariants & invars, std::ostream & out, 
          int indent = 0, std::string prefix = "") const;

      const std::vector<int> & get_sequence() const { return this->sequence; }
      const std::vector<int> & get_strands() const { return this->strands; }
      const std::vector<int> & get_fake_sequence() const { return this->f_sequence; }

      std::string get_structure_string(int i_struc = 0) const;
      std::string get_sequence_string(const NupackInvariants & invars) const;

      int size() const { return this->sequence.size(); }
      int get_max_depth() const { return this->tree.get_max_depth(); }

      void replace_node(const StructureResult & other, int k,
          const PhysicalParams & params, const NupackInvariants & invars);

      void update_defects();

      void print_strucfiles(const StructureSpec & spec, const PhysicalParams & params,
          const NupackInvariants & invars, std::string file_prefix) const;

      static DBL_TYPE get_energy(const std::vector<int> & struc, 
          const std::vector<int> & breaks, 
          const std::vector<int> & seq, int sym, const PhysicalParams & params,
          const NupackInvariants & invars);

      DBL_TYPE get_eval_time() const { return this->eval_time; }

      friend class StructureSpec;

    protected:
      NodeResult tree;

    private:
      void copy(const StructureResult & other);
      std::vector<int> sequence;
      std::vector<int> f_sequence;
      std::vector<int> nuc_ids;
      std::vector<std::vector<int> > structures;
      std::vector<int> struc_ids;
      std::vector<int> breaks;
      std::vector<int> strands;
      std::vector<int> strand_map;
      std::vector<int> domain_map;

      DBL_TYPE defect;
      std::vector<DBL_TYPE> defects;
      std::vector<std::vector<DBL_TYPE> > pos_defects;
      std::vector<Map> nuc_defects;

      std::vector<DBL_TYPE> struc_energies;
      DBL_TYPE pfunc;
      PairProbs ppairs;
      PairProbs target;
      int symmetry;
      
      DBL_TYPE eval_time;
  };
}
