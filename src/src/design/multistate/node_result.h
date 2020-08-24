#pragma once

#include "pair_probabilities.h"
#include "physical_spec.h"

#include <vector>
#include <string>
#include <memory>

namespace nupack {
    class NodeResult;
    class StructureSpec;
    class StructureResult;
    
    /**
      * This class is used to build the decomposition tree for a single complex
      * the leaf nodes are used to perform the pair probability and partition
      * function evaluations. Parental nodes recursively merge the resulting
      * leaf properties to estimate the root-node properties.
      */
      
  using ChildPair = std::pair<NodeResult, NodeResult>;
  class NodeResult {
    public:
      NodeResult() {}

      DBL_TYPE get_pfunc() const { return pfunc_corrected; }

      DBL_TYPE collect_eval_times() const;
      PairProbs collect_pair_probs(const PhysicalParams & params,
          const NupackInvariants & invars) const;
      void map_defects(std::vector<DBL_TYPE> &nuc_defects);

      DBL_TYPE evaluate_dummy_pfuncs(const NodeSpec & spec,
          const StructureSpec & struc, const PhysicalParams & params, 
          const NupackInvariants & invars);
      void evaluate_leaf(const NodeSpec & spec, const SequenceState& seqs,
          const StructureSpec & strucspec,
          const StructureResult & strucres, const PhysicalParams & params,
          const NupackInvariants & invars);
      void evaluate(const NodeSpec & spec, const SequenceState & seqs,
          const StructureResult & struc, 
          const StructureSpec & strucspec, const PhysicalParams & params, 
          const NupackInvariants & invars);
      void init(const NodeSpec & spec);
      void serialize(std::ostream & out, int indent, const std::string & prefix, int & id) const;

      int get_n_leaves() const;
      int get_n_children() const { return children.size(); }

      void replace_node(const NodeResult & other, int k, 
          const PhysicalParams & params, const NupackInvariants & invars);

      /* Static utilities */
      static DBL_TYPE merge_pfuncs(const NodeResult & left,
          const NodeResult & right, const PhysicalParams & params, 
          const NupackInvariants & invars);

      int get_max_depth() const;

      friend class NodeSpec;

    protected:
      std::vector<NodeResult> children;
      std::vector<ChildPair> paired_children;
      int native_index(int i) const;
      void clear();
      void clear_children();

    private:
      void copy(const NodeResult & other);
      std::vector<int> eval_sequence;
      std::vector<int> nuc_ids;
      std::vector<int> to_full;
      std::vector<int> breaks;
      std::vector<bool> native;
      std::vector<DBL_TYPE> nuc_defects;

      PairProbs ppairs;
      DBL_TYPE pfunc_corrected;
      DBL_TYPE eval_time;
  };
}
