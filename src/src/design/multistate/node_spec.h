#pragma once

#include <memory>
#include <vector>
#include <utility>
#include <iostream>

namespace nupack {
  class NupackInvariants;
  class SplitSet;
  class SequenceState;
  class StructureSpec;
  class StructureResult;
  class PhysicalParams;
  class NodeResult;
 
  using Limits = std::pair<int, int>;
  using Assume = std::pair<int, int>;
  
  class NodeSpec {  
    public:
      NodeSpec(const std::vector<Limits> & lims, const std::vector<Assume> assume = {}) : 
           lims(lims), assume(assume) {}

      void split(const SplitSet & splits);
      
      void decompose(const StructureSpec & spec, const NupackInvariants & invars);
      void decompose_ppair(const NodeResult & res, const SequenceState & seqs,
          StructureSpec & spec, const StructureResult & strucres, 
          const PhysicalParams & params, const NupackInvariants & invars);
      void decompose_ppair_at(int k, const NodeResult & res, const SequenceState & seqs,
          StructureSpec & strucspec, const StructureResult & strucres,
          const PhysicalParams & params, const NupackInvariants & invars);

      void forbid_children(StructureSpec & strucspec, const NodeResult & res,
          const PhysicalParams & params, const NupackInvariants & invars) const;
      
      std::vector<int> get_nuc_ids(const StructureSpec & struc, 
          const NupackInvariants & invars) const;
      std::vector<int> get_breaks(const StructureSpec & struc, 
          const NupackInvariants & invars) const;
      std::vector<int> get_native_map(const StructureSpec & struc, 
          const NupackInvariants & invars) const;
      std::vector<bool> get_native(const StructureSpec & struc,
          const NupackInvariants & invars) const;
      std::vector<int> get_to_node(const StructureSpec & struc,
          const NupackInvariants & invars,
          std::vector<int> to_full = {}) const;
      std::vector<std::vector<int> > get_structures(const StructureSpec & struc,
          const NupackInvariants & invars) const;

      NodeSpec get_depth(int depth) const;
      const NodeSpec & get_child(int i) const;
      
      int get_n_nodes() const;
      int get_n_leaves() const;
      int get_n_children() const { return children.size(); }
      int get_max_depth() const;
      
      const std::vector<Limits> & get_lims() const { return lims; }
      const std::vector<Assume> & get_assume() const { return assume; }

      void print_leaves() const;
      void print_decomposition() const;

    protected:
      std::vector<NodeSpec> & get_children() { return children; }
      void print_node_decomp(std::ostream & out=std::cerr) const;

    private:
      void clear();
      void clear_children() { children.clear(); }
      
      std::vector<Limits> lims;
      std::vector<Assume> assume;
      std::vector<NodeSpec> children;
  };
}
