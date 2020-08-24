#pragma once

#include "node_spec.h"
#include "named_spec.h"
#include "complex_spec.h"
#include "pair_probabilities.h"

#include "types.h"

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <iostream>

namespace nupack {
  class NupackInvariants;
  class PhysicalParams;
  class StructureResult;
  class SequenceSpec;
  class SequenceState;
  
    
  class StructureSpec : public NamedSpec {
    public:
      StructureSpec(const std::string & name, const std::string & struc);
      StructureSpec(const std::string & name, const std::vector<int> & struc, 
          const std::vector<int> & breaks, const std::vector<int> & strands);
      StructureSpec(const StructureSpec & other);
      
      StructureSpec(StructureSpec &&) = default;
      

      virtual void set_id(int id);
      StructureSpec & operator=(const StructureSpec & other);
      StructureSpec & operator=(StructureSpec && other) = default;
      
      void set_strand_names(const std::vector<std::string> & strand_names) { 
        this->strand_names = strand_names;
      }
      void set_no_target() { strucs.clear(); }

      void resolve_strand_names(const std::map<std::string, int> & namemap);
      void resolve_nuc_ids(const SequenceSpec & seqs);

      void set_strands(const std::vector<int> & strand_ids);
  
      void decompose(const NupackInvariants & invars) { tree->decompose(*this, invars); }
      void decompose_ppair(int k, const StructureResult & res, const SequenceState & seqs,
          const PhysicalParams & params, const NupackInvariants & invars);

      void merge(const StructureSpec & other);
      void compute_target_matrix();

      StructureSpec get_depth(int depth) const;
      
      int get_symmetry() const { return symmetry; }

      int size() const;

      const std::vector<int> & get_strands() const { return strands; }
      const std::vector<int> & get_nuc_ids() const { return nuc_ids; }
      const std::vector<int> & get_breaks() const { return breaks; }

      const std::vector<std::string> & get_strand_names() const { return strand_names; }
      const std::vector<vec_structure> & get_structures() const { return strucs; }
      const PairProbs & get_target() const { return target_matrix; }
      const std::vector<std::string> & get_names() const { return struc_names; }

      int get_n_nodes() const { return tree->get_n_nodes(); }
      int get_n_leaves() const { return tree->get_n_leaves(); }
      int get_max_depth() const { return tree->get_max_depth(); }
      bool is_forbidden(int i, int j) const;
      void add_forbidden(int i, int j) { forbidden.emplace_back(i, j); }

      void serialize(std::ostream & out = std::cout, int indent = 0,
          std::string prefix = "") const {target_matrix.serialize(out, target_matrix.get_n());}

      OrderSpec make_order() const;

      const NodeSpec & get_tree() const { return *(tree); }
      const std::vector<int> & get_domain_map() const { return domain_map; }
      const std::vector<int> & get_strand_map() const { return strand_map; }
      const std::vector<int> & get_struc_ids() const { return struc_ids; }
      

  
    private:
      void clone(const StructureSpec & other);
      void rotate();

      // Used in pre-resolved phase
      std::vector<std::string> strand_names;

      std::vector<vec_structure> strucs;
      PairProbs target_matrix;
      std::vector<int> strands;
      std::vector<int> breaks;
      
      std::vector<std::string> struc_names;
      std::vector<int> struc_ids;
      std::vector<int> nuc_ids;

      std::vector<int> domain_map;
      std::vector<int> strand_map;

      std::vector<std::pair<int, int> > forbidden;

      int symmetry;

      std::shared_ptr<NodeSpec> tree;
  };
  

}
