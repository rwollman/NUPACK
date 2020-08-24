#pragma once

#include "named_spec.h"
#include "nupack_invariants.h"

#include "structure_spec.h"
#include "tube_spec.h"

#include <thermo.h>

#ifdef JSONCPP_FOUND
#include <json/json.h>
#endif

#include <vector>
#include <string>
#include <utility>

namespace nupack {
  class PhysicalParams {
    public:
      DBL_TYPE temperature;
  };

  // Forward declare these so we can use them during decomposition
  class SingleParamResult;
  class PhysicalResult;
  class SequenceState;

  class SingleParamSpec : public NamedSpec {
    public:
      SingleParamSpec() : NamedSpec("-default", "specset"), off_targets(false) {}

      void set_parameters(const PhysicalParams & params) { this->params = params; }

      void add_tube(const TubeSpec & spec) { this->tubes.push_back(spec); }
      std::pair<int, int> add_structure(const StructureSpec & spec);

      void enumerate_off_targets(const std::map<std::string, int> & strands);
      void compute_target_matrices() { 
        for (auto & s : strucs) {
          s.compute_target_matrix(); 
          // s.serialize();
        }
      }

      StructureSpec & get_struc(std::string name);
      TubeSpec & get_tube(std::string name);

      std::pair<int, int> get_struc_id(std::string name) const;
      int get_tube_id(std::string name) const;

      const std::vector<StructureSpec> & get_strucs() const { return this->strucs; }
      const std::vector<OrderSpec> & get_orders() const { return this->orders; }

      const std::vector<int> & get_struc_ord_map() const { return this->struc_ord_map; }
      const std::vector<int> & get_ord_struc_map() const { return this->ord_struc_map; }
      const std::vector<TubeSpec> & get_tubes() const { return this->tubes; }
      const PhysicalParams & get_params() const { return this->params; }

      int get_max_depth() const;

      void decompose(const NupackInvariants & invars) { for (auto & s : strucs) s.decompose(invars); }
      void decompose_ppair(int i, int k, const SingleParamResult & res,
          const SequenceState & seqs, const NupackInvariants & invars);
      void add_off_target(int i, const SequenceSpec & seqs);

      SingleParamSpec get_depth(int depth, bool off_targets=false) const;
      bool eval_off_targets() const { return this->off_targets; }

    protected:
      std::vector<StructureSpec> & get_mod_strucs() const;

    private:
      bool off_targets;
      PhysicalParams params;

      std::vector<TubeSpec> tubes;
      std::vector<StructureSpec> strucs;
      std::vector<OrderSpec> orders;

      std::vector<int> struc_ord_map;
      std::vector<int> ord_struc_map;
  };

  class PhysicalSpec {
    public:
      PhysicalSpec() {}

      void add_param_spec(const SingleParamSpec & specification) { 
        specifications.push_back(specification);
      }
      const std::vector<SingleParamSpec> & get_param_specs() const { 
        return specifications; 
      }
      void add_off_target(int i_phys, int i_ord, const SequenceSpec & seqs) { 
        specifications[i_phys].add_off_target(i_ord, seqs);
      }

      void decompose(const NupackInvariants & invars) { for (auto & s : specifications) s.decompose(invars); }
      void decompose_ppair(int s, int i, int k, const PhysicalResult & res,
        const SequenceState & seqs, const NupackInvariants & invars);

      PhysicalSpec get_depth(int depth, bool off_targets=false) const;
      int get_max_depth() const;

#ifdef JSONCPP_FOUND
      Json::Value make_json_structures(const SequenceSpec & seqs,
          const NupackInvariants & invars) const;
      Json::Value make_json_tubes() const;
#endif

    private:
      std::vector<SingleParamSpec> specifications;
  };
}
