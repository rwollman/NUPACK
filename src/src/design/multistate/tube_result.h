#pragma once

#include "physical_spec.h"

#include "equilibrium_concentrations.h"

#include "types.h"

#include <vector>
#include <string>
#include <memory>

namespace nupack {
  class StructureResult;    
  class OrderResult;    
  
  struct ConcSolverParams {
    static constexpr int n_points = 1;
    static constexpr int max_iters = 10000;
    static constexpr double tol = 1e-8;
    static constexpr double delta_bar = TRUST_REGION_DELTABAR;
    static constexpr double eta = TRUST_REGION_ETA;
    static constexpr double min_delta = 1e-12;
    static constexpr int max_trial = 1000;
    static constexpr int perturb_scale = 100;
    static constexpr int quiet = 1;
    static constexpr int write_log_file = 0;
    static constexpr char *log_file = NULL;
    static constexpr int seed = 0; // Don't reseed the RNG!
  };
  
  /**
   * This class holds the physical results for a tube of interacting
   * nucleic acids.
   */
  class TubeResult {
    public:
      TubeResult() : defect(DBL_MAX), nuc_conc(1.0) {};

      void evaluate(const std::vector<StructureResult> & strucs, 
          const std::vector<OrderResult> & orders, 
          const std::vector<int> & ord_struc_map,
          const TubeSpec & spec, const PhysicalParams & params,
          const NupackInvariants & invars);
      const std::vector<DBL_TYPE> & get_concentrations() const { return this->x; }

      const Map get_nuc_defects() const { 
        return this->nucleotide_defects;
      }
      DBL_TYPE get_nuc_conc() const { return this->nuc_conc; }

      void serialize(const TubeSpec & spec, 
          const std::vector<StructureSpec> & strucspecs,
          const std::vector<int> & order_specs,
          const PhysicalParams & params,
          const NupackInvariants & invars, 
          std::ostream & out,
          int indent = 0, std::string prefix = "") const;
      DBL_TYPE get_defect() const { return this->defect; }
      
      void clear_state();

    private:
      std::vector<std::vector<DBL_TYPE> > target_x;
      std::vector<DBL_TYPE> x;
      Map nucleotide_defects;
      DBL_TYPE defect;
      DBL_TYPE nuc_conc;
  };

  
}
