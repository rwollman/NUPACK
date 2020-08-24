#pragma once

#include "symmetry_calc.h"
#include "physical_spec.h"
#include "sequence_spec.h"

#include "nupack_invariants.h"

#include <vector>
#include <string>

namespace nupack {
  class NupackInvariants;
  class PhysicalSpec;
  class SequenceSpec;
  class SymmetrySpec; 

  class EvalSpec {
    public:
      // Use default constructors/destructors
      EvalSpec() {}
      
      NupackInvariants options;
      PhysicalSpec physical;
      SequenceSpec sequences;
      SymmetrySpec symmetries;
    
    private:
      void copy(const EvalSpec & other);
  };
}