#pragma once

#include "physical_result.h"
#include "symmetry_calc.h"

#include <vector>
#include <string>

namespace nupack {
  struct EvalResult {
    PhysicalResult physical;
    SequenceState sequences;
    
    SymmetryResult symmetries;
  };
}
