#pragma once

#include "physical_spec.h"
#include "objective_handler.h"

#include <string>
#include <vector>

namespace nupack {
  class SymmetryObjective : public Objective_CRTP<SymmetryObjective> {
    public:
      SymmetryObjective(std::string name, DBL_TYPE stop) :
          Objective_CRTP(name, "SSM", stop), id(0) {}

      DBL_TYPE get_defect(const EvalSpec & spec,
          const EvalResult & res, WeightSpec & weightspec);
      void get_variable_mut_weights(const EvalSpec & spec,
          const EvalResult & res, WeightSpec & weightspec,
          std::vector<DBL_TYPE> & weights);

      bool satisfied(const EvalSpec & spec, const EvalResult & res,
          WeightSpec & weightspec);

    private:
      int id;
  };
};

