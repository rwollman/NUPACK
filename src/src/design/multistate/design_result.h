#pragma once

#include "eval_result.h"

#include <vector>

namespace nupack {
  class DesignSpec;
  
  class DesignResult {
    public:
      using Tabu = std::vector<std::vector<std::pair<int, int> > >;

      void init_random(const DesignSpec & spec);

      bool is_tabu(const std::vector<int> & vars);
      void add_tabu(const std::vector<int> & vars);
      void clear_tabu() { this->tabu.clear(); }

      std::vector<std::pair<int, int> > get_mutations(const std::vector<int> & vars);

      EvalResult eval;
      std::vector<DBL_TYPE> objectives;
      std::vector<bool> satisfied;

    private:
      Tabu tabu;
  };
}

