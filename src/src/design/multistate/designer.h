#pragma once

#include "design_spec.h"
#include "design_result.h"

#include <vector>
#include <string>

namespace nupack {

  class Designer {
    using Results = std::vector<DesignResult>;
    public:
      Designer(const DesignSpec & spec) : spec(spec) {};

      void evaluate_defect();
      void optimize();
      void optimize_tubes();
      void optimize_forest(PhysicalSpec & spec, Results & res);
      void optimize_leaves(PhysicalSpec & spec, Results & res);

      void evaluate(PhysicalSpec & spec, Results & res);
      void mutate_leaves(PhysicalSpec & spec, Results & res);
      void make_offspring(PhysicalSpec & spec, Results & res, int n=1);
      void redecompose(PhysicalSpec & spec, Results & parent, Results & child);
      void update_tabu(Results & old_gen, Results & new_gen);
      
      bool pick_best_par(Results & res);

      bool parent_satisfied(PhysicalSpec & curspec, Results & parent, Results & child);
      bool tubes_satisfied(PhysicalSpec & curspec, Results & tuberes, Results & treeres);
      bool all_satisfied(PhysicalSpec & curspec, Results & res);

      void add_off_targets(PhysicalSpec & curspec, Results & tuberes, Results & treeres);

      const Results & get_results() {return this->results; }
      
      void print_leaf_root(PhysicalSpec & spec, Results & res, std::string type, int index);
      void serialize_defects(std::ostream & out, PhysicalSpec & spec, 
          Results & res, int indent, std::string prefix);

    private:
      DesignSpec spec;
      Results results;
      
      void copy_sequences(Results & a, Results & b);
  };
}
