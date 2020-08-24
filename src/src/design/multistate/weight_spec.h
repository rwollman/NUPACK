#pragma once

#include "types.h"

#include <shared.h>

#include <map>
#include <vector>
#include <string>

namespace nupack {
  class WeightSpec {
    public:
      WeightSpec();

      void add_weight(const std::vector<std::string> & spec, DBL_TYPE weight);
      /*
       * Resolve the names using a vector of name-maps. In our case,
       * this is currently:
       * parameterset, tube, structures, strands, domains
       *
       * This is written to be generic across index depth in case
       * further levels of indexing are necessary. The resolution is
       * required for efficient weight lookups.
       */
      void resolve_names(const std::vector<std::map<std::string, int> > & maps);

      void serialize(std::ostream & out) const;

      DBL_TYPE get_weight(const std::vector<int> & spec);

    protected:
      /* Raw specifications */
      std::vector<std::vector<std::string> > specs;
      std::vector<DBL_TYPE> weight;

      /* Resolved ids */
      std::vector<std::vector<int> > resolved_ids;

      /* Cached weights */
      Map cached_weight;
  };
};

