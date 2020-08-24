#pragma once

#include "constraint_handler.h"
#include "objective_handler.h"

#include <vector>
#include <string>

namespace nupack {

  class DesignSpec {
    public:
      DesignSpec() {}
      
      EvalSpec eval;
      ConstraintHandler constraints;
      ObjectiveHandler objectives;

#ifdef JSONCPP_FOUND
      void serialize_json(std::ostream & out, bool pprint = false);

      Json::Value make_json_value() const;
#endif
  };
}

