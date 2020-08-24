#include "design_spec.h"

namespace nupack {

#ifdef JSONCPP_FOUND
void DesignSpec::serialize_json(std::ostream & out, bool pprint) {
  Json::Value root;

  root = make_json_value();

  if (pprint) {
    out << root;
  } else {
    Json::FastWriter writer;
    out << writer.write(root);
  }
}
#endif

#ifdef JSONCPP_FOUND
Json::Value DesignSpec::make_json_value() const {
  Json::Value val;

  val["options"] = this->eval.options.make_json_value();
  val["structures"] = this->eval.physical.make_json_structures(
      this->eval.sequences, this->eval.options);
  val["tubes"] = this->eval.physical.make_json_tubes();


  std::vector<int> poss_nucs = this->constraints.get_possible_nucleotides();

  val["strands"] = this->eval.sequences.make_json_strands(
      poss_nucs, this->eval.options);
  val["domains"] = this->eval.sequences.make_json_domains(poss_nucs, 
      this->eval.options);
  return val;
}
#endif
}
