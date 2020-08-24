#include "objective_handler.h"
#include "pathway_utils.h"
#include "design_debug.h"
#include "algorithms.h"

#include <algorithm>

namespace nupack {
DBL_TYPE StructureObjective::get_defect(const EvalSpec & spec, 
    const EvalResult & res, WeightSpec & weightspec) {
  if (id == -1) {
    const auto cid = spec.physical.get_param_specs()[0].get_struc_id(name);

    id = cid.first;
    targ_id = cid.second;
  }

  NUPACK_CHECK(res.physical.get_results().size() > 0, 
      "Unevaluated result or invalid result in defect calculation");
  NUPACK_CHECK((int)res.physical.get_results()[0].get_strucs().size() > id,
      "Unevaluated result or invalid result in defect calculation");
  
  const auto & struc = res.physical.get_results().at(0).get_strucs().at(id);
  auto nuc_defects = struc.get_nuc_defects(targ_id);
  auto weighted_nuc_defects = weight_defects(nuc_defects, weightspec);

  DBL_TYPE defect = 0;
  for (auto & w : weighted_nuc_defects) defect += w.second;
  return defect / struc.size();
}

// ID = {physical id, tube id, structure id, strand id, domain id, index in structure, nucleotide id}

void get_variable_mut_weights(const std::vector<int> & var_ids, 
    WeightSpec & weightspec, std::vector<DBL_TYPE> & weights, std::vector<int> prefix, 
    const Map & defects, DBL_TYPE scale) {
  
  for (const auto & def : defects) {
    std::vector<int> nuc_ind = prefix;
    append(nuc_ind, def.first);
    
    auto var_id = var_ids.at(get_nuc_id(def.first));
    auto weight = weightspec.get_weight(nuc_ind);
    
    weights[var_id] += weight * def.second / scale;
  }
}

Map weight_defects(const Map nuc_defects, WeightSpec & weights, std::vector<int> f) {
  Map weighted;
  for (auto & def : nuc_defects) {
    const auto & id = def.first;
    const auto & defect = def.second;
    std::vector<int> full_id = f;
    append(full_id, id);
    DBL_TYPE weight = weights.get_weight(full_id);

    weighted[full_id] = weight * defect;
  }
  return weighted;
}

Map StructureObjective::weight_defects(Map nuc_defects, WeightSpec & weights) {
  return ::nupack::weight_defects(nuc_defects, weights, {get_physical_id(), -1});
}


void StructureObjective::get_variable_mut_weights(const EvalSpec & spec, 
    const EvalResult & res, WeightSpec & weightspec, std::vector<DBL_TYPE> & weights) {
  if (!Objective::satisfied(spec, res, weightspec)) {
    const auto & struc = res.physical.get_results().at(0).get_strucs().at(id);
    ::nupack::get_variable_mut_weights(spec.sequences.get_variables(), 
        weightspec, weights, {get_physical_id(), -1},
        struc.get_nuc_defects(targ_id), struc.size());
  }
}

DBL_TYPE TubeObjective::get_defect(const EvalSpec & spec, const EvalResult & res, 
    WeightSpec & weightspec) {
  if (id == -1) id = spec.physical.get_param_specs()[0].get_tube_id(name);

  NUPACK_CHECK(res.physical.get_results().size() > 0, 
      "Unevaluated result or invalid result in defect calculation");
  
  const auto & tube = res.physical.get_results().at(0).get_tubes().at(id);
  const auto & nuc_defects = tube.get_nuc_defects();
  auto weighted_nuc_defects = weight_defects(nuc_defects, weightspec);

  DBL_TYPE defect = 0;
  for (auto & def : weighted_nuc_defects) defect += def.second;
  return defect / tube.get_nuc_conc();
}


Map TubeObjective::weight_defects(Map nuc_defects, WeightSpec & weights) {
  return ::nupack::weight_defects(nuc_defects, weights, {get_physical_id()});
}

void TubeObjective::get_variable_mut_weights(const EvalSpec & spec,
    const EvalResult & res, WeightSpec & weightspec, std::vector<DBL_TYPE> & weights) {
  if (!Objective::satisfied(spec, res, weightspec)) {
    const auto & tube = res.physical.get_results().at(0).get_tubes().at(this->id);
    ::nupack::get_variable_mut_weights(spec.sequences.get_variables(), 
        weightspec, weights, {get_physical_id()},
        tube.get_nuc_defects(), tube.get_nuc_conc());
  }
}

std::vector<int> TubeObjective::get_tube_id() const {
  return std::vector<int>(1, this->id);
}

void CombinedObjective::resolve_names(const EvalSpec & spec) {
  for (auto & e_name : element_names) {
    auto struc_id = spec.physical.get_param_specs()[0].get_struc_id(e_name);

    if (struc_id.first >= 0) {
      // This is a structure
      struc_ids.push_back(struc_id.first);
      target_ids.push_back(struc_id.second);
    } else {
      // This is a tube, throws an exception otherwise
      // That behavior should change so that we can
      // actually search through to find names.
      // Not necessary for now.
      int tube_id = spec.physical.get_param_specs()[0].get_tube_id(e_name);
      tube_ids.push_back(tube_id);
    }
  }
}

DBL_TYPE CombinedObjective::get_defect(const EvalSpec & spec, const EvalResult & res,
    WeightSpec & weightspec) {
  if (!element_names.empty() && struc_ids.empty() && tube_ids.empty()) resolve_names(spec);

  NUPACK_CHECK(res.physical.get_results().size() > 0, 
      "Unevaluated result or invalid result in defect calculation");

  // calculate structure defects
  std::vector<Map> nuc_defects;
  std::vector<const StructureResult *> strucs;
  nuc_defects.reserve(this->struc_ids.size() + this->tube_ids.size());
  strucs.reserve(this->struc_ids.size());

  for (auto i = 0; i < this->struc_ids.size(); i++) {
    const auto & struc = 
        res.physical.get_results().at(0).get_strucs().at(this->struc_ids[i]);
    const auto & cur_nuc_defects = 
        struc.get_nuc_defects(this->target_ids[i]);
    nuc_defects.push_back(cur_nuc_defects);
    strucs.push_back(&struc);
  }
  
  // calculate tube defects
  std::vector<const TubeResult *> tubes;
  tubes.reserve(this->tube_ids.size());  
  for (auto & tube_id : tube_ids) {
    NUPACK_CHECK(res.physical.get_results()[0].get_tubes().size() > tube_id,
        "Unevaluated result or invalid result in defect calculation");

    const auto & tube = res.physical.get_results().at(0).get_tubes().at(tube_id);
    const auto & cur_nuc_defects = tube.get_nuc_defects();
    nuc_defects.push_back(cur_nuc_defects);
    tubes.push_back(&tube);
  }

  auto weighted_nuc_defects = weight_defects(nuc_defects, weightspec);
  
  // sum and normalize weighted defects
  DBL_TYPE total_defect = 0;
  for (auto i = 0; i < weighted_nuc_defects.size(); i++) {
    DBL_TYPE cur_defect = 0;
    for (auto & w : weighted_nuc_defects[i]) cur_defect += w.second;
    
    if (i < this->struc_ids.size()) {
      cur_defect /= (*strucs[i]).size(); 
    } else {
      cur_defect /= (*tubes[i - this->struc_ids.size()]).get_nuc_conc();
    }
    
    total_defect += cur_defect;
  }

  total_defect /= (this->tube_ids.size() + this->struc_ids.size());

  return total_defect;
}

std::vector<Map> CombinedObjective::weight_defects(
    std::vector<Map> nuc_defects, WeightSpec & weights) {
  std::vector<Map> weighted(nuc_defects.size());
  int physical_id = get_physical_id();

  for (auto it = nuc_defects.begin(); it != nuc_defects.end(); it++) {
    int pos = it - nuc_defects.begin();
    for (auto & m : *it) {
      const auto & id = m.first;
      DBL_TYPE & defect = m.second;
      
      std::vector<int> full_id;
      if (6 == id.size()) { // tube
        full_id = {physical_id};
      } 
      else if (5 == id.size()) { // structure
        full_id = {physical_id, -1};
      } 
      // else invalid
      append(full_id, id);
      
      DBL_TYPE weight = weights.get_weight(full_id);
      weighted[pos][full_id] = weight * defect;
    }
  }
  return weighted;
}

void CombinedObjective::get_variable_mut_weights(const EvalSpec & spec,
    const EvalResult & res, WeightSpec & weightspec, std::vector<DBL_TYPE> & weights) {
  if (!element_names.empty() && struc_ids.empty() && tube_ids.empty()) resolve_names(spec);

  const auto & var_ids = spec.sequences.get_variables();
  
  if (!Objective::satisfied(spec, res, weightspec)) {
    for (auto i = 0; i < struc_ids.size(); i++) {
      const auto & struc = res.physical.get_results().at(0).get_strucs().at(struc_ids[i]);
      ::nupack::get_variable_mut_weights(var_ids, weightspec, weights, 
          {get_physical_id(), -1}, struc.get_nuc_defects(target_ids[i]), struc.size());
      
    }

    for (auto & tube_id : tube_ids) {
      const auto & tube = res.physical.get_results().at(0).get_tubes().at(tube_id);
      ::nupack::get_variable_mut_weights(var_ids, weightspec, weights, 
          {get_physical_id()}, tube.get_nuc_defects(), tube.get_nuc_conc());
    }
  }
}

void CombinedObjective::serialize(const EvalSpec & spec, const EvalResult & res, 
    WeightSpec & weightspec, std::ostream & out, int indent, 
    std::string prefix) {
  Objective::serialize(spec, res, weightspec, out, indent, prefix);
  
  std::string ind_string = prefix + std::string(indent, ' ');
  
  out << ind_string << "including: ";
  out << to_string(element_names);
  out << std::endl;
}

ObjectiveHandler & ObjectiveHandler::operator=(const ObjectiveHandler & other) {
  this->copy(other);
  return *this;
}

void ObjectiveHandler::copy(const ObjectiveHandler & other) {
  if (&other != this) {
    this->weights = other.weights;
    for (auto & ob : other.objectives) objectives.push_back(ob->clone());
  }
}

std::vector<DBL_TYPE> ObjectiveHandler::get_defects(EvalSpec & spec, EvalResult & res) {
  std::vector<DBL_TYPE> defects;
  for (auto & ob : objectives)
    defects.push_back(ob->get_defect(spec, res, this->weights));
  return defects;
}

void ObjectiveHandler::get_variable_mut_weights(const EvalSpec & spec,
    const EvalResult & res, std::vector<DBL_TYPE> & weights) {
  const std::vector<int> & vars = spec.sequences.get_variables();
  weights.assign(vars.size(), 0);
  for (auto & ob : objectives)
    ob->get_variable_mut_weights(spec, res, this->weights, weights);
}

bool ObjectiveHandler::all_satisfied(const EvalSpec & spec, const EvalResult & res) {
  return all_of(satisfied(spec, res), [](const bool a) { return a; });
}

std::vector<bool> ObjectiveHandler::satisfied(const EvalSpec & spec, const EvalResult & res) {
  std::vector<bool> sat;
  for (auto & ob : objectives) 
    sat.push_back(ob->satisfied(spec, res, this->weights));
  return sat;
}

void ObjectiveHandler::short_serial(const EvalSpec & spec, const EvalResult & res,
    std::ostream & out, int indent, std::string prefix) {
  std::string ind_string(indent, ' ');

  ind_string = prefix + ind_string; 

  DBL_TYPE cur_defect = 0;
  out << ind_string;
  for (auto & ob : objectives) {
    cur_defect = ob->get_defect(spec, res, this->weights);
    out << cur_defect << ": " << ob->satisfied(spec, res, this->weights) << ",   ";
  }
  out << std::endl;
}

void ObjectiveHandler::serialize(const EvalSpec & spec, const EvalResult & res,
    std::ostream & out, int indent, std::string prefix) {
  std::string ind_string(indent, ' ');
  ind_string = prefix + ind_string; 

  DBL_TYPE mean_objective = 0;
  for (auto & ob : objectives) mean_objective += ob->get_defect(spec, res, this->weights);

  mean_objective /= this->objectives.size();
  out << ind_string << "mean objective: " << flt_format << mean_objective << std::endl;

  out << ind_string << "objectives:" << std::endl;
  for (auto & ob : objectives) 
    ob->serialize(spec, res, this->weights, out, indent + 4, prefix);
};

std::vector<std::string> ObjectiveHandler::get_names() {
  std::vector<std::string> names;
  for (auto & ob : objectives) names.push_back(ob->get_name());
  return names;
}

std::vector<DBL_TYPE> ObjectiveHandler::get_stops(const EvalSpec & spec,
    const EvalResult & res) {
  std::vector<DBL_TYPE> stops;
  for (auto & ob : objectives) stops.push_back(ob->get_stop(spec, res));
  return stops;
}
}
