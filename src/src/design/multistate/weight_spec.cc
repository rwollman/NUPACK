#include "weight_spec.h"

namespace nupack {
WeightSpec::WeightSpec() {}

void WeightSpec::add_weight(const std::vector<std::string> & spec, DBL_TYPE weight) {
  this->specs.push_back(spec);
  this->weight.push_back(weight);
  this->cached_weight.clear();
}

void WeightSpec::resolve_names(const std::vector<std::map<std::string, int> > & maps) {
  int j;
  int k;
  int m_depth = maps.size();
  int n_specs = this->specs.size();
  this->resolved_ids.resize(n_specs);

  for (auto i = 0; i < n_specs; i++) {
    auto spec_depth = this->specs[i].size();
    for (k = 0, j = 0; k < spec_depth && j < m_depth; j++) {
      auto it = maps[j].find(this->specs[i][k]);
      auto cur_id = -1;
      if (it != maps[j].end()) {
        cur_id = it->second;
        k++;
      }
      this->resolved_ids[i].push_back(cur_id);
    }

    for (; j < m_depth; j++) {
      this->resolved_ids[i].push_back(-1);
    }
  }
  this->cached_weight.clear();
}

DBL_TYPE WeightSpec::get_weight(const std::vector<int> & spec) {
  auto cache_it = this->cached_weight.find(spec);
  DBL_TYPE weight;
  if (cache_it != this->cached_weight.end()) {
    return cache_it->second;
  } else {
    auto n_specs = this->resolved_ids.size();
    weight = 1.0;
    for (auto i = 0; i < n_specs; i++) {
      auto match = true;
      for (auto j = 0; j < this->resolved_ids[i].size(); j++) {
        if (this->resolved_ids[i][j] != -1 &&
            this->resolved_ids[i][j] != spec[j]) {
          match = false;
          break;
        }
      }

      if (match) {
        weight *= this->weight[i];
      }
    }
    this->cached_weight[spec] = weight;
  }

  return weight;
}
}
