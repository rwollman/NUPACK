#include "physical_spec.h"
#include "sequence_spec.h"
#include "pair_probabilities.h"
#include "structure_utils.h"

#include "physical_result.h"
#include "design_debug.h"
#include "utils.h"
#include "algorithms.h"

#include <json/json.h>

#include <sys/time.h>

#include <set>
#include <map>
#include <algorithm>
#include <list>
#include <ctime>
#include <cfloat>

namespace nupack {
std::pair<int, int> SingleParamSpec::add_structure(const StructureSpec & spec) {
  OrderSpec temp_ord = spec.make_order();

  unsigned int ord_id = 0;
  while (ord_id < orders.size() && orders[ord_id] != temp_ord) {
    ord_id++;
  }
  int struc_id = strucs.size();
  std::pair<int, int> ret(-1, -1);

  if (ord_id < orders.size() && ord_struc_map[ord_id] >= 0) {
    struc_id = ord_struc_map[ord_id];
    strucs[struc_id].merge(spec);
    NUPACK_DEBUG("Merging: " << spec.get_name() << " at " << struc_id);
    ret.first = ord_id;
    ret.second = strucs[struc_id].get_structures().size() - 1;
  } else {
    if (ord_id == orders.size()) {
      ord_struc_map.push_back(strucs.size());
      orders.push_back(temp_ord);
    } else {
      ord_struc_map[ord_id] = strucs.size();
    }

    struc_ord_map.push_back((int)ord_id);
    strucs.push_back(spec);

    ret.first = (int)ord_id;
    ret.second = 0;
  }

  return ret;
}

SingleParamSpec SingleParamSpec::get_depth(int depth, bool off_targets) const {
  SingleParamSpec res(*this);
  res.strucs.clear();
  for (auto & str : strucs) {
    StructureSpec temp = str.get_depth(depth);
    res.strucs.push_back(std::move(temp));
  }

  res.off_targets = off_targets;
  return res;
}

int SingleParamSpec::get_max_depth() const {
  using str = typename decltype(strucs)::value_type;
  return std::max_element(strucs.begin(), strucs.end(), smaller_depth<str>)->get_max_depth();
}

void SingleParamSpec::enumerate_off_targets(
    const std::map<std::string, int> & strands) {
  std::map<OrderSpec, int> off_target_map;
  int n_on_target = this->orders.size();
  NUPACK_DEBUG("N on targets " << n_on_target);
  
  for (auto i = 0; i < n_on_target; i++) {
    auto & strands = this->orders[i].get_strands();
    NUPACK_CHECK(!contains(strands, off_target_map), "Duplicate strand orderings.")
    off_target_map[strands] = i;
  }

  for (auto & tube : tubes) {
    tube.enumerate_off_targets(off_target_map);
    tube.resolve_listed_offtargets(strands, off_target_map);
  }

  this->orders.resize(off_target_map.size(), std::vector<int>());
  this->ord_struc_map.resize(off_target_map.size(), -1);

  for (auto & ot : off_target_map) {
    this->orders[ot.second] = ot.first;
    if (ot.second >= n_on_target) this->ord_struc_map[ot.second] = -1;
  }
  NUPACK_DEBUG("Total orders: " << this->ord_struc_map.size());
}

StructureSpec & SingleParamSpec::get_struc(std::string name) {
  for (auto & str : strucs) {
    if (name == str.get_name()) return str;
  }
  NUPACK_ERROR("Invalid structure " + name);
}

std::pair<int, int> SingleParamSpec::get_struc_id(std::string name) const {
  int i = 0; 
  int j = 0;
  for ( ; i < this->strucs.size(); i++) {
    const auto & cnames = this->strucs[i].get_names();
    for (j = 0; j < cnames.size(); j++) {
      if (name == cnames[j]) break;
    }

    if (j < cnames.size() && name == cnames[j]) break;
  }

  std::pair<int, int> ret(-1, -1);
  if (i < this->strucs.size()) {
    ret.first = i;
    ret.second = j;
  }
  return ret;
}

int SingleParamSpec::get_tube_id(std::string name) const {
  for (auto i = 0; i < this->tubes.size(); i++) {
    if (this->tubes[i].get_name() == name) return i;
  }
  NUPACK_ERROR("Tube " + name + " not found");
}

void SingleParamSpec::decompose_ppair(int i, int k, const SingleParamResult & res, 
    const SequenceState & seqs, const NupackInvariants & invars) {
  NUPACK_CHECK(i < this->strucs.size(), "Invalid index " + to_string(i) + 
      " in decompose_ppair");
  this->strucs[i].decompose_ppair(k, res.strucs[i], seqs, this->params, invars);
}

void SingleParamSpec::add_off_target(int i_ord, const SequenceSpec & seqs) {
  std::vector<int> strands = this->orders.at(i_ord).get_strands();
  std::vector<int> breaks;
  std::string name = "S-"; 
  unsigned int n_nucs = 0;
  for (auto s : strands) {
    auto & str = seqs.get_strands()[s];
    n_nucs += str.size();
    name += str.get_name();
    breaks.push_back(n_nucs);
  }
  std::vector<int> struc(n_nucs, -1);
  breaks.pop_back();

  StructureSpec newstruc(name, struc, breaks, strands);
  newstruc.set_no_target();
  newstruc.resolve_nuc_ids(seqs);

  NUPACK_CHECK(ord_struc_map[i_ord] < 0, "Reading structure that is present: " + 
      to_string(i_ord));

  this->strucs.push_back(newstruc);
  this->struc_ord_map.push_back(i_ord);
  this->ord_struc_map[i_ord] = this->strucs.size() - 1;
}

int PhysicalSpec::get_max_depth() const {
  using spec = typename decltype(specifications)::value_type;
  return std::max_element(specifications.begin(), specifications.end(), 
      smaller_depth<spec>)->get_max_depth();
}

PhysicalSpec PhysicalSpec::get_depth(int depth, bool off_targets) const {
  PhysicalSpec res(*this);
  res.specifications.clear();
  for (auto & spec : specifications) {
    SingleParamSpec temp = spec.get_depth(depth, off_targets);
    res.specifications.push_back(temp);
  }
  return res;
}

void PhysicalSpec::decompose_ppair(int s, int i, int k, const PhysicalResult & res,
    const SequenceState & seqs, const NupackInvariants & invars) {
  auto spec_size = specifications.size();
  NUPACK_CHECK(s < spec_size, "Invalid specification index "
      + to_string(s) + " in decompose_ppair");
  NUPACK_CHECK(spec_size == res.results.size(),
      "Invalid result size, doesn't match spec size " + to_string(spec_size)
      + " != " + to_string(res.results.size()));

  specifications[s].decompose_ppair(i, k, res.results[s], seqs, invars);
}

#ifdef JSONCPP_FOUND
Json::Value PhysicalSpec::make_json_structures(const SequenceSpec & seqs,
    const NupackInvariants & invars) const {
  // Make json structures
  Json::Value root(Json::arrayValue);
  for (auto & spec : specifications) {
    const auto & strucs = spec.get_strucs();
    for (auto & struc : strucs) {
      Json::Value curstruc;
      const auto & names = struc.get_names();
      const auto & pairs = struc.get_structures();
      const auto & breaks = struc.get_breaks();
      
      auto n_size = names.size();
      auto p_size = pairs.size();
      NUPACK_DEBUG("names: " << n_size << " pairs: " << p_size);
      NUPACK_CHECK(n_size == p_size,
          "Name structure size mismatch: " + to_string(n_size) 
          + " != " + to_string(p_size));

      Json::Value strands(Json::arrayValue);
      auto strand_names = struc.get_strand_names();

      for (auto & sn : strand_names) strands.append(sn);

      for (auto j = 0; j < n_size; j++) {
        curstruc["name"] = names[j];
        curstruc["structure"] = StructureUtils::pairs_to_dpp(pairs[j], breaks);
        curstruc["strands"] = strands;
        root.append(curstruc);
      }
    }
  }
  return root;
}

Json::Value PhysicalSpec::make_json_tubes() const {
  Json::Value root(Json::arrayValue);
  for (auto & spec : specifications) {
    DBL_TYPE water_conc = water_density(spec.get_params().temperature - ZERO_C_IN_KELVIN);
    const auto & tubes = spec.get_tubes();
    const auto & ord_struc_map = spec.get_ord_struc_map();
    for (auto & t : tubes) {
      Json::Value curtube;

      curtube["name"] = t.get_name();
      const auto & complexes = t.get_complexes();

      Json::Value strucs(Json::arrayValue);
      
      for (auto & comp : complexes) {
        auto struc_ind = ord_struc_map[comp.order_ind];
        if (struc_ind >= 0) {
          auto & names = spec.get_strucs()[struc_ind].get_names();

          for (auto k = 0; k < comp.target_inds.size(); k++) {
            Json::Value curtarget;

            curtarget["name"] = names[comp.target_inds[k]];
            curtarget["concentration[M]"] = 
                static_cast<double>(water_conc * comp.target_concs[k]);
            strucs.append(curtarget);
          }
        }
      }
            
      root.append(curtube);
    }
  }
  return root;
}
#endif // JSONCPP_FOUND
}
