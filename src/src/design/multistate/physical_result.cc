
#include "physical_result.h"
#include "sequence_spec.h"

#include "design_debug.h"
#include "algorithms.h"

#include <sys/time.h>
#include <cfloat>
#include <sstream>
#include <algorithm>
#include <numeric>

namespace nupack {

void SingleParamResult::print_res_files(const SingleParamSpec & spec,
    const NupackInvariants & invars, std::string file_prefix) const {
  const auto & strucspecs = spec.get_strucs();
  for (auto i = 0; i < this->strucs.size(); i++) {
    strucs[i].print_strucfiles(strucspecs[i], this->params, invars, file_prefix);
  }
}

void SingleParamResult::serialize(const SingleParamSpec & spec,
    const NupackInvariants & invars, std::ostream & out,
    int indent, std::string prefix) const {
  const auto & params = spec.get_params();
  const auto & ss = spec.get_strucs();
  const auto & ts = spec.get_tubes();
  
  std::string ind_str(indent, ' ');
  std::string name_ind(indent, ' ');
  ind_str = prefix + ind_str;
  if (indent > 2) name_ind = prefix + std::string(indent-2, ' ') + "- ";

  out << name_ind << "name: " << spec.get_name() << std::endl;

  if (this->strucs.size() > 0) {
    out << ind_str << "structures: " << std::endl;
    
    // Serialize structures
    auto struc_it = this->strucs.cbegin();
    auto ss_it = ss.cbegin();
    for ( ; struc_it != this->strucs.end() && ss_it != ss.end(); ++struc_it, ++ss_it)
      struc_it->serialize(*ss_it, params, invars, out, indent + 4, prefix);
    
  } else {
    out << ind_str << "structures: []" << std::endl;
  }

  if (this->tubes.size() > 0) {
    out << ind_str << "tubes: " << std::endl;
    
    // Serialize tubes
    auto tube_it = tubes.cbegin();
    auto ts_it = ts.cbegin();
    for ( ; tube_it != this->tubes.end() && ts_it != ts.end(); ++tube_it, ++ts_it) {
      tube_it->serialize(*ts_it, ss, spec.get_ord_struc_map(), params, invars,
          out, indent + 4, prefix);
    }
  } else {
    out << ind_str << "tubes: []" << std::endl;
  }
}

void SingleParamResult::replace_node(const SingleParamResult & other, 
    const SingleParamSpec & spec, int i, int k, const NupackInvariants & invars) {
  NUPACK_CHECK(this->strucs.size() == other.strucs.size(), 
      "Structure sizes don't match in replace node");
  NUPACK_CHECK(this->tubes.size() == other.tubes.size(),
      "Tube sizes don't match in replace node");
  NUPACK_CHECK(i < this->strucs.size(), to_string(i)
      + " is not a valid structure index "+ "# strucs: " 
      + to_string(this->strucs.size()));

  this->strucs[i].replace_node(other.strucs[i], k, spec.get_params(), invars);

  const auto & order_to_struc = spec.get_ord_struc_map();
  const auto & tubes = spec.get_tubes();
  for (auto i_tube = 0; i_tube < this->tubes.size(); i_tube++) {
    NUPACK_DEBUG("Evaluating tube: " << to_string(i_tube));
    this->tubes[i_tube].evaluate(this->strucs, this->orders,
        order_to_struc, tubes[i_tube], spec.get_params(), invars);
  }
}

DBL_TYPE SingleParamResult::get_eval_time() const {
  DBL_TYPE tot_time = accumulate(orders, (DBL_TYPE)(0.0),
      [](const decltype(tot_time) & t, decltype(orders)::const_reference ord) {
        return ord.get_evaluated() ? t + ord.get_eval_time() : 0;
      });
  tot_time = accumulate(strucs, tot_time, 
      [](const decltype(tot_time) & t, decltype(strucs)::const_reference str) {
        return t + str.get_eval_time(); 
      });

  return tot_time;
}

void SingleParamResult::evaluate(const SequenceState & seqs, 
    const SingleParamSpec & spec, const NupackInvariants & invars) {
  const auto & strucs = spec.get_strucs();
  const auto & tubes = spec.get_tubes();
  const auto & orders = spec.get_orders();
  const auto & order_to_struc = spec.get_ord_struc_map();

  int n_strucs = strucs.size();
  int n_orders = orders.size();
  int n_tubes = tubes.size();

  this->params = spec.get_params();
  this->strucs.resize(n_strucs);
  this->orders.resize(n_orders);
  this->tubes.resize(n_tubes);

  for (auto i_struc = 0; i_struc < n_strucs; i_struc++) {
    this->strucs[i_struc].evaluate(seqs, strucs[i_struc], spec.get_params(), invars);
  }

  bool eval_off_targets = spec.eval_off_targets();
  for (auto i_ord = 0; i_ord < n_orders; i_ord++) {
    if (eval_off_targets && -1 == order_to_struc[i_ord]) {
      this->orders[i_ord].evaluate(seqs, orders[i_ord], spec.get_params(), invars);
    } else {
      this->orders[i_ord].clear_evaluated();
    }
  }

  for (auto i_tube = 0; i_tube < n_tubes; i_tube++) {
    this->tubes[i_tube].evaluate(this->strucs, this->orders,
        order_to_struc, tubes[i_tube], spec.get_params(), invars);
  }
}

int SingleParamResult::get_max_depth() const {
  if (strucs.size() == 0) return 0;
  using str = typename decltype(strucs)::value_type;
  return std::max_element(strucs.begin(), strucs.end(), 
      smaller_depth<str>)->get_max_depth();
}

DBL_TYPE PhysicalResult::get_eval_time() const {
  using res = typename decltype(results)::value_type;
  return accumulate(results, (DBL_TYPE)(0.0), 
      [](const DBL_TYPE a, const res b) { return a + b.get_eval_time(); });
}

void PhysicalResult::evaluate(const SequenceState & seqs, const PhysicalSpec & spec,
    const NupackInvariants & invars) {
  const auto & specs = spec.get_param_specs();
  int n_specs = specs.size();
  results.resize(n_specs);
  for (auto i_spec = 0; i_spec < n_specs; i_spec++) {
    results[i_spec].evaluate(seqs, specs[i_spec], invars);
  }
}

int PhysicalResult::get_max_depth() const {
  if (results.size() == 0) return 0;
  using res = typename decltype(results)::value_type;
  return std::max_element(results.begin(), results.end(), 
      smaller_depth<res>)->get_max_depth();
}

void PhysicalResult::replace_node(const PhysicalResult & other, const PhysicalSpec & spec,
    int s, int i, int k, const NupackInvariants & invars) {
  const auto & specs = spec.get_param_specs();
  NUPACK_CHECK(results.size() == other.results.size(), 
      "Result sizes don't match in PhysicalResult::replace_node");
  NUPACK_CHECK(results.size() == specs.size(),
      "Result sizes don't match spec sizes in PhysicalResult::replace_node");
  NUPACK_CHECK(s < results.size(),
      to_string(s) + " is not a valid specification index");

  results[s].replace_node(other.results[s], specs[s], i, k, invars);
}

void PhysicalResult::serialize(const PhysicalSpec & spec, 
    NupackInvariants & invars, std::ostream & out, int indent, 
    std::string prefix) const {
  std::string ind_str(indent, ' ');
  ind_str = prefix + ind_str;

  DBL_TYPE root_time = this->get_eval_time();
  int n_on_targets = 0;
  int n_off_targets = 0;
  int n_orders = 0;
  
  const auto & sps = spec.get_param_specs();
  for (const auto & sp : sps) {
    for (const auto & struc : sp.get_strucs()) {
      if (struc.get_structures().size() > 0) n_on_targets++;
    }
    n_orders += sp.get_orders().size();
  }
  n_off_targets = n_orders - n_on_targets;

  out << ind_str << "design properties:" << std::endl;
  out << ind_str << "    evaluation time: " << flt_format << root_time << std::endl;
  out << ind_str << "    on-targets: " << n_on_targets << std::endl;
  out << ind_str << "    off-targets: " << n_off_targets << std::endl;
  out << ind_str << "physical results:" << std::endl;

  auto spr_it = results.begin(); 
  auto sps_it = sps.begin();
  for ( ; spr_it != this->results.end() && sps_it != sps.end(); ++spr_it, ++sps_it) {
    spr_it->serialize(*sps_it, invars, out, indent + 4, prefix);
  }
}

void PhysicalResult::print_res_files(const PhysicalSpec & spec,
    const NupackInvariants & invars, std::string file_prefix) const {
  const std::vector<SingleParamSpec> & param_specs = spec.get_param_specs();
  auto res_it = results.begin(); 
  auto sp_it = param_specs.begin();
  for ( ; res_it != results.end() && sp_it != param_specs.end(); ++res_it, ++sp_it) {
    res_it->print_res_files(*sp_it, invars, file_prefix);
  }
}

}
