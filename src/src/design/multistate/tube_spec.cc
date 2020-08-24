#include "tube_spec.h"
#include "structure_spec.h"

#include "pathway_utils.h"
#include "design_debug.h"
#include "algorithms.h"

#include <set>
#include <algorithm>

namespace nupack {
TubeSpec & TubeSpec::operator=(const TubeSpec & other) {
  if (this != &other) this->clone(other);
  return *this;
}

void TubeSpec::clone(const TubeSpec & other) {
  NamedSpec::clone(other);
  
  struc_names = other.struc_names;
  complexes = other.complexes;
  
  whitelist = other.whitelist;
  blacklist = other.blacklist;

  resolved = other.resolved;
  maxsize = other.maxsize;
  passive_frac = other.passive_frac;
}

void TubeSpec::add_target(std::string name, DBL_TYPE conc) {
  struc_names.push_back(name);
  complexes.emplace_back(TubeComplex{-1, {0}, {conc}});
}

void TubeSpec::set_target_conc(std::string name, DBL_TYPE conc) {
  NUPACK_CHECK(!this->resolved, "Trying to change resolved spec");
  
  auto it = std::find(struc_names.begin(), struc_names.end(), name);
  NUPACK_CHECK(contains_it(it, struc_names), 
      "Couldn't find structure " + name + " in tube " + this->name)

  complexes[it - struc_names.begin()].target_concs[0] = conc;
}

void TubeSpec::resolve_structures(const std::map<std::string, std::pair<int, int> > & map) {
  std::vector<TubeComplex> comps;

  std::map<int, int> order_map;
  for (auto i_struc = 0; i_struc < this->struc_names.size(); i_struc++) {
    const std::string & curname = this->struc_names[i_struc];
    auto cur_id = map.find(curname);

    NUPACK_CHECK(contains_it(cur_id, map), "Structure " + curname + " not given a mapping");

    int cur_order = cur_id->second.first;
    int cur_target = cur_id->second.second;

    auto cur_ind_it = order_map.find(cur_order);
    if (!contains_it(cur_ind_it, order_map)) {
      order_map[cur_order] = comps.size();
      
      comps.emplace_back(TubeComplex{cur_order, {cur_target}, 
          {complexes.at(i_struc).target_concs[0]}});
    } else {
      auto & comp = comps[cur_ind_it->second];
      comp.target_inds.emplace_back(cur_target);
      comp.target_concs.emplace_back(complexes[i_struc].target_inds[0]);
    }
  }
  
  complexes = comps;
  this->resolved = true;
}

void TubeSpec::add_ordering(int order_ind) {
  NUPACK_CHECK(this->resolved || this->complexes.empty(), "Conflict resolving structures");

  complexes.emplace_back(TubeComplex{order_ind});
}

OrderSpec TubeSpec::make_order(const std::vector<int> & a) {
  auto n = a.size();
  NUPACK_DEBUG("order size: " + to_string(n));
  std::vector<int> tempa(n, 0);
  for (auto k = 0; k < n; k++) {
    NUPACK_DEBUG("a[k]: " + to_string(a[k]));
    tempa[k] = strand_map.at(a[k]);
  }
  return {tempa};
}

// create OrderSpec based on some necklace specification
void TubeSpec::update_orders(std::map<OrderSpec, int> & order_map, const OrderSpec & ord) {
  if (!contains(ord, order_map)) {
    auto n_ord = order_map.size();
    order_map[ord] = n_ord;
    add_ordering(n_ord);
  } else {
    auto n_ord = order_map[ord];
    if (std::none_of(complexes.begin(), complexes.end(), [&] (const TubeComplex & c) {
        return c.order_ind == n_ord; })) add_ordering(n_ord);
  }
}

void TubeSpec::generate_strand_map(const std::map<OrderSpec, int> & order_map) {
  std::map<int, OrderSpec> revmap;
  for (auto & om : order_map) revmap.emplace(om.second, om.first);

  auto n_strucs = complexes.size();
  for (auto i_str = 0; i_str < n_strucs; i_str++) {
    auto it_rm = revmap.find(complexes[i_str].order_ind);
    NUPACK_CHECK(contains_it(it_rm, revmap), "On target generation failed. "
        + to_string(i_str) + " " + to_string(complexes[i_str].order_ind));
    const std::vector<int> & strands = it_rm->second.get_strands();
    for (auto & str : strands) {
      if (!lin_contains(str, strand_map)) strand_map.push_back(str);
    }
  }
  NUPACK_DEBUG("strand_map size: " + to_string(strand_map.size()));
}

void TubeSpec::enumerate_off_targets(std::map<OrderSpec, int> & order_map) {
  generate_strand_map(order_map);

  std::vector<int> a; 
  int n_strands = strand_map.size();

  NUPACK_DEBUG("N complexes before: " << order_map.size());
  
  for (auto n = 1; n <= this->maxsize; n++) {
    a.assign(n, 0);

    auto ord = make_order(a);
    update_orders(order_map, ord);
    int i = n;
    for ( ; i > 0; i--) {
      if (a[i-1] != n_strands - 1) break;
    }

    while (i > 0) {
      a[i-1] = a[i-1] + 1;
      for (auto j = 1; j < n - i + 1; j++) a[i + j - 1] = a[j - 1]; 

      if (n % i == 0) {
        NUPACK_DEBUG("n: " + to_string(n) + ", i: " + to_string(i));
        ord = make_order(a);
        update_orders(order_map, ord);
      }
      
      for (i = n; i > 0; i--) {
        if (a[i-1] != n_strands - 1) break;
      }
    }
  }
  NUPACK_DEBUG("N complexes after: " << order_map.size());
}

void TubeSpec::resolve_listed_offtargets(const std::map<std::string, int> & strands, 
    std::map<OrderSpec, int> & order_map) {
  if (!whitelist.empty()) {
    auto wlist = resolve_list(strands, whitelist);
    add_whitelist_complexes(order_map, wlist);
  }
  if (!blacklist.empty()) {
    auto blist = resolve_list(strands, blacklist);
    remove_blacklist_complexes(order_map, blist);
  }
}

std::vector<Order> TubeSpec::resolve_list(const std::map<std::string, int> & strands, NameList & list) {
  std::vector<Order> order_list;
  Order cur_order;
  
  for (auto el : list) {
    cur_order.clear();
    for (auto strand_name : el) {
      auto it = strands.find(strand_name);
      NUPACK_CHECK(lin_contains(it->second, strand_map), 
          strand_name + " contained in explicit list for tube " + name + 
          " is not a strand in an on-target complex for tube");
      cur_order.emplace_back(it->second);
    }
    
    order_list.emplace_back(cur_order);
  }
  
  return order_list;
}

void TubeSpec::remove_blacklist_complexes(std::map<OrderSpec, int> & order_map,
    std::vector<Order> & blacklist) {
  for (auto & b : blacklist) {
    auto it = order_map.find({b});
    // either error, as below, or ignore and continue (possibly alert via output)
    // if (!contains_it(b, order_map)) break;
    NUPACK_CHECK(contains_it(it, order_map),
        "Attempting to remove complex present in no ensembles");
    
    auto ord = it->second;
    auto comp = std::find_if(complexes.begin(), complexes.end(), [&] (const TubeComplex & c) {
        return c.order_ind == ord; });
    NUPACK_CHECK(contains_it(comp, complexes), 
        "Attempting to remove non-enumerated complex from tube " + name);
    NUPACK_CHECK(comp->is_off_target(), "Attempting to remove on-target complex");
    
    complexes.erase(comp);    
  }
}

void TubeSpec::add_whitelist_complexes(std::map<OrderSpec, int> & order_map,
    std::vector<Order> & whitelist) {
  for (auto & w : whitelist) update_orders(order_map, {w});
}
}
