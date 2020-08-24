#include "sequence_spec.h"
#include "sequence_utils.h"

#include "design_debug.h"
#include "nupack_invariants.h"
#include "algorithms.h"

namespace nupack {
void SingleSequenceSpec::serialize(const NupackInvariants & invars, std::ostream &out) const{
  int n_nucs = this->nuc_ids.size();
  out << this->name << ": "; 
  for (auto i = 0; i < n_nucs - 1; i++) out << this->nuc_ids[i] << ", ";
  out << this->nuc_ids[n_nucs - 1];
}

void SequenceSpec::add_strand_by_nuc_ids(std::string name, std::vector<int> nucs,
    std::vector<int> domain_ids) {
  SingleSequenceSpec tempspec(name, nucs);
  tempspec.set_domain_ids(domain_ids);

  this->strand_name_to_ind_map[name] = this->strands.size();
  this->strands.push_back(tempspec);
}

std::vector<int> SequenceSpec::get_strand_ids(std::vector<std::string> names) {
  std::vector<int> out;
  for (auto n : names) out.push_back(get_strand_id(n)); 
  return out;
}

int SequenceSpec::get_strand_id(std::string name) {
  for (auto i_str = 0; i_str < strands.size(); i_str++) {
    if (name == strands[i_str].get_name()) return i_str;
  }
  NUPACK_ERROR("Invalid strand name " + name);
}

void SequenceSpec::add_strand_by_domains(std::string name, std::vector<std::string> dom_names) {
  unsigned int n_doms = dom_names.size();
  unsigned int m_doms = this->domains.size();
  std::vector<int> ids;
  
  // build list of domains by ids
  for (auto i_dom = 0; i_dom < n_doms; i_dom++) {
    const std::string i_name = dom_names[i_dom];
    for (auto j_dom = 0; j_dom < m_doms; j_dom++) {
      const std::string j_name = this->domains[j_dom].get_name();
      if (i_name == j_name) {
        ids.push_back(j_dom);
        break;
      }
    }

    NUPACK_CHECK(ids.size() == i_dom + 1, "Domain name " + i_name + " in strand " 
        + name + " not found");
  }

  this->add_strand_by_domain_ids(name, ids);
}

void SequenceSpec::add_strand_by_domain_ids(std::string name, std::vector<int> dom_ids) {
  std::vector<int> nuc_ids;

  for (auto dom : dom_ids) {
    NUPACK_CHECK(dom >= 0 && dom < this->domains.size(),
        "Invalid domain index " + to_string(dom));

    const auto & ids = this->domains[dom].get_nuc_ids();
    for (auto id : ids) nuc_ids.push_back(id);
  }

  add_strand_by_nuc_ids(name, nuc_ids, dom_ids);
}

void SequenceSpec::set_variable(int id, int var) {
  NUPACK_CHECK(id < nucleotides.size(), "ID isn't valid in set_variable");
  NUPACK_CHECK(nucleotides.size() == variables.size(),
      "nucleotides and variables out of sync in set_variable");
  variables[id] = var;
}

const SingleSequenceSpec & SequenceSpec::get_element(std::string name) const {
  auto it = domain_name_to_ind_map.find(name);
  if (contains_it(it, domain_name_to_ind_map)) return domains.at(it->second);
  
  it = strand_name_to_ind_map.find(name);
  if (contains_it(it, strand_name_to_ind_map)) return strands.at(it->second);
  
  NUPACK_ERROR("Invalid element name: " + name);
}

const SingleSequenceSpec & SequenceSpec::get_domain(std::string name) const {
  auto it = domain_name_to_ind_map.find(name);
  NUPACK_CHECK(contains_it(it, domain_name_to_ind_map),
      "Invalid domain supplied " + name);
  
  // Should check/throw exception as necessary
  return domains.at(it->second);
}

const SingleSequenceSpec & SequenceSpec::get_strand(std::string name) const {
  auto it = strand_name_to_ind_map.find(name);
  NUPACK_CHECK(contains_it(it, strand_name_to_ind_map),
      "Strand " + name + " not found / created");
  
  // Should check/throw exception as necessary
  return strands.at(it->second);
}

void SequenceSpec::add_domain_by_nuc_ids(std::string name, std::vector<int> nucs) {
  domain_name_to_ind_map[name] = domains.size();
  domains.emplace_back(name, nucs);
}

void SequenceSpec::serialize_ids(const NupackInvariants & invars, std::ostream & out, 
    int indent, std::string prefix) const {
  for (const auto dom : domains) dom.serialize(invars, out);
  for (const auto str : strands) str.serialize(invars, out);
}

void SequenceSpec::serialize(const NupackInvariants & invars, std::ostream & out, 
    int indent, std::string prefix) const {
  std::vector<int> c_seq;
  std::string c_seq_str;

  std::string ind_str(indent, ' ');
  ind_str += prefix;
  auto n_doms = this->domains.size();
  auto n_strs = this->strands.size();
  
  out << ind_str;
  out << "N Domains: " << n_doms << std::endl;
  try {
    for (auto dom : domains) {
      const auto & name = dom.get_name();
      const auto & nuc_ids = dom.get_nuc_ids();

      c_seq.clear();
      for (auto n : nuc_ids) c_seq.push_back(nucleotides.at(n));
      c_seq_str = SequenceUtils::nuc_to_str(c_seq, invars.material);
      
      out << ind_str;
      // Align the domain and complements
      auto colon = (name.back() == '*') ? ": " : " : ";
      out << "    " << name << colon << c_seq_str << std::endl;
    }

    out << "N Strands: " << n_strs << std::endl;
    for (auto str : strands) {
      const auto & name = str.get_name();
      const auto & nuc_ids = str.get_nuc_ids();

      c_seq.clear();
      for (auto n : nuc_ids) c_seq.push_back(nucleotides.at(n));
      c_seq_str = SequenceUtils::nuc_to_str(c_seq, invars.material);

      out << ind_str;
      out << "    " << name << ": " << c_seq_str << std::endl;
    }
  } catch (NupackException & e) {
    e.print_message(std::cerr);
    NUPACK_ERROR("Error serializing");
  }
}

void SequenceSpec::add_domain_by_constraints(std::string name, std::string nuc_str,
    const NupackInvariants & invars) {
  add_domain_by_constraints(name, SequenceUtils::str_to_nuc(nuc_str), invars);
}

// TODO fix this method to be more clear and less duplicative
void SequenceSpec::add_domain_by_constraints(std::string name, std::vector<int> nuc_cons,
    const NupackInvariants & invars) {
  unsigned int n_nucs = nuc_cons.size();
  unsigned int c_max_id = nucleotides.size();
  std::string canonical_name;
  std::string comp_name;
  if (name.back() == '*') {
    canonical_name = std::string(name.begin(), name.end() - 1);
    comp_name = name;
  } else {
    canonical_name = name;
    comp_name = name + "*";
  }
  int cnamesize = canonical_name.size();

  std::vector<int> comp_cons;
  SequenceUtils::get_complement(nuc_cons, comp_cons, invars);

  bool found_comp = false;
  bool found_forw = false;
  for (auto & dom : domains) {
    auto & cur_name = dom.get_name();
    if (0 == canonical_name.compare(0, cnamesize, cur_name, 0, cnamesize)) {
      const auto & nuc_ids = dom.get_nuc_ids();
      if (comp_name == cur_name) {
        NUPACK_CHECK(n_nucs == nuc_ids.size(),
            "Domain and complement sizes unequal " + canonical_name);
        for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++)
          nucleotides[nuc_ids[i_nuc]] = comp_cons[i_nuc];
        
        found_comp = true;
        break;
      } else if (cur_name == canonical_name) {
        NUPACK_CHECK(n_nucs == nuc_ids.size(),
            "Domain and complement sizes unequal " + canonical_name);
        for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) 
          nucleotides[nuc_ids[i_nuc]] = nuc_cons[i_nuc];
        
        found_forw = true;
        break;
      }
    }
  }

  std::vector<int> nuc_ids(n_nucs, 0);
  if (!found_forw) {
    nucleotides.resize(c_max_id + n_nucs, 0);
    variables.resize(c_max_id + n_nucs, -1);
    
    for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      nuc_ids[i_nuc] = c_max_id + i_nuc;
      nucleotides[i_nuc + c_max_id] = nuc_cons[i_nuc];
    }
    add_domain_by_nuc_ids(name, nuc_ids);
  }

  if (!found_comp) {
    nucleotides.resize(c_max_id + 2 * n_nucs);
    variables.resize(c_max_id + 2 * n_nucs, -1);

    for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      nuc_ids[i_nuc] = c_max_id + n_nucs + i_nuc;
      nucleotides[c_max_id + n_nucs + i_nuc] = comp_cons[i_nuc];
    }
    add_domain_by_nuc_ids(comp_name, nuc_ids);
  }
}

std::string SequenceSpec::partially_specified_domain() const {
  std::vector<int> domain_nucs;
  for (const auto & dom : domains) {
    auto & name = dom.get_name();
    
    if (name.find('*') == std::string::npos) {
      // auto start = domain_name_to_ind_map[name];
      const auto & nuc_ids = dom.get_nuc_ids();
      auto cur_len = nuc_ids.size();
      domain_nucs.resize(cur_len);
      for (auto i = 0; i < cur_len; ++i) {
        domain_nucs[i] = nucleotides[nuc_ids[i]];
      }
      if (!SequenceUtils::all_are_nucleotides(domain_nucs)) return name; 
    }
  }
  return std::string("");
}

#ifdef JSONCPP_FOUND
Json::Value SequenceSpec::make_json_domains(
    const std::vector<int> & poss, 
    const NupackInvariants & invars) const {
  Json::Value root(Json::arrayValue);
  unsigned int i, j, n;

  n = this->domains.size();
  for (i = 0; i < n; i++) {
    Json::Value d;
    d["name"] = this->domains[i].get_name();

    std::vector<int> nuc_ids = this->domains[i].get_nuc_ids();
    std::vector<int> new_nucs;
    for (j = 0; j < nuc_ids.size(); j++) {
      new_nucs.push_back(poss[nuc_ids[j]]);
    }

    std::string tmp = SequenceUtils::nuc_to_str(new_nucs, invars.material);
    d["sequence"] = tmp;

    root.append(d);
  }
  return root;
}

Json::Value SequenceSpec::make_json_strands(
    const std::vector<int> & poss,
    const NupackInvariants & invars) const {
  Json::Value root(Json::arrayValue);

  int i, n;
  unsigned int j;
  n = this->strands.size();
  for (i = 0; i < n; i++) {
    Json::Value s;
    s["name"] = this->strands[i].get_name();
    std::vector<int> dom_ids = this->strands[i].get_domain_ids();
    Json::Value doms(Json::arrayValue);
    for (j = 0; j < dom_ids.size(); j++) {
      int c_dom = dom_ids[j];
      doms.append(this->domains[c_dom].get_name());
    }
    s["domains"] = doms;

    std::vector<int> nuc_ids = this->strands[i].get_nuc_ids();
    std::vector<int> new_nucs;
    for (j = 0; j < nuc_ids.size(); j++) {
      new_nucs.push_back(poss[nuc_ids[j]]);
    }

    s["sequence"] = SequenceUtils::nuc_to_str(new_nucs, invars.material);

    root.append(s);
  }
  return root;
}
#endif
}
