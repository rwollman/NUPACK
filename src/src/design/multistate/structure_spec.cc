#include "structure_spec.h"
#include "sequence_spec.h"
#include "structure_result.h"

#include "structure_utils.h"
#include "pathway_utils.h"
#include "design_debug.h"

#include "algorithms.h"

#include <algorithm>

namespace nupack {

StructureSpec::StructureSpec(const std::string & name, const std::string & struc) :
    StructureSpec(name, StructureUtils::dpp_to_pairs(struc).first, 
        StructureUtils::dpp_to_pairs(struc).second, {}) {}

StructureSpec::StructureSpec(const std::string & name, 
    const vec_structure & struc, const std::vector<int> & breaks, 
    const std::vector<int> & strands) : NamedSpec(name, "structure"), 
    strucs{struc}, strands(strands), breaks(breaks), struc_names{name},
    struc_ids{id}, symmetry(OrderSpec(strands).get_symmetry()) {
  std::vector<Limits> lims;
  lims.emplace_back(0, struc.size());
  tree = std::make_shared<NodeSpec>(lims);
}

StructureSpec::StructureSpec(const StructureSpec & other) : tree(nullptr) {
  clone(other);
}

StructureSpec & StructureSpec::operator=(const StructureSpec & other) {
  clone(other);
  return *this;
}

void StructureSpec::decompose_ppair(int k, const StructureResult & res, 
    const SequenceState & seqs, const PhysicalParams & params, 
    const NupackInvariants & invars) {
  tree->decompose_ppair_at(k, res.tree, seqs, *this, res, params, invars);
#ifndef NDEBUG
  NUPACK_DEBUG("Decomposition of " + name);
  for (auto & str : strucs) 
    NUPACK_DEBUG(StructureUtils::pairs_to_dpp(str, breaks));
  NUPACK_DEBUG("Leaves");
  tree->print_leaves();
  NUPACK_DEBUG("Full decomp");
  tree->print_decomposition();
#endif
}

void StructureSpec::set_id(int id) {
    this->id = id;
    this->struc_ids[0] = id;
}

// right rotation
void StructureSpec::rotate() {
  if (breaks.size() > 0) {
    int n_nucs = size();
    int rot_size = breaks.back();

    NUPACK_DEBUG("Before1:" << strands.size());
    NUPACK_DEBUG("Before2:" << nuc_ids.size());
    NUPACK_DEBUG("Before3:" << strucs[0].size());
    
    // nuc_ids are absolute    
    if (nuc_ids.size() > 0) {
      std::rotate(nuc_ids.rbegin(), nuc_ids.rbegin() + rot_size, nuc_ids.rend());
    }
    
    // structure map specification is relative
    vec_structure tmp_struc(n_nucs, -1);
    for (auto & str : strucs) {
      for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        if (str[i_nuc] < 0) {
          tmp_struc[(i_nuc + rot_size) % n_nucs] = -1;
        } else {
          tmp_struc[(i_nuc + rot_size) % n_nucs] = (str[i_nuc] + rot_size) % n_nucs;
        }
      }
      str = tmp_struc;
    }

    // this actually assumes the strands are strands, not domains.
    // must have actually resolved that by this point
    NUPACK_CHECK(strands.size() == strand_names.size(),
        "Can only rotate after resolving strand names");
    
    // strand ids and names are absolute
    std::rotate(strands.rbegin(), strands.rbegin() + 1, strands.rend());
    std::rotate(strand_names.rbegin(), strand_names.rbegin() + 1, strand_names.rend());
    
    // breaks are relative to lengths
    std::vector<int> new_breaks;
    new_breaks.push_back(rot_size);
    for (auto i_br = 0; i_br < breaks.size() - 1; i_br++) {
      new_breaks.push_back(breaks[i_br] + rot_size);
    }
    breaks = new_breaks;

    forbidden.clear();
  }
}

void StructureSpec::merge(const StructureSpec & other) {
  NUPACK_CHECK(nuc_ids.size() == other.nuc_ids.size(),
      "Size mismatch when merging " + other.name + " to " + name);
  if (strands != other.strands) {
    for (auto i_r = 0; i_r < strands.size() && strands != other.strands; i_r++)
      rotate();
  }
  NUPACK_CHECK(strands == other.strands,
      "Strand mismatch when merging " + other.name + " to " + name);

  for (const auto t : struc_names)
    for (const auto o : other.struc_names)
      NUPACK_CHECK(t != o, "Structures " + t + " and " + o + " shouldn't be merged twice")
    
  append(strucs, other.strucs);
  append(struc_names, other.struc_names);
  append(struc_ids, other.struc_ids);
}

void StructureSpec::compute_target_matrix() {
  double num_strucs = strucs.size();
  for (auto & struc : strucs) {
    PairProbs curprobs(struc);
    target_matrix.merge(curprobs, 1.0 / num_strucs, 1.0);
  }
}

bool StructureSpec::is_forbidden(int i, int j) const {
  for (const auto f : forbidden) {
    if ((i == f.first && j == f.second) || (i == f.second && j == f.first))
        return true;
  }
  for (const auto b : breaks) {
    if (i == b || i == b - 1 || j == b || j == b - 1) return true;
  }
  return false;
}

void StructureSpec::set_strands(const std::vector<int> & strand_ids) {
  OrderSpec temporder(strands);

  NUPACK_CHECK(strands.size() == strand_ids.size(),
      "Invalid: setting " + to_string(strand_ids.size()) + " to structure with "
      + to_string(strands.size()) + " strands");

  strands = strand_ids;
  symmetry = temporder.get_symmetry();
}

void StructureSpec::resolve_strand_names(const std::map<std::string, int> & name_map) {
  this->strands.clear();
  for (auto & s : strand_names) {
    NUPACK_CHECK(contains(s, name_map), "Unmapped strand name " + s);
    this->strands.push_back(name_map.find(s)->second);
  }
}

/**
 * Resolve the nucleotide ids based on the provided strand specification
 *
 * This sets the unique nucleotide ids, the domain ids, and the strand
 * ids for each structure.
 */
void StructureSpec::resolve_nuc_ids(const SequenceSpec & seqs) {
  const auto & strands = seqs.get_strands();
  const auto & domains = seqs.get_domains();
  nuc_ids.clear();
  domain_map.clear();
  strand_map.clear();
  breaks.clear();

  /* Resolve the actual nucleotide ids */
  int i_nuc = 0;
  for (auto & c_str : this->strands) {
    auto & s_nuc_ids = strands[c_str].get_nuc_ids();
    append(nuc_ids, s_nuc_ids);
    i_nuc += strands[c_str].size();
    breaks.push_back(i_nuc);
  }
  breaks.pop_back();


  for (auto & c_str : this->strands) {
    strand_map.insert(strand_map.end(), strands[c_str].get_nuc_ids().size(), c_str);
  }

  for (auto & c_str : this->strands) {
    const auto & cur_doms = strands[c_str].get_domain_ids();

    for (auto c_dom : cur_doms) {
      domain_map.insert(domain_map.end(), domains[c_dom].get_nuc_ids().size(), c_dom);
    }
  }

  NUPACK_CHECK(domain_map.size() == strand_map.size(),
          "Domain map size: " + to_string(domain_map.size())
          + " Strand map size: " + to_string(strand_map.size()));
  NUPACK_CHECK(domain_map.size() == i_nuc,
          "Nuc count: " + to_string(i_nuc)
          + " Strand map size: " + to_string(strand_map.size()));
}

StructureSpec StructureSpec::get_depth(int depth) const {
  StructureSpec res(*this);
  res.tree = std::make_shared<NodeSpec>(this->tree->get_depth(depth));
  return res;
}

OrderSpec StructureSpec::make_order() const {
  OrderSpec temp(this->strands);
  return temp;
}

void StructureSpec::clone(const StructureSpec & other) {
  if (this != &other) {

    this->NamedSpec::clone(other);
    this->strand_names = other.strand_names;
    this->strucs = other.strucs;
    this->target_matrix = other.target_matrix;
    this->struc_names = other.struc_names;
    this->struc_ids = other.struc_ids;
    this->domain_map = other.domain_map;
    this->strand_map = other.strand_map;
    this->strands = other.strands;
    this->nuc_ids = other.nuc_ids;
    this->breaks = other.breaks;
    this->forbidden = other.forbidden;
    this->symmetry = other.symmetry;
    this->tree = std::make_shared<NodeSpec>(*(other.tree));
  }
}

int StructureSpec::size() const {
  if (!nuc_ids.empty()) {
    return nuc_ids.size();
  } else if (!strucs.empty()) {
    return strucs[0].size();
  }
  NUPACK_ERROR("Size called on empty structure");
}


}
