#include "node_spec.h"
#include "node_result.h"

#include "physical_spec.h"
#include "split_set.h"
#include "pair_probabilities.h"

#include "structure_utils.h"
#include "pathway_utils.h"
#include "design_debug.h"

#include <iostream>
#include <algorithm>

namespace nupack {

void NodeSpec::clear() {
  lims.clear();
  assume.clear();
  children.clear();
}

NodeSpec NodeSpec::get_depth(int depth) const {
  NodeSpec rval(lims, assume);
  if (depth > 0) {
    for (const auto & c : children) {
      NodeSpec child = c.get_depth(depth - 1);
      rval.children.push_back(child);
    }
  }
  return rval;
}

void NodeSpec::print_leaves() const {
  if (children.empty()) {
    print_node_decomp();
  } else {
    for (const auto & c : children) c.print_leaves();
  }
}

void NodeSpec::print_decomposition() const {
  print_node_decomp();
  for (const auto & c : children) c.print_decomposition();
}

void NodeSpec::print_node_decomp(std::ostream & out) const {
  using std::max;
  int maxsize = 0;
  const auto & m = *std::max_element(lims.begin(), lims.end(), 
      [](const Limits & a, const Limits & b) {
      return max(a.first, a.second) < max(b.first, b.second); } );
  maxsize = max(m.first, m.second);
  
  std::string tmpstruc(maxsize + 1, ' ');
  for (auto & l : lims) {
    for (auto i_nuc = l.first; i_nuc < l.second; ++i_nuc) {
      tmpstruc[i_nuc] = '-';
    }
  }
  out << tmpstruc << std::endl;
}

// decompose a particular leaf node
void NodeSpec::decompose_ppair_at(int k, 
    const NodeResult & res, const SequenceState & seqs,
    StructureSpec & strucspec, const StructureResult & strucres,
    const PhysicalParams & params, const NupackInvariants & invars) {
  int c_k = k;
  if (res.get_n_children() == 0) {
    NUPACK_CHECK(k == 0, "Invalid k reached: " + to_string(k));
    decompose_ppair(res, seqs, strucspec, strucres, params, invars);
  } else {
    NUPACK_CHECK(res.get_n_children() == get_n_children(),
        "PPair decomp on invalid node " + to_string(res.get_n_children()) + " != " + 
        to_string(get_n_children()));
    int n_c = res.get_n_children();
    for (auto i_c = 0; i_c < n_c; i_c++) {
      auto cur_n_leaves = res.children[i_c].get_n_leaves();
      if (c_k >= cur_n_leaves) {
        c_k -= cur_n_leaves;
      } else {
        children[i_c].decompose_ppair_at(c_k, res.children[i_c],
            seqs, strucspec, strucres, params, invars);
        break;
      }
    }
  }
}

void NodeSpec::decompose_ppair(const NodeResult & res, const SequenceState & seqs,
    StructureSpec & strucspec, const StructureResult & strucres,
    const PhysicalParams & params, const NupackInvariants & invars) {
  forbid_children(strucspec, res, params, invars);

  auto strucs = get_structures(strucspec, invars);
  auto splits = StructureUtils::get_consistent_splits(strucs, invars);
  auto to_full = get_native_map(strucspec, invars);
  auto to_node = get_to_node(strucspec, invars, to_full);

  auto raw_ppairs = res.collect_pair_probs(params, invars);

  raw_ppairs.clear_forbidden(strucspec);
  PairProbs mapped_probs(raw_ppairs, to_node);

  auto minsplit = StructureUtils::get_minimal_splits(strucspec, splits, 
      mapped_probs, to_full.size(), invars);

  SplitSet native_minsplit(minsplit, to_full);

  split(native_minsplit);
  for (auto & c : children) {
    NodeResult tempres;
    tempres.init(c);
    tempres.evaluate_leaf(c, seqs, strucspec, strucres, params, invars);
    c.decompose_ppair(tempres, seqs, strucspec, strucres, params, invars);
  }
}

void NodeSpec::decompose(const StructureSpec & strucspec, const NupackInvariants & invars) {
  clear_children();

  auto strucs = get_structures(strucspec, invars);
  auto splits = StructureUtils::get_consistent_splits(strucs, invars);
  auto to_full = get_native_map(strucspec, invars);
  auto to_node = get_to_node(strucspec, invars, to_full);

  PairProbs raw_ppairs;
  double num_strucs = strucs.size();
  for (auto & struc : strucs) {
    PairProbs curprobs;
    for (auto i_nuc = 0; i_nuc < struc.size(); i_nuc++) {
      if (struc[i_nuc] >= 0) {
        curprobs.push_back(to_full[i_nuc], to_full[struc[i_nuc]], 
            1.0 / num_strucs);
      } else {
        curprobs.push_back(to_full[i_nuc], -1, 1.0 / num_strucs);
      }
    }
    raw_ppairs.merge(curprobs, 1.0, 1.0);
  }

  raw_ppairs.clear_forbidden(strucspec);
  PairProbs mapped_probs(raw_ppairs, to_node);

  auto minsplit = StructureUtils::get_minimal_splits(strucspec, splits, 
      mapped_probs, to_full.size(), invars);

  SplitSet native_minsplit(minsplit, to_full);

  split(native_minsplit);
  for (auto & c : children) c.decompose(strucspec, invars);
}

void NodeSpec::forbid_children(StructureSpec & strucspec,
    const NodeResult & res, const PhysicalParams & params,
    const NupackInvariants & invars) const {
  std::vector<Assume> assumed;
  if (children.size() > 0) {
    for (auto it = children.begin(); it != children.end(); it +=2) {
      auto & c = *it;
      for (auto & a : c.assume) {
        int c_i = a.first;
        int c_j = a.second;
        bool found = false;
        for (auto this_a : assume) {
          int p_i = this_a.first;
          int p_j = this_a.second;
          found = ((p_i == c_i && p_j == c_j) || (p_i == c_j && p_j == c_i));
        }
        if (!found) assumed.emplace_back(c_i, c_j);
      }
    }

    DBL_TYPE total_prob = 0.0;
    PairProbs ppairs = res.collect_pair_probs(params, invars);
    for (auto & a : assumed) {  
      std::vector<int> h_i;
      std::vector<int> h_j;
      for (auto j = - invars.H_split + 1; j < invars.H_split + 1; j++) {
        h_i.push_back(a.first + j);
        h_j.push_back(a.second - j);
      }
      std::vector<DBL_TYPE> h_ppairs = ppairs.get_pair_probs(h_i, h_j);
      DBL_TYPE min_ppair = *std::min_element(h_ppairs.begin(), h_ppairs.end());
      total_prob += min_ppair;
    }
    if (total_prob > invars.f_split) {
      NUPACK_DEBUG("Probabilities failed. Forbidding splits.");
      for (auto & a : assumed) strucspec.add_forbidden(a.first, a.second); 
    } else {
      NUPACK_DEBUG("Probabilities not maintained " << total_prob);
    }
  }
}

int NodeSpec::get_n_nodes() const {
  int n_nodes = 1; // For myself
  for (const auto & c : children) n_nodes += c.get_n_nodes();
  return n_nodes;
}

int NodeSpec::get_n_leaves() const {
  int n_leaves = 0;
  if (children.empty()) return 1;
  for (const auto & c : children) n_leaves += c.get_n_leaves();
  return n_leaves;
}

int NodeSpec::get_max_depth() const {
  if (children.empty()) return 0;
  return std::max_element(children.begin(), children.end(), 
      smaller_depth<NodeSpec>)->get_max_depth() + 1;
}

const NodeSpec & NodeSpec::get_child(int i) const {
  NUPACK_CHECK(i < children.size() && i >= 0, 
      "Invalid child " + to_string(i) + " number of children: " + to_string(i)); 
  return children[i];
}

void NodeSpec::split(const SplitSet & splits) {
  const auto & points = splits.get_points();
  clear_children();
  for (const auto & point : points) {
    int i_nuc = point.left;
    int j_nuc = point.right;
    std::vector<Limits> left_lims;
    std::vector<Limits> right_lims;

    int seg_start, seg_stop;
    unsigned int i_seg;
    for (i_seg = 0; i_seg < lims.size(); ++i_seg) {
      seg_start = lims[i_seg].first;
      seg_stop = lims[i_seg].second;
      if (seg_stop > i_nuc) {
        break;
      } else {
        left_lims.emplace_back(seg_start, seg_stop);
      }
    }

    NUPACK_CHECK(seg_start < i_nuc, "Invalid start/end " + to_string((int)i_seg) + 
        " ss " + to_string(seg_start) + " " + to_string(seg_stop));
    left_lims.emplace_back(seg_start, i_nuc);

    if (seg_stop <= j_nuc) {
      right_lims.emplace_back(i_nuc, seg_stop);
      i_seg++;
    } else {
      right_lims.emplace_back(i_nuc, j_nuc + 1);
    }

    for (; i_seg < lims.size(); i_seg++) {
      seg_start = lims[i_seg].first;
      seg_stop = lims[i_seg].second;
      if (seg_stop > j_nuc) {
        break;
      } else {
        right_lims.emplace_back(seg_start, seg_stop);
      }
    }

    NUPACK_CHECK(seg_stop > j_nuc + 1, "Invalid right split point");
    if (seg_start > i_nuc) {
      right_lims.emplace_back(seg_start, j_nuc + 1);
    }
    left_lims.emplace_back(j_nuc + 1, seg_stop);

    i_seg++;
    for (; i_seg < lims.size(); i_seg++) {
      seg_start = lims[i_seg].first;
      seg_stop = lims[i_seg].second;
      left_lims.emplace_back(seg_start, seg_stop);
    }
    
    std::vector<Assume> left_assume;
    std::vector<Assume> right_assume;
    for (auto & a : assume) {
      auto a_i = a.first;
      auto a_j = a.second;
      if (a_i < i_nuc || a_i > j_nuc) {
        NUPACK_CHECK(a_j < i_nuc || a_j > j_nuc, "Invalid split, crosses previous assumption");
        left_assume.emplace_back(a_i, a_j);
      } else {
        NUPACK_CHECK(a_j >= i_nuc && a_j <= j_nuc, "Invalid split, crosses previous assumption");
        right_assume.emplace_back(a_i, a_j);
      }
    }
    left_assume.emplace_back(i_nuc - 1, j_nuc + 1);
    right_assume.emplace_back(j_nuc, i_nuc);
    
    children.emplace_back(left_lims, left_assume);
    children.emplace_back(right_lims, right_assume);
  }
}

std::vector<int> NodeSpec::get_nuc_ids(const StructureSpec & struc, 
    const NupackInvariants & invars) const {
  const std::vector<int> & nucs = struc.get_nuc_ids();
  auto native_map = get_native_map(struc, invars);
  
  std::vector<int> final_nucs;
  final_nucs.reserve(native_map.size());
  for (auto n : native_map) final_nucs.emplace_back(nucs[n]);
  return final_nucs;
}

std::vector<int> NodeSpec::get_breaks(const StructureSpec & struc,
    const NupackInvariants & invars) const {
  auto full_breaks = struc.get_breaks();
  auto to_node = get_to_node(struc, invars);
  std::vector<int> breaks;

  for (auto i = 1; i < lims.size(); i++) {
    auto start = lims[i].first;
    NUPACK_DEBUG_CHECK(start >= invars.H_split,
        "Segment starts too close to beginning of structure: " 
        + to_string(i) + " " + to_string(start));
    if (invars.include_dummies) {
      breaks.push_back(to_node[start] - invars.H_split);
    } else {
      breaks.push_back(to_node[start]);
    }
  }

  for (auto br : full_breaks) {
    if (to_node[br] > 0) breaks.push_back(to_node[br]);
  }
  std::sort(breaks.begin(), breaks.end());

  return breaks;
}

std::vector<int> NodeSpec::get_native_map(const StructureSpec & struc,
    const NupackInvariants & invars) const {
  std::vector<int> native_map;
  for (auto & l : lims) {
    int seg_start = l.first;
    int seg_stop = l.second;
    if (invars.include_dummies) {
      if (seg_start >= invars.H_split) seg_start -= invars.H_split;
      if (seg_stop < struc.size() - invars.H_split) seg_stop += invars.H_split;
    }
    for (auto i_nuc = seg_start; i_nuc < seg_stop; i_nuc++) {
      native_map.push_back(i_nuc);
    }
  }
  return native_map;
}

std::vector<std::vector<int> > NodeSpec::get_structures(const StructureSpec & struc,
    const NupackInvariants & invars) const {
  auto to_full = get_native_map(struc, invars);
  auto to_node = get_to_node(struc, invars, to_full);
  const auto & pairings = struc.get_structures();
  
  std::vector<vec_structure> res;

  const int n = to_full.size();
  for (auto & p : pairings) {
    auto match = std::all_of(assume.begin(), assume.end(), [&p](const decltype(assume)::value_type & a) {
        return p.at(a.first) == a.second;});
    
    if (match) {
      vec_structure curstruc(n, -5);
      
      // what is this?
      int i;
      int j = -1;
      for (i = 0; i < n; i++) {
        int m_j = p[to_full[i]];
        
        j = -1;
        if (m_j >= 0) {
          j = to_node[m_j];
          curstruc[j] = i;
        }
        curstruc[i] = j;
      }

      NUPACK_CHECK(j >= -1, "Invalid nucleotide mapping " + to_string(i));
      res.push_back(curstruc);
    }
  }
  return res;
}

// reverses get_native_map
std::vector<int> NodeSpec::get_to_node(const StructureSpec & struc,
    const NupackInvariants & invars, std::vector<int> to_full) const {
  std::vector<int> to_node(struc.size(), 0);
  if (to_full.empty()) to_full = get_native_map(struc, invars);
  for (auto i = 0; i < to_full.size(); i++) to_node[to_full[i]] = i;
  return to_node;
}

// only relevant if include_dummies is used
std::vector<bool> NodeSpec::get_native(const StructureSpec & struc,
    const NupackInvariants & invars) const {
  std::vector<bool> native(get_native_map(struc, invars).size(), true);

  if (invars.include_dummies) {
    auto to_node = get_to_node(struc, invars);
    for (auto & l : lims) {
      auto seg_start = l.first;
      auto seg_stop = l.second;
      if (seg_start >= invars.H_split) {
        for (auto i_nuc = seg_start - invars.H_split; i_nuc < seg_start; i_nuc++) {
          native[to_node[i_nuc]] = false;
        }
      }
      if (seg_stop < struc.size() - invars.H_split) {
        for (auto i_nuc = seg_stop; i_nuc < seg_stop + invars.H_split; i_nuc++) {
          native[to_node[i_nuc]] = false;
        }
      }
    }
  }
  return native;
}
}
