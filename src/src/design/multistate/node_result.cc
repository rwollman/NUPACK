#include "node_result.h"
#include "sequence_state.h"

#include "design_debug.h"
#include "algorithms.h"

#include <memory>

namespace nupack {

int NodeResult::get_max_depth() const {
  if (children.size() == 0) return 0;
  return (*std::max_element(children.begin(), children.end(), 
      smaller_depth<NodeResult>)).get_max_depth() + 1;
}

DBL_TYPE NodeResult::merge_pfuncs(const NodeResult & left,
    const NodeResult & right, const PhysicalParams & params, 
    const NupackInvariants & invars) {
  auto split_len = (invars.include_dummies) ? invars.H_split : 0;
  auto i = right.to_full[split_len];
  auto j = right.to_full[right.eval_sequence.size() - split_len - 1];
  auto d = i - 1;
  auto e = j + 1;
  auto a0 = left.eval_sequence[left.native_index(d)];
  auto a1 = right.eval_sequence[split_len];
  auto b0 = right.eval_sequence[right.eval_sequence.size() - split_len - 1];
  auto b1 = left.eval_sequence[left.native_index(e)];

  DBL_TYPE int_ene = HelixEnergy(a0, b1, a1, b0);
  return EXP_FUNC(-int_ene / (kB * params.temperature)) * left.get_pfunc() * right.get_pfunc();
}

void NodeResult::replace_node(const NodeResult & other, int k,
    const PhysicalParams & params, const NupackInvariants & invars) {
  int num_children = get_n_children();
  int other_num_children = other.get_n_children();
  int c_k = k;
  if (other_num_children == 0 || num_children == 0) {
    NUPACK_CHECK(k == 0, "Reached leaf with child index: " + to_string(other_num_children) + 
        ", " + to_string(num_children) + " k = " + to_string(k));
    NUPACK_DEBUG("Current pfunc: " << -kB * params.temperature * LOG_FUNC(pfunc_corrected));
    NUPACK_DEBUG("Other pfunc: " << -kB * params.temperature * LOG_FUNC(other.pfunc_corrected));
    *this = other;
  } else {
    NUPACK_CHECK(other_num_children == num_children, 
        "Alternate decomposition must be sourced from the same tree");
    for (auto i = 0; i < other_num_children; ++i) {
      auto child_n_leaves = other.children[i].get_n_leaves();
      child_n_leaves = std::min(child_n_leaves, children[i].get_n_leaves());
      if (c_k >= child_n_leaves) {
        c_k -= child_n_leaves;
      } else {
        children[i].replace_node(other.children[i], c_k, params, invars);
        c_k = -1;
        break;
      }
    }

    pfunc_corrected = 0;
    for (auto it = children.begin(); it != children.end(); it += 2) {
      pfunc_corrected += merge_pfuncs(*(it), *(it + 1), params, invars);
      NUPACK_CHECK(it != children.end(), "children incorrectly formatted.");
    }
    if (c_k < 0) {
      NUPACK_DEBUG("Merged 1: " << -kB * params.temperature * LOG_FUNC(pfunc_corrected))
    }
  }
}

int NodeResult::get_n_leaves() const {
  if (children.size() == 0) return 1;
  using child = typename decltype(children)::value_type;
  return accumulate(children, 0, [](const int a, const child & b) {
      return a + b.get_n_leaves();
  });
}

int NodeResult::native_index(int i) const {
  // Native indices are sorted by construction => binary search
  auto it = std::lower_bound(to_full.begin(), to_full.end(), i);
  NUPACK_CHECK(it != to_full.end(), "Error finding native index");
  return it - to_full.begin();
}

// TODO make NodeSpec::get_children() public
void NodeResult::init(const NodeSpec & spec) {
  clear();
  for (auto i = 0; i < spec.get_n_children(); ++i) {
    auto tmp = NodeResult();
    tmp.init(spec.get_child(i));
    children.emplace_back(tmp);
  }
}

DBL_TYPE NodeResult::evaluate_dummy_pfuncs(const NodeSpec & spec,
    const StructureSpec & struc, const PhysicalParams & params, 
    const NupackInvariants & invars) {
  int n_nucs = 2 * (invars.H_split + 1);
  const auto & assumed = spec.get_assume();
  const std::vector<int> & sequence = this->eval_sequence;
  std::vector<int> to_node = spec.get_to_node(struc, invars);

  // One for break, one for termination
  std::vector<int> curseq(n_nucs + 2, 0);
  DBL_TYPE res = 1;

  // TODO make the dummy sequences, should be from assumed_i to assumed_i + H_split then
  // strand break, then assumed_j - H_split + 1 through assumed_j evaluate the partition function
  // and take QB of the upper right hand corner
  std::vector<DBL_TYPE> qb(n_nucs * n_nucs, 0.0);
  EXTERN_QB = qb.data();

  if (invars.include_dummies) {
    for (auto & a : assumed) {
      auto j_nuc = 0;
      auto nuc_start = a.first;
      auto nuc_end = nuc_start + invars.H_split + 1;
      NUPACK_DEBUG_CHECK(nuc_end <= to_node.size(),
          "Invalid split point encountered");

      for (auto i_nuc = nuc_start; i_nuc < nuc_end; ++i_nuc) {
        curseq[j_nuc] = sequence[to_node[i_nuc]];
        ++j_nuc;
      }

      curseq[j_nuc] = STRAND_PLUS;
      ++j_nuc;

      nuc_end = a.second + 1;
      nuc_start = nuc_end - invars.H_split - 1;
      
      for (auto i_nuc = nuc_start; i_nuc < nuc_end; ++i_nuc) {
        curseq[j_nuc] = sequence[to_node[i_nuc]];
        ++j_nuc;
      }
      curseq[j_nuc] = -1;
      ++j_nuc;

      NUPACK_CHECK(j_nuc == (n_nucs + 2), "Invalid sequence length");

      pfuncFull(curseq.data(), 3, invars.material, invars.dangle_type, 
          params.temperature - ZERO_C_IN_KELVIN, 1, 
          invars.sodium, invars.magnesium,
          invars.use_long_helix);

      res *= qb[pf_index(0, 2*(invars.H_split + 1) - 1, 
          2*(invars.H_split + 1))];
    }
  } else {
    for (auto & a : assumed) {
      if (sequence[to_node[a.first]] != BASE_C && sequence[to_node[a.second]] != BASE_C) {
        res *= EXP_FUNC(-AT_PENALTY / (kB * params.temperature));
      }
    }
  }

  //needed because global garbage 
  EXTERN_QB = NULL;
  return res;
}


void NodeResult::evaluate_leaf(const NodeSpec & spec, const SequenceState& seqs,
    const StructureSpec & strucspec, const StructureResult & struc,
    const PhysicalParams & params,
    const NupackInvariants & invars) {
  DBL_TYPE bonus = invars.dG_clamp;
  DBL_TYPE start_time = get_current_time();
  int j_nuc;
  std::vector<int> to_node = spec.get_to_node(strucspec, invars);
  const auto & strucs = strucspec.get_structures();

  // Set the permanent info on the node
  this->nuc_ids = spec.get_nuc_ids(strucspec, invars);
  this->eval_sequence = seqs.get_sequence(this->nuc_ids);
  this->breaks = spec.get_breaks(strucspec, invars);
  this->to_full = spec.get_native_map(strucspec, invars);
  this->native = spec.get_native(strucspec, invars);

  int n_nucs = this->eval_sequence.size();
  NUPACK_CHECK((int)this->to_full.size() == n_nucs, 
      "to_full size does not match eval_sequence size");

  // Construct eval sequence with breaks
  std::vector<int> fullseq;
  auto i_nuc = 0;
  for (auto & br : this->breaks) {
    fullseq.insert(fullseq.end(), this->eval_sequence.begin() + 
        i_nuc, this->eval_sequence.begin() + br);
    fullseq.push_back(STRAND_PLUS);
    i_nuc = br;
  }
  int n_brs = this->breaks.size();
  if (n_brs > 0) {
    fullseq.insert(fullseq.end(), this->eval_sequence.begin() + 
        this->breaks[n_brs - 1], this->eval_sequence.end());
  } else {
    fullseq.insert(fullseq.end(), this->eval_sequence.begin(), 
        this->eval_sequence.end());
  }
  fullseq.push_back(-1);
  NUPACK_CHECK(fullseq.size() == n_brs + 1 + eval_sequence.size(),
      "Fullseq size: " + to_string(fullseq.size()) + " evalseq size: " 
      + to_string(eval_sequence.size()) + " n breaks: " + to_string(n_brs));

  std::vector<DBL_TYPE> pp((n_nucs + 1) * n_nucs, 0.0);
  std::vector<DBL_TYPE> q(n_nucs * n_nucs, 0.0);
  // global passing mechanism using these variable names
  pairPr = pp.data();
  EXTERN_Q = q.data();

  const auto & assumed = spec.get_assume();
  DBL_TYPE pfunc = 1;
  if (invars.include_dummies) {
    pfuncFull(fullseq.data(), 3, invars.material, invars.dangle_type,
      params.temperature - ZERO_C_IN_KELVIN, 1, invars.sodium,
      invars.magnesium, invars.use_long_helix);
    pfunc = q[pf_index(0, n_nucs - 1, n_nucs)];
  } else {
    std::vector<DBL_TYPE> bonuses(n_nucs * n_nucs, 1);
    DBL_TYPE bonus_per = bonus;
    DBL_TYPE total_bonus = 1.0;
    for (auto & a : assumed) {
      i_nuc = to_node[a.first];
      auto j_nuc = to_node[a.second];
      if (!(i_nuc < j_nuc)) std::swap(i_nuc, j_nuc); 
      bonuses[i_nuc * n_nucs + j_nuc] *= EXP_FUNC(-bonus_per / (kB * params.temperature));
      total_bonus *= EXP_FUNC(-bonus_per / (kB * params.temperature));
    }
    pfuncFullWithBonuses(fullseq.data(), 3, invars.material, 
        invars.dangle_type, params.temperature - ZERO_C_IN_KELVIN, 1, 1,
        invars.sodium, invars.magnesium, invars.use_long_helix, bonuses.data());
    pfunc = q[pf_index(0, n_nucs - 1, n_nucs)];
    pfunc /= total_bonus;
  }

  this->ppairs.clear();
  for (i_nuc = 0; i_nuc < n_nucs; ++i_nuc) {
    if (this->native[i_nuc]) {
      std::vector<bool> saved_inds(n_nucs + 1, false);
      for (auto i_str = 0; i_str < strucs.size(); ++i_str) {
        auto m_j_nuc = strucs[i_str][to_full[i_nuc]];
        if (m_j_nuc >= 0) {
          j_nuc = to_node[m_j_nuc];
          if (j_nuc >= 0) saved_inds[j_nuc] = true; 
        } else {
          saved_inds[n_nucs] = true;
        }
      }
      for (j_nuc = 0; j_nuc < n_nucs; ++j_nuc) {
        if (this->native[j_nuc]) {
          if (pp[i_nuc * (n_nucs + 1) + j_nuc] > invars.min_ppair 
              || saved_inds[j_nuc]) {
            this->ppairs.push_back(to_full[i_nuc], to_full[j_nuc], 
                pp[i_nuc * (n_nucs + 1) + j_nuc]);
          }
        }
      }
      if (pp[i_nuc * (n_nucs + 1) + n_nucs] > invars.min_ppair
          || saved_inds[n_nucs]) {
        this->ppairs.push_back(to_full[i_nuc], -1, 
            pp[i_nuc * (n_nucs + 1) + j_nuc]);
      }
    }
  }

  DBL_TYPE min_ppair = 1.0;
  for (auto & a : assumed) {
    auto i_nuc = to_node[a.first];
    auto j_nuc = to_node[a.second];
    auto ppair = pp[i_nuc * (n_nucs + 1) + j_nuc];
    min_ppair = std::min(min_ppair, ppair);
  }
  min_ppair = std::max<DBL_TYPE>(min_ppair, 0.0001);

  pfunc *= min_ppair;
  DBL_TYPE dummy_pfuncs = this->evaluate_dummy_pfuncs(spec, strucspec, params, invars);
  DBL_TYPE orig_pfunc = pfunc;
  if (dummy_pfuncs < 1.0 && invars.include_dummies) dummy_pfuncs = 1.0;
  pfunc = orig_pfunc / dummy_pfuncs;
  this->pfunc_corrected = pfunc;
  if (!(pfunc > 0)) {
    NUPACK_DEBUG("PFUNC: " << exp_format << orig_pfunc);
    NUPACK_DEBUG("DUMMY: " << exp_format << dummy_pfuncs );
    NUPACK_DEBUG("lPFUNC: " << flt_format << -LOG_FUNC(orig_pfunc) );
    NUPACK_DEBUG("lPFUNC: " << flt_format << -LOG_FUNC(pfunc) );
  }

  DBL_TYPE end_time = get_current_time();
  this->eval_time = end_time - start_time;

  //needed because global garbage 
  EXTERN_Q = NULL;
  pairPr = NULL;
}

void NodeResult::evaluate(const NodeSpec & spec, const SequenceState & seqs,
    const StructureResult & struc, 
    const StructureSpec & strucspec, const PhysicalParams & params, 
    const NupackInvariants & invars) {
  int n_children = this->children.size();
  if (n_children != spec.get_n_children()) {
    this->clear();
    this->init(spec);
    n_children = this->children.size();
  }
  
  if (n_children == 0) {
    std::vector<int> nuc_ids = spec.get_nuc_ids(strucspec, invars);
    std::vector<int> breaks = spec.get_breaks(strucspec, invars);
    std::vector<int> eval_seq = seqs.get_sequence(nuc_ids);
    if (eval_seq != this->eval_sequence || breaks != this->breaks ||
        nuc_ids != this->nuc_ids) {
      this->evaluate_leaf(spec, seqs, strucspec, struc, params, invars);
    }
  } else {
    this->eval_time = 0;
    for (auto i = 0; i < n_children; ++i) {
      this->children[i].evaluate(spec.get_child(i), seqs, struc,
          strucspec, params, invars);
    }
    this->pfunc_corrected = 0;
    for (auto i = 0; i < n_children; i+=2) {
      this->pfunc_corrected += merge_pfuncs((this->children[i]),
          (this->children[i + 1]), params, invars);
    }
  }

  this->to_full = spec.get_native_map(strucspec, invars);
  this->nuc_ids = spec.get_nuc_ids(strucspec, invars);
  this->breaks = spec.get_breaks(strucspec, invars);
  this->eval_sequence = seqs.get_sequence(nuc_ids);
}

void NodeResult::serialize(std::ostream & out, int indent, const std::string & prefix, int & id) const {
  std::string ind_str(indent, ' ');
  ind_str = prefix + ind_str;
  int n_tot = *std::max_element(to_full.begin(), to_full.end());
  std::string outstr(n_tot + 1, ' ');
  std::string ids("0123456789abcdefghijklmnopqrstuvwxyz");
  auto id_char = ids[id % ids.size()];
  for (auto tf : to_full) outstr[tf] = id_char;
  out << ind_str << outstr << std::endl;
  
  ++id;
  if (this->children.size() > 0) {
    for (auto & c : children) c.serialize(out, indent, prefix, id);
  } 
}

DBL_TYPE NodeResult::collect_eval_times() const {
  if (this->children.size() == 0) return this->eval_time;
   
  DBL_TYPE ret = 0.0;
  for (auto & c : children) ret += c.collect_eval_times();
  return ret;
}

PairProbs NodeResult::collect_pair_probs(const PhysicalParams & params,
    const NupackInvariants & invars) const {
  if (this->children.size() == 0) return this->ppairs;
  
  PairProbs ppairs_tot;
  DBL_TYPE pfunc = 0;
  for (auto i = 0; i < this->children.size(); i += 2) {
    DBL_TYPE cur_pfunc = merge_pfuncs((this->children[i]),
        (this->children[i + 1]), params, invars);
    if (cur_pfunc > 0) {
      PairProbs ppairs1 = this->children[i].collect_pair_probs(params, invars);
      PairProbs ppairs2 = this->children[i + 1].collect_pair_probs(params, invars);
      ppairs1.merge(ppairs2, 1.0, 1.0);
      ppairs_tot.merge(ppairs1, cur_pfunc / (cur_pfunc + pfunc), pfunc / (cur_pfunc + pfunc));
      pfunc += cur_pfunc;
    }
  }
  return ppairs_tot; 
}

void NodeResult::clear_children() { children.clear(); }

void NodeResult::clear() {
  clear_children();
  this->eval_sequence.clear();
  this->nuc_ids.clear();
  this->to_full.clear();
  this->breaks.clear();
  this->nuc_defects.clear();
  this->ppairs.clear();
  this->pfunc_corrected = 0;
}  

}
