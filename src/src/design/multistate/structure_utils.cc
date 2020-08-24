#include "structure_utils.h"
#include "split_set.h"
#include "pair_probabilities.h"

#include "nupack_invariants.h"
#include "design_debug.h"

#include "types.h"

#include <thermo.h>

namespace nupack { namespace StructureUtils {
structure_pair dpp_to_pairs(const std::string & struc) {
  std::vector<int> breaks;
  std::vector<int> pairs;
  std::vector<int> stack;
  
  int i_nuc = 0;
  for (char curchar : struc) {
    int j_nuc = -1;
    switch (curchar) {
      case '+': 
        breaks.push_back(i_nuc);
        break;
      case '.':
        pairs.push_back(-1);
        i_nuc++;
        break;
      case ')':
        NUPACK_CHECK(stack.size() > 0, "Unbalanced dot-parens");
        j_nuc = stack.back();
        pairs.push_back(j_nuc);
        pairs[j_nuc] = i_nuc;
        i_nuc++;
        stack.pop_back();
        break;
      case '(':
        pairs.push_back(-5);
        stack.push_back(i_nuc);
        i_nuc++;
        break;
      default:
        NUPACK_ERROR("Invalid character in dot-parens structure: " + struc);
    }
  }
  NUPACK_CHECK(stack.size() == 0, "Unbalanced dot-parens");
  return {pairs, breaks};
} 

std::string pairs_to_dpp(const vec_structure & pairs, 
    const std::vector<int> & breaks) {
  std::string struc;
  int n_nucs = pairs.size();
  int n_breaks = breaks.size();
  int i_break = 0;
  for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    if (i_break < n_breaks && breaks[i_break] == i_nuc) {
      struc.push_back('+');
      i_break++;
    }
    if (pairs[i_nuc] == -1) {
      struc.push_back('.');
    } else if (pairs[i_nuc] > i_nuc) {
      struc.push_back('(');
    } else if (pairs[i_nuc] < i_nuc) {
      struc.push_back(')');
    }
  }
  return struc;
}

void get_strand_lengths(const std::string & dpp, std::vector<int> & lengths) {
  std::vector<int> breaks;
  std::vector<int> pairs;
  std::tie(pairs, breaks) = dpp_to_pairs(dpp);
  get_strand_lengths(pairs, breaks, lengths);
}

void get_strand_lengths(const vec_structure & pairs,
    const std::vector<int> & breaks, std::vector<int> & lengths) {
  lengths.clear();
  int prev = 0;
  for (auto br : breaks) {
    lengths.push_back(br - prev);
    prev = br;
  }

  if (breaks.size() > 0) {
    lengths.push_back(pairs.size() - breaks.back());
  } else {
    lengths.push_back(pairs.size());
  }
}

std::vector<SplitSet> get_consistent_splits(
    const std::vector<vec_structure> & strucs, const NupackInvariants & invars) {
  std::vector<SplitSet> splits;
  splits.emplace_back();
  std::vector<SplitSet> newsplits;

  for (auto struc : strucs) {
    newsplits.clear();
    for (auto i_nuc = 0; i_nuc < struc.size(); i_nuc++) {
      auto j_nuc = struc.at(i_nuc);
      if (j_nuc > i_nuc) {
        for (auto s : splits) {
          if (s.crosses_all(i_nuc, j_nuc) && 
              SplitSet::allowed(i_nuc, j_nuc, struc.size(), invars)) {
            SplitSet newset = s;
            newset.push_back(i_nuc, j_nuc, 0);
            newsplits.push_back(newset);
          }
        }
      }
    }
    splits = newsplits;
  }
  return splits;
}

SplitSet get_minimal_splits(const StructureSpec & spec, 
    const std::vector<SplitSet> & set, const PairProbs & ppairs, int n,
    const NupackInvariants & invars) {
  auto ppair = ppairs.get_mat(n);
  auto tmp_inds = ppairs.get_inds();

  std::vector<DBL_TYPE> helix_probs(n * n, 0.0);
  std::vector<std::pair<int, int>> poss;
  std::vector<DBL_TYPE> poss_pp;

  DBL_TYPE cur_ppair;
  DBL_TYPE min_ppair;
  for (auto & t : tmp_inds) {
    auto i_nuc = t.first;
    auto j_nuc = t.second;
    if (i_nuc >= invars.H_split && i_nuc < n - invars.H_split
        && j_nuc >= invars.H_split && j_nuc < n - invars.H_split) {
      min_ppair = 1.0;
      for (auto i = -invars.H_split; i < invars.H_split; i++) {
        auto d_nuc = i_nuc + i;
        auto e_nuc = j_nuc - i;
        cur_ppair = ppair[d_nuc * (n + 1) + e_nuc];
        if (cur_ppair < min_ppair) min_ppair = cur_ppair;
      }
      helix_probs[i_nuc * n + j_nuc] = min_ppair;
    }
  }

  for (auto & t : tmp_inds) {
    auto i_nuc = t.first;
    auto j_nuc = t.second;
    if (i_nuc >= 0 && j_nuc >= 0) {
      if (helix_probs[i_nuc * n + j_nuc] > 0 && SplitSet::allowed(i_nuc, j_nuc, n, invars)) {
        poss.emplace_back(i_nuc, j_nuc);
        poss_pp.push_back(helix_probs[i_nuc * n + j_nuc]);
      }
    }
  }

  std::vector<SplitSet> pset = set;
  DBL_TYPE best_cost = DBL_MAX;
  SplitSet best_split;
  for (auto & ss : pset) {
    const auto & points = ss.get_points();
    ss.set_prob(0.0);
    for (const auto & point : points) {
      auto i_nuc = point.left;
      auto j_nuc = point.right;
      ss.add_prob(helix_probs[i_nuc * n + j_nuc]);
    }
    if (ss.get_prob() > invars.f_split && ss.get_cost(n) < best_cost) {
      best_cost = ss.get_cost(n);
      best_split = ss;
    }
  }

  std::vector<DBL_TYPE> new_costs(poss.size(), 0);
  std::vector<DBL_TYPE> split_prob(poss.size(), 0);

  for (auto i = 0; i < poss.size(); i++) {
    SplitSet tmp;
    tmp.push_back(poss[i].first, poss[i].second, 0);
    new_costs[i] = tmp.get_cost(n);
    split_prob[i] = helix_probs[poss[i].first * n + poss[i].second];
  }

  for (auto & ss : pset) {
    int i_stack = 0;

    std::vector<bool> split_allowed(poss.size(), false);

    for (auto i = 0; i < poss.size(); i++) {
      split_allowed[i] = ss.crosses_all(poss[i].first, poss[i].second)
        && SplitSet::allowed(poss[i].first, poss[i].second, n, invars);
    }

    std::vector<SplitSet> split_stack(1, ss);
    std::vector<DBL_TYPE> cost_stack(1, ss.get_cost(n));;
    std::vector<std::vector<bool> > split_allowed_stack;
    split_allowed_stack.push_back(split_allowed);
    int max_depth = 6;
    while (i_stack >= 0) {
      NUPACK_CHECK(i_stack == cost_stack.size() - 1, "Out of sync i_stack");
      DBL_TYPE max_ppair = 0;
      int j_pos = -1;
      for (auto i = 0; i < poss.size(); i++) {
        if (split_allowed_stack[i_stack][i] && max_ppair < split_prob[i] &&
            new_costs[i] + cost_stack[i_stack] < best_cost) {
          j_pos = i;
          max_ppair = split_prob[i];
        }
      }

      if (j_pos >= 0) {
        auto i_nuc = poss[j_pos].first;
        auto j_nuc = poss[j_pos].second;

        SplitSet tmpsplit = split_stack[i_stack];
        tmpsplit.push_back(i_nuc, j_nuc, max_ppair);
        split_allowed_stack[i_stack][j_pos] = false;

        if (tmpsplit.get_prob() > invars.f_split) {
          // Satisfy probability condition
          // TODO check that split isn't forbidden
          NUPACK_CHECK(tmpsplit.get_cost(n) < best_cost, 
              "Cost condition not satisfied on verification");
          best_cost = tmpsplit.get_cost(n);
          best_split = tmpsplit;
        } else if (i_stack < max_depth - 1) {
          split_stack.push_back(tmpsplit);
          cost_stack.push_back(tmpsplit.get_cost(n));
          split_allowed_stack.push_back(split_allowed_stack[i_stack]);
          for (auto i = 0; i < poss.size(); i++) {
            split_allowed_stack[i_stack+1][i] = split_allowed_stack[i_stack][i] &&
              tmpsplit.crosses_all(poss[i].first, poss[i].second);
          }
          
          i_stack++;
        } else {
          split_stack.pop_back();
          cost_stack.pop_back();
          split_allowed_stack.pop_back();

          i_stack --;
        }
      } else {
        split_stack.pop_back();
        cost_stack.pop_back();
        split_allowed_stack.pop_back();

        i_stack --;
      }
    }
  }

  return best_split;
}
}}
