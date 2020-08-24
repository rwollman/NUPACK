#include "pair_probabilities.h"

#include "pathway_utils.h"
#include "design_debug.h"

#include "algorithms.h"

#include <algorithm>

namespace nupack {
PairProbs::PairProbs(const PairProbs & other, const std::vector<int> & map) {
  for (auto & p : other.probs) {
    int d = (p.i >= 0 && p.i < map.size()) ? map.at(p.i) : -1;
    int e = (p.j >= 0 && p.j < map.size()) ? map.at(p.j) : -1;
    if (d >= 0 && (p.j == -1 || e >= 0)) probs.emplace_back(d, e, p.prob);
  }
}

PairProbs::PairProbs(const std::vector<int> & structure) {
  for (auto i = 0; i < structure.size(); ++i) {
    auto j = structure[i];
    if (j >= 0 || j == -1) push_back(i, j, 1);
  }
}

std::vector<DBL_TYPE> PairProbs::get_pair_probs(std::vector<int> i_v, 
    std::vector<int> j_v) const {
  NUPACK_CHECK(i_v.size() == j_v.size(), "Invalid indices provided");

  std::vector<DBL_TYPE> res(i_v.size(), 0.0);
  for (auto & p : probs) {
    for (auto j = 0; j < i_v.size(); ++j) {
      if (i_v[j] == p.i && j_v[j] == p.j) res[j] = p.prob;
    }
  }
  return res;
}

std::vector<DBL_TYPE> PairProbs::get_mat(int n) const {
  std::vector<DBL_TYPE> res(n * (n + 1), 0);
  for (auto & p : probs) {
    NUPACK_DEBUG_CHECK(p.i < n && p.j < n,
        "Invalid indices: (" + to_string(p.i) + ", " + to_string(p.j) 
        + ") out of seq length: " + to_string(n));
    auto j = (p.j < 0) ? n : p.j;
    res[p.i * (n + 1) + j] = p.prob;
  }
  return res;
}

int PairProbs::get_n() const {
  using el = typename decltype(probs)::value_type;
  if (probs.size() == 0) return 1;
  auto n = std::max_element(probs.begin(), probs.end(), [](const el & a, const el & b) {
    return a.max_index() < b.max_index();
  })->max_index();
  return n + 1;
}

void PairProbs::clear_forbidden(const StructureSpec & spec) {
  using el = typename decltype(probs)::value_type;
  std::remove_if(probs.begin(), probs.end(), [&](const el & a) {
    return spec.is_forbidden(a.i, a.j);
  });
}

std::vector<std::pair<int, int> > PairProbs::get_inds() const {
  std::vector<std::pair<int, int> > res;
  res.reserve(probs.size());
  for (auto & p : probs) res.emplace_back(p.i, p.j);
  return res;
}

void PairProbs::clear() {
  probs.clear();
}

void PairProbs::push_back(int i, int j, DBL_TYPE ppair) {
  probs.emplace_back(i, j, ppair);
}

void PairProbs::merge(const PairProbs & other, DBL_TYPE other_scale, DBL_TYPE this_scale) {
  NUPACK_CHECK(other_scale >= 0 && other_scale <= 1, 
      "Invalid other_scale " + to_string(other_scale));
  NUPACK_CHECK(this_scale >= 0 && this_scale <= 1, 
      "Invalid this_scale " + to_string(this_scale));
  
  int i_pos = 0;
  int j_pos = 0;
  int n = std::max(this->get_n() + 1, other.get_n() + 1);
  std::vector<PairProbTriple> ppt;
  // iterate through both vectors of ppair triples
  while (i_pos < probs.size() || j_pos < other.probs.size()) {
    // this
    int i = (i_pos < probs.size()) ? probs[i_pos].i : n;
    int j = (i_pos < probs.size()) ? probs[i_pos].j : n;
    
    // other
    int k = (j_pos < other.probs.size()) ? other.probs[j_pos].i : n;
    int l = (j_pos < other.probs.size()) ? other.probs[j_pos].j : n;
    
    int j_ord = (j + n) % n;
    int l_ord = (l + n) % n;
    
    if (i < k || (i == k && j_ord < l_ord)) {
      ppt.emplace_back(i, j, probs[i_pos].prob * this_scale);
      ++i_pos;
    } else if (k < i || ((i == k) && l_ord < j_ord)) {
      ppt.emplace_back(k, l, other.probs[j_pos].prob * other_scale);
      ++j_pos;
    } else {
      NUPACK_CHECK(i == k && j == l, "Invalid assumption on k");
      ppt.emplace_back(i, j, 
          probs[i_pos].prob * this_scale + other.probs[j_pos].prob * other_scale);
      ++i_pos;
      ++j_pos;
    }
  }
  probs = ppt;
}

std::vector<DBL_TYPE> PairProbs::get_nuc_defects(std::vector<int> pairing) const {
  std::vector<DBL_TYPE> res(pairing.size(), 1.0);
  for (auto & p : probs) {
    NUPACK_CHECK(p.i >= 0 && p.i < pairing.size(),
        "Invalid nucleotide for current pairing");

    if (pairing[p.i] == p.j || (p.j >= 0 && pairing[p.j] == p.i)) {
      res[p.i] = 1.0 - p.prob;
      if (p.j >= 0) res[p.j] = 1.0 - p.prob;
    }
  }
  return res;
}

std::vector<DBL_TYPE> PairProbs::get_nuc_defects(const PairProbs & target) const {
  std::vector<DBL_TYPE> defects(target.get_n(), 0.0);
  sort(probs);
  sort(target.probs);
  
  auto p = probs.begin();
  auto t = target.probs.begin();
  while (p != probs.end() && t != target.probs.end()) {
    if (*p < *t) {
      ++p;
    }
    else if (*p == *t) {
      defects[p->i] += std::max<DBL_TYPE>(*t - *p, 0);
      ++p, ++t;
    }
    else if (*t < *p) {
      defects[t->i] += t->prob;
      ++t;
    }
    else {
      NUPACK_ERROR("problem with generalized defect. " 
          + to_string(p->i + 1) + ", " + to_string(p->j + 1) 
          + ", " + to_string(t->i + 1) + ", " + to_string(t->j + 1));
    }
  }
  
  for (; t != target.probs.end(); ++t) defects[t->i] += t->prob;
  
  return defects;
}


void PairProbs::serialize(std::ostream & out, int n) const {
  for (auto & p : probs) {
    int j = (p.j == -1) ? n : p.j;
    out << p.i + 1 << "  " << j + 1 << "  " << flt_format << p.prob << std::endl; 
  }
}
}
