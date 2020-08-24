#include "split_set.h"

#include "nupack_invariants.h"
#include "design_debug.h"

namespace nupack {
bool SplitSet::allowed(int i, int j, int n, const NupackInvariants & invars) {
  bool allow = true;
  int l1 = j - i + 1;
  int l2 = n - j + i - 1;
  if ((i < invars.H_split) || (j < i) || (j + invars.H_split >= n) ||
      (l1 < invars.N_split) || (l2 < invars.N_split)) {
    allow = false;
  }
  return allow;
}

SplitSet::SplitSet(const SplitSet & other, const std::vector<int> & map) {
  this->clear();
  try {
    for (auto & point : other.points) {
      points.emplace_back(map.at(point.left), map.at(point.right));
    }
  } catch (...) {
    this->clear();
    throw NupackException("Invalid map being used to create split set");
  }
  this->prob = other.prob;
}

void SplitSet::clear() {
  this->points.clear();
  this->prob = 0;
}

void SplitSet::serialize(std::ostream & out) {
  out << "SplitSet" << std::endl;
  for (auto & point : points) {
    out << point.left << " " << point.right << std::endl;
  }
  out << prob << std::endl;
}

void SplitSet::clone(const SplitSet & other) {
  this->points = other.points;
  this->prob = other.prob;
}

void SplitSet::push_back(int i, int j, DBL_TYPE prob) {
  bool found = false;
  for (auto & point : points) {
    if ((i == point.left && j == point.right) ||
        (i == point.right && j == point.left)) {
      found = true;
    }
  }
  if (!found) {
    this->prob += prob;
    points.emplace_back(i, j);
  }
}

DBL_TYPE SplitSet::get_cost(int n) {
  DBL_TYPE cost = 0;
  for (auto & point : points) {
    auto l1 = n - point.right + point.left - 1;
    auto l2 = point.right - point.left + 1;
    cost += pow((double) l1, 3) + pow((double) l2, 3);
  }
  return cost;
}

bool SplitSet::crosses_all(int i, int j) const {
  if (i > j) std::swap(i, j);
  for (auto & point : points) {
    auto d = point.left;
    auto e = point.right;
    if (d > e) std::swap(d, e);
    if (d == i && e == j) return false; // redundant
    if ((d > i && d < j) && (e > i && e < j)) return false; // inside
    if ((d < i || d > j) && (e < i || e > j)) return false; // outside
  }
  return true;
}
}
