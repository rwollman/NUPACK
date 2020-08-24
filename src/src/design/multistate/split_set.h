#pragma once

#include <thermo.h>

#include <iostream>
#include <vector>
#include <utility>

namespace nupack {
  class NupackInvariants;
  
  struct Split {
    Split(int l, int r) : left(l), right(r) {};
    
    int left;
    int right;
  };
  
  class SplitSet {
    public:
      SplitSet() {}
      ~SplitSet() {}
      SplitSet(const SplitSet & other, const std::vector<int> & map);

      void push_back(int i, int j, DBL_TYPE prob);
      bool crosses_all(int i, int j) const;

      void add_prob(DBL_TYPE prob) { this->prob += prob; }
      void set_prob(DBL_TYPE prob) { this->prob = prob; }
      DBL_TYPE get_prob() const { return this->prob; }

      const std::vector<Split> & get_points() const { return this->points; }

      DBL_TYPE get_cost(int n);

      static bool allowed(int i, int j, int n, const NupackInvariants & invars);

      void serialize(std::ostream & out);
    private:
      void clear();
      void clone(const SplitSet & other);

      std::vector<Split> points;
      DBL_TYPE prob {0};
  };
}
