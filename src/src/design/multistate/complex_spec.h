#pragma once 


#include <vector>


namespace nupack {
  
  class OrderSpec {
    public:
      friend bool operator< (const OrderSpec & o1, const OrderSpec & o2) {
        return o1.strands < o2.strands;
      }
      friend bool operator> (const OrderSpec & o1, const OrderSpec & o2) {
        return o1.strands > o2.strands;
      }
      friend bool operator==(const OrderSpec & o1, const OrderSpec & o2) {
        return o1.strands == o2.strands;
      }
      friend bool operator!= (const OrderSpec & o1, const OrderSpec & o2) {
        return o1.strands != o2.strands;
      }

      OrderSpec(const OrderSpec & o) : strands(o.strands), symmetry(o.symmetry) {}

      OrderSpec(const std::vector<int> & order);

      int get_symmetry() const { return this->symmetry; }
      const std::vector<int> & get_strands() const { return this->strands; }

    protected:
      std::vector<int> strands;

    private:
      int symmetry;
  };
}
