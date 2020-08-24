#pragma once

#include "named_spec.h"

#include <thermo.h>

#include <vector>
#include <string>
#include <map>

namespace nupack {
  class OrderSpec;
  class TubeComplex;
  
  using Complex = std::vector<std::string>;
  using NameList = std::vector<Complex>;
  using Order = std::vector<int>;
  
  class TubeSpec : public NamedSpec {
    public:
      TubeSpec() : NamedSpec("", "tube"), resolved(false), passive_frac(0.0), maxsize(0) {}
      TubeSpec(const TubeSpec & other) { this->clone(other); }
      TubeSpec & operator=(const TubeSpec & other);

      void clone(const TubeSpec & other);
      // Used in pre-resolution phase
      void add_target(std::string name, DBL_TYPE conc=0.0);
      void set_target_conc(std::string name, DBL_TYPE conc);
      void set_maxsize(int maxsize) { this->maxsize = maxsize; }

      // Used to enter the resolved phase
      void resolve_structures(const std::map<std::string, std::pair<int, int> > & map);

      void add_ordering(int order_ind);
      void enumerate_off_targets(std::map<OrderSpec, int> & order_map);
      void resolve_listed_offtargets(const std::map<std::string, int> & strands, 
          std::map<OrderSpec, int> & order_map);

      const std::vector<std::string> & get_struc_names() const { return this->struc_names; }
      const std::vector<TubeComplex> & get_complexes() const { return complexes; }

      void set_passive(DBL_TYPE passive) { this->passive_frac = passive; }
      DBL_TYPE get_passive() const { return this->passive_frac; }

      NameList & get_whitelist() { return whitelist; };
      NameList & get_blacklist() { return blacklist; };
    private:
      OrderSpec make_order(const std::vector<int> & a);
      void update_orders(std::map<OrderSpec, int> & order_map, const OrderSpec & ord);
      
      void generate_strand_map(const std::map<OrderSpec, int> & order_map);
      void remove_blacklist_complexes(std::map<OrderSpec, int> & order_map, std::vector<Order> & blacklist);
      void add_whitelist_complexes(std::map<OrderSpec, int> & order_map, std::vector<Order> & whitelist);
      
      std::vector<Order> resolve_list(const std::map<std::string, int> & strands, NameList & list);
      
      // Used in pre-resolution phase
      std::vector<std::string> struc_names;
      NameList whitelist;
      NameList blacklist;
      
      bool resolved;

      // Used in post-resolution phase
      std::vector<int> strand_map;
      std::vector<TubeComplex> complexes;
      DBL_TYPE passive_frac; // Passive mass fraction for tree optimization

      int maxsize;
  };
  
  class TubeComplex {
    public:
      TubeComplex(int ind, const std::vector<int> & t_ind = {}, 
          const std::vector<DBL_TYPE> & t_concs = {}) :
          order_ind(ind), target_inds{t_ind}, target_concs{t_concs} {}
      
      bool is_off_target() { return target_concs.empty(); }
      
      int order_ind;
      std::vector<int> target_inds;
      std::vector<DBL_TYPE> target_concs;
  };
}
