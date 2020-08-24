#pragma once

#include "named_spec.h"
#include "algorithms.h"
#include "sequence_utils.h"

#ifdef JSONCPP_FOUND
#include <json/json.h>
#endif

#include <vector>
#include <map>
#include <string>
#include <iostream>

namespace nupack { 
  class NupackInvariants;
  
  class SingleSequenceSpec : public NamedSpec {
    public:
      SingleSequenceSpec(std::string name, const std::vector<int> & nuc_ids) :
          NamedSpec(name, "seqspec"), nuc_ids(nuc_ids) {}
      void serialize(const NupackInvariants & invars, std::ostream & out) const;
      const std::string & get_name() const { return this->name; }
      const std::vector<int> & get_nuc_ids() const { return this->nuc_ids; }
      int size() const { return this->nuc_ids.size(); }
      const std::vector<int> & get_domain_ids() const { return this->domain_ids; }
      void set_domain_ids(const std::vector<int> & ids) { this->domain_ids = ids; }
      
    private:
      std::vector<int> nuc_ids;
      std::vector<int> domain_ids;
  };

  // Contains the complete specification linking to underlying set of
  // nucleotide assignment variables
  class SequenceSpec {
    public:
      SequenceSpec() {}

      std::vector<int> get_strand_ids(std::vector<std::string> names);
      int get_strand_id(std::string names);

      void serialize(const NupackInvariants & invars, std::ostream & out, 
          int indent = 0, std::string prefix = "") const;
      void serialize_ids(const NupackInvariants & invars, std::ostream & out, 
          int indent = 0, std::string prefix = "") const;

      void add_domain_by_constraints(std::string name, std::string nuc_cons, 
          const NupackInvariants & invars);
      void add_domain_by_constraints(std::string name, std::vector<int> nuc_cons,
          const NupackInvariants & invars);
      void add_domain_by_nuc_ids(std::string name, std::vector<int> nucs);

      void add_strand_by_nuc_ids(std::string name, std::vector<int> nucs, std::vector<int> domain_ids = {});

      void add_strand_by_domains(std::string name, std::vector<std::string> names);
      void add_strand_by_domain_ids(std::string name, std::vector<int> dom_ids);

      void set_variable(int id, int var);

      const std::vector<SingleSequenceSpec> & get_domains() const { return this->domains; }
      const std::vector<SingleSequenceSpec> & get_strands() const { return this->strands; }
      const std::vector<int> & get_nucleotides() const { return this->nucleotides; }
      const std::vector<int> & get_variables() const { return this->variables; }

      const SingleSequenceSpec & get_domain(std::string name) const;
      const SingleSequenceSpec & get_strand(std::string name) const;
      const SingleSequenceSpec & get_element(std::string name) const;

      int get_specified_id(std::string name, const std::map<std::string, int> & map) const;
      int get_domain_id(std::string name) const {return get_specified_id(name, domain_name_to_ind_map); }
      int get_strand_id(std::string name) const {return get_specified_id(name, strand_name_to_ind_map); }

      const std::map<std::string, int> & get_strand_map() const { 
        return this->strand_name_to_ind_map;
      }
      const std::map<std::string, int> & get_domain_map() const { 
        return this->domain_name_to_ind_map;
      }
      
      std::string partially_specified_domain() const;

#ifdef JSONCPP_FOUND
      Json::Value make_json_strands(
          const std::vector<int> & nucs, const NupackInvariants & invars) const;
      Json::Value make_json_domains(
          const std::vector<int> & nucs, const NupackInvariants & invars) const;
#endif

    protected:
      std::vector<int> nucleotides;
      std::vector<int> variables;

      std::vector<SingleSequenceSpec> domains;
      std::vector<SingleSequenceSpec> strands;

      std::map<std::string, int> domain_name_to_ind_map;
      std::map<std::string, int> strand_name_to_ind_map;
  };
  
  inline int SequenceSpec::get_specified_id(std::string name, 
      const std::map<std::string, int> & map) const {
    auto it = map.find(name);
    if (contains_it(it, map)) return it->second;
    
    return -1;
  }
}
