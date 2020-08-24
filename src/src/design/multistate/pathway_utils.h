#pragma once

#include "physical_spec.h"
#include "sequence_spec.h"
#include "pair_probabilities.h"
#include "named_spec.h"

#include "types.h"

#include <thermo.h>

#include <vector>
#include <sstream>
#include <type_traits>
#ifdef HAVE_CSTDINT_H 
  #include <cstdint>
#endif

namespace nupack {
  class DesignSpec;

  int parse_design(std::string filename, DesignSpec & spec); 
  int sample_weighted_int(const std::vector<DBL_TYPE> & weights);
  DBL_TYPE get_current_time();

  class SourceSpec : public NamedSpec {
    public:
      SourceSpec(const std::string & name, const std::vector<int> & nucs) : 
          NamedSpec(name, "source"), 
          nucs(nucs) {}

      const std::vector<int> & get_nucs() const { return this->nucs; }

    private:
      PairProbs ppairs;
      std::vector<int> nucs;
  };

  class WindowSpec : public NamedSpec {
    public:
      WindowSpec(const std::string & name, const std::vector<SingleSequenceSpec> & domains);

      void set_source(std::vector<SourceSpec *>);
      void allow_similar(const SourceSpec & osource, DBL_TYPE min_sim, DBL_TYPE max_sim);
      void exclude_range(int min_int, int max_int);

      const std::string & get_source() const { return this->source_names[0]; }
      std::vector<int> get_nuc_ids() const { return this->nuc_ids; }
      std::vector<std::vector<int> > get_constraints() const;
      
      void disallow_exclusion() { exclude_allowed = false; }
      bool is_exclusion_allowed() { return exclude_allowed; }
    
    private:
      std::vector<int> nuc_ids;
      std::vector<std::string> domains;
      // need vectors for each of these
      std::vector<std::string> source_names;
      std::vector<int> source_nucs; 
      std::vector<bool> allowed;
      
      bool exclude_allowed = true;
  };

  class DomainListSpec : public NamedSpec {
    public:
      // type_str either "identical", "wobble", or "complementary"
      DomainListSpec(std::string type_str, const std::vector<std::string> & dlist1, 
          const std::vector<std::string> & dlist2, int line=0) :
          NamedSpec("", type_str), dlist1(dlist1), dlist2(dlist2), line(line) {}
          
      std::vector<int> get_nuc_ids1(const SequenceSpec & seqspec) const { 
        return this->get_nuc_ids_helper(this->dlist1, seqspec); 
      }
      std::vector<int> get_nuc_ids2(const SequenceSpec & seqspec) const {
        return this->get_nuc_ids_helper(this->dlist2, seqspec);
      }
      int get_line() const { return this->line; }
    
    private:
      std::vector<int> get_nuc_ids_helper(const std::vector<std::string> & dnames, 
          const SequenceSpec & seqspec) const;

      // Just the lists of domain names
      std::vector<std::string> dlist1; 
      std::vector<std::string> dlist2;
      int line;
  };
  
  /*
   * This is used for general formatting of floats
   */
  std::ostream & scale_format(std::ostream & out);
  std::ostream & exp_format(std::ostream & out);
  std::ostream & flt_format(std::ostream & out);
  std::ostream & longflt_format(std::ostream & out);
  std::ostream & bool_format(std::ostream & out);
  std::ostream & param_format(std::ostream & out);
  
  template <class T, int_if<!is_iterable<T>::value && !is_pair<T>::value> = 0>
  std::string to_string(T t);
  template <class T, int_if<is_pair<T>::value> = 0>
  std::string to_string(T t);
  template <class T, int_if<is_iterable<T>::value> = 0> 
  std::string to_string(T t);
  std::string to_string(bool b);
  std::string to_string(std::string b);
  
  template <class T, int_if<!is_iterable<T>::value && !is_pair<T>::value>>
  inline std::string to_string(T t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
  }
  
  template <class T, int_if<is_pair<T>::value>>
  inline std::string to_string(T t) {
    std::stringstream ss;
    ss << '{' << to_string(t.first) << ": " << to_string(t.second) << "}";
    return ss.str();
  }
  
  template <class T, int_if<is_iterable<T>::value>> 
  inline std::string to_string(T t) {
    std::stringstream ss;
    ss << '(';
    auto i = t.size();
    for (auto & el : t) {
      ss << (--i == 0 ? to_string(el) : to_string(el) + ", ");
    }
    ss << ')';
   return ss.str();
  }
  
  inline std::string to_string(bool b) {
    return b ? "true" : "false";
  };

  inline std::string to_string(std::string str) {
    return str;
  };

};

