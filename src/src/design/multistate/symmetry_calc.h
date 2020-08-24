#pragma once

#include "named_spec.h"

#include <thermo.h>

#include <string>
#include <vector>
#include <map>

namespace nupack {
  class SequenceSpec;
  class SequenceState;
  
  class SingleSymmetrySpec : public NamedSpec {
    public:
      SingleSymmetrySpec(const std::string & name, const std::vector<std::string> & domains);

      void resolve(const SequenceSpec & seqs);

      void set_word_len(std::vector<int> lengths);

      void set_weights(DBL_TYPE weight);
      void set_weights(std::vector<DBL_TYPE> weights);

      const std::vector<std::vector<int> > & get_nuc_ids() const { return this->nuc_ids; }
      const std::vector<int> & get_word_len() const { return this->word_len; }
      const std::vector<DBL_TYPE> & get_weights() const { return this->weights; }

    private:
      std::vector<std::string> names;
      std::vector<std::vector<int> > nuc_ids;
      std::vector<int> word_len;
      std::vector<DBL_TYPE> weights;
      DBL_TYPE weightscale;

      const static DBL_TYPE default_weight;
  };

  class SingleSymmetryResult {
    public:
      void evaluate(const SequenceState & seqs, 
          const SingleSymmetrySpec & spec);

      DBL_TYPE get_defect() const { return this->defect; }
      std::map<int, DBL_TYPE> get_nuc_defects() const { return this->nuc_defects; }

    private:
      std::vector<std::vector<DBL_TYPE> > symmetries;
      std::vector<std::vector<int> > last_seqs;
      std::map<int, DBL_TYPE> nuc_defects;
      DBL_TYPE defect;
  };

  class SymmetrySpec {
    public:
      void add_objective(const SingleSymmetrySpec & spec) { this->specs.push_back(spec); }
      const std::vector<SingleSymmetrySpec> & get_specs() const { return this->specs; }
    private:
      std::vector<SingleSymmetrySpec> specs;
  };

  class SymmetryResult {
    public:
      void evaluate(const SequenceState & seqs, const SymmetrySpec & spec);
      const std::vector<SingleSymmetryResult> & get_results() const { return this->results; }
    private:
      std::vector<SingleSymmetryResult> results;
  };

}

