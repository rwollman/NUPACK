#pragma once

#include "physical_spec.h"

#include <vector>
#include <string>
#include <memory>


namespace nupack {
  class SequenceSpec;
  class SingleSequenceSpec;
      /**
       * This class holds a fully defined sequence for a domain or strand
       */
  class SingleSequenceState {
    public:
      
      SingleSequenceState() {}
      SingleSequenceState(const std::vector<int> & nucs, 
          const SingleSequenceSpec & nuc_ids) {
        set_seqs(nucs, nuc_ids);
      }

      void serialize(const SingleSequenceSpec & spec, 
          const NupackInvariants & invars, std::ostream & out) const;
      void set_seqs(const std::vector<int> & nucs, 
          const SingleSequenceSpec & nuc_ids);

      const std::vector<int> & get_nuc_ids() const { return this->nuc_ids; }
      const std::vector<int> & get_nucs() const { return this->nucs; }

    protected:
      std::vector<int> nucs;
      std::vector<int> nuc_ids;
  };

      /**
       * This class holds a fully defined sequence for a set of domains 
       * and strands.
       */
  class SequenceState {
    public:
      SequenceState() {};

      // int init_random(const SequenceSpec & spec);
      const std::vector<int> & get_variables() const { return this->variables; }
      void set_variables(const std::vector<int> & vars, const SequenceSpec & spec);

      void serialize(const SequenceSpec & spec, 
          const NupackInvariants & invars, std::ostream & out,
          int indent = 0, std::string prefix = "") const;

      const std::vector<SingleSequenceState> & get_strands() const { return this->strands; }

      std::vector<int> get_sequence(const std::vector<int> & nuc_ids) const;

      int n_nucs() const;

      static int random_nucleotide(int con, int cur);
      static std::vector<int> random_nucleotides(const std::vector<int> & nucleotides);

    protected:
      void fill_in_sequences(const SequenceSpec & spec);

      std::vector<int> variables;
      std::vector<int> nucleotides;
      std::vector<SingleSequenceState> domains;
      std::vector<SingleSequenceState> strands;
  };    

}

