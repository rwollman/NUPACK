#pragma once

#include "types.h"

#include <thermo.h>

#include <vector>
#include <string>

namespace nupack {
  class NupackInvariants;
  
  namespace SequenceUtils {
    int bool_to_nuc(std::vector<trinary> in);
    std::vector<int> bool_to_nuc(std::vector<std::vector<trinary> > in);
    std::vector<trinary> nuc_to_bool(int in);
    AllowTable nucs_to_bools(const std::vector<int> & in);
    
    std::string nuc_to_str(const std::vector<int> & in, const int material = RNA);
    BASES char_to_nuc(const char nuc);
    std::vector<int> str_to_nuc(const std::string & in);
    char nuc_to_char(const int nuc, const int material);
    
    BASES get_complement(const int base, const bool allow_wobble);
    void get_complement(const std::vector<int> & in, 
        std::vector<int> & out, const NupackInvariants & invars);
    void get_complement(const std::string & in, 
        std::string & out, const NupackInvariants & invars);
    bool all_are_nucleotides(std::vector<int> sequence);
  };
}
