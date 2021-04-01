#pragma once

#include "types.h"

#include "adapter.h"

#include <vector>
#include <string>

namespace nupack {
namespace custom_csp {

using namespace design;

// Nucleotides
enum BASES {
    BASE_NONE = -1, // No base is possible
    BASE_A = 0,
    BASE_C = 1,
    BASE_G = 2,
    BASE_T = 3,
    BASE_U = 3,
    BASE_R = 4, // AG
    BASE_M = 5, // AC
    BASE_S = 6, // CG
    BASE_W = 7, // AU
    BASE_K = 8, // GU
    BASE_Y = 9, // CU
    BASE_V = 10, // ACG
    BASE_H = 11, // ACU
    BASE_D = 12, // AGU
    BASE_B = 13, // CGU
    BASE_N = 14, //ACTG
    STRAND_PLUS = 15, // Strand break
};

struct NupackInvariants;

namespace SequenceUtils {

int bool_to_nuc(vec<trinary> const & in);
vec<int> bool_to_nuc(vec<vec<trinary> > const & in);
vec<trinary> nuc_to_bool(int in);
AllowTable nucs_to_bools(const vec<int> & in);

string nuc_to_str(const vec<int> & in, const int material = RNA);
BASES char_to_nuc(const char nuc);
vec<int> str_to_nuc(const string & in);

char nuc_to_char(const int nuc, const int material);

BASES get_complement(const int base, const bool allow_wobble);
void get_complement(const vec<int> & in,
                    vec<int> & out, const NupackInvariants & invars);
void get_complement(const string & in,
                    string & out, const NupackInvariants & invars);
bool all_are_nucleotides(vec<int> sequence);

};

}
}
