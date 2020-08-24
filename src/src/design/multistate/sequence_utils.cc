#include "sequence_utils.h"

#include "pathway_utils.h"
#include "design_debug.h"
#include "nupack_invariants.h"

namespace nupack { namespace SequenceUtils {
constexpr bool NUCLEOTIDE_CONSTRAINTS_BOOL[][4] = {
  { true, true,  true,  true }, // BASE_N
  { true, false, false, false}, // BASE_A
  {false,  true, false, false}, // BASE_C
  {false, false,  true, false}, // BASE_G
  {false, false, false,  true}, // BASE_U
  { true, false,  true, false}, // BASE_R
  { true,  true, false, false}, // BASE_M
  {false,  true,  true, false}, // BASE_S
  { true, false, false,  true}, // BASE_W
  {false, false,  true,  true}, // BASE_K
  {false,  true, false,  true}, // BASE_Y
  { true,  true,  true, false}, // BASE_V
  { true,  true, false,  true}, // BASE_H
  { true, false,  true,  true}, // BASE_D
  {false,  true,  true,  true}, // BASE_B
};

std::vector<trinary> nuc_to_bool(int in) {
  NUPACK_CHECK(in >= 0 && in < 15, 
      "Invalid nucleotide constraint: " + to_string(in));
  return std::vector<trinary>(NUCLEOTIDE_CONSTRAINTS_BOOL[in], 
      NUCLEOTIDE_CONSTRAINTS_BOOL[in] + 4);
}

AllowTable nucs_to_bools(const std::vector<int> & in) {
  std::vector<std::vector<trinary>> bools;
  for (auto & i : in) {
    bools.push_back(nuc_to_bool(i));
  }
  return bools;
}

int bool_to_nuc(std::vector<trinary> in) {
  constexpr int m = 15, n = 4;
  NUPACK_CHECK(in.size() == n, 
      "Invalid nucleotide booleans provided, must be size " + to_string((int)n));
  bool match = true;
  
  int i = 0;
  for ( ; i < m; i++) {
    match = true;
    for (auto j = 0; j < n; j++) {
      if ((NUCLEOTIDE_CONSTRAINTS_BOOL[i][j] && !in[j])
          || (!NUCLEOTIDE_CONSTRAINTS_BOOL[i][j] && in[j])) {
        match = false;
      }
    }

    if (match) break;
  }
  if (!match) i = BASE_NONE; 
  return i;
}

std::vector<int> bool_to_nuc(std::vector<std::vector<trinary> > in) {
  const unsigned int num_bases = 4;
  std::vector<int> out;
  for (auto elem : in) {
    out.push_back((elem.size() == num_bases) ? bool_to_nuc(elem): BASE_NONE);
  }
  return out;
}

std::string nuc_to_str(const std::vector<int> & in, const int material) {
  int n_nucs = in.size();
  std::string out(n_nucs, 'N');

  for (int i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    out[i_nuc] = nuc_to_char(in[i_nuc], material);
  }
  return out;
}

char nuc_to_char(const int nuc, const int material) {
    char U_char = 'U';
    if (material == DNA) {
      U_char = 'T';
    }

    switch (nuc) {
      case BASE_N: return 'N';
      case BASE_A: return 'A';
      case BASE_C: return 'C';
      case BASE_G: return 'G';
      case BASE_U: return U_char;
      case BASE_R: return 'R';
      case BASE_Y: return 'Y';
      case BASE_M: return 'M';
      case BASE_K: return 'K';
      case BASE_S: return 'S';
      case BASE_W: return 'W';
      case BASE_V: return 'V';
      case BASE_B: return 'B';
      case BASE_H: return 'H';
      case BASE_D: return 'D';
      case BASE_NONE: return '-'; break;
      case STRAND_PLUS: return '+'; break;
      default:
        NUPACK_ERROR("Invalid nucleotide encountered");
  }
}

BASES char_to_nuc(const char nuc) {
  switch (toupper(nuc)) {
    case 'N': return BASE_N;
    case 'A': return BASE_A;
    case 'C': return BASE_C;
    case 'G': return BASE_G;
    case 'T': return BASE_U;
    case 'U': return BASE_U;
    case 'R': return BASE_R;
    case 'Y': return BASE_Y;
    case 'M': return BASE_M;
    case 'K': return BASE_K;
    case 'S': return BASE_S;
    case 'W': return BASE_W;
    case 'V': return BASE_V;
    case 'B': return BASE_B;
    case 'H': return BASE_H;
    case 'D': return BASE_D;
    default:
      NUPACK_ERROR("Invalid nucleotide encountered " + to_string((int) nuc) + ".");
  }
}

bool all_are_nucleotides(std::vector<int> sequence) {
  for (auto s : sequence) {
    if (s == BASE_A || s == BASE_C || s == BASE_G || s == BASE_U)
      continue;
    return false;
  }
  return true;
}

std::vector<int> str_to_nuc(const std::string & in) {
  int n_nucs = in.size();
  std::vector<int> out(n_nucs, 0);
  for (int i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    try { out[i_nuc] = char_to_nuc(in[i_nuc]); }
    catch (NupackException & e) {
      e.print_message();
      NUPACK_ERROR("Encountered in string " + in + ".");
    }
  }
  return out;
}

BASES get_complement(const int base, const bool allow_wobble) {
  switch (base) {
    case BASE_N: return BASE_N;
    case BASE_A: return BASE_U;
    case BASE_C: return BASE_G;
    case BASE_G: 
      if (allow_wobble) {
        return BASE_Y;
      } else {
        return BASE_C;
      }
    case BASE_U: 
      if (allow_wobble) {
        return BASE_R;
      } else {
        return BASE_A;
      }
    case BASE_R: return BASE_Y;
    case BASE_M: return BASE_K;
    case BASE_S: return BASE_S;
    case BASE_W: return BASE_W;
    case BASE_K: return BASE_N;
    case BASE_Y: return BASE_R;
    case BASE_V: return BASE_B;
    case BASE_H: return BASE_D;
    case BASE_D: return BASE_N;
    case BASE_B: return BASE_N;
    default:
      NUPACK_ERROR("Invalid nucleotide encountered");
      break;
  }
}

void get_complement(const std::vector<int> & in, 
    std::vector<int> & out, const NupackInvariants & invars) {
  int n_nucs = in.size();
  out.resize(n_nucs, 0);
  for (int i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    out[n_nucs - i_nuc - 1] = get_complement(in[i_nuc], invars.allow_wobble);
  }
}

void get_complement(const std::string & in, 
    std::string & out, const NupackInvariants & invars) {
  std::vector<int> tempin;
  std::vector<int> tempout;

  try {
    tempin = str_to_nuc(in);
    get_complement(tempin, tempout, invars);
    out = nuc_to_str(tempout, invars.material);
  } catch (NupackException & e) {
    e.print_message(std::cerr);
    NUPACK_ERROR("Error finding complement");
  }
}
}}
