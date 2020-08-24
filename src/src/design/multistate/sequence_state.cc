#include "sequence_state.h"

#include "design_debug.h"
#include "pathway_utils.h"

#include <array>

namespace nupack {  
void SingleSequenceState::serialize(const SingleSequenceSpec & spec, 
    const NupackInvariants & invars, std::ostream & out) const {
  const std::string & name = spec.get_name();
  std::string seq;
  NUPACK_EXC_CHECK(seq = SequenceUtils::nuc_to_str(this->nucs, invars.material),
      "Error converting nucleotides");

  out << name;
  // TODO use this after validation
  // if (name.back() != '*') out << " ";
  // out << ": " << seq;
  out << " : " << seq; // TODO and remove this line
}

void SingleSequenceState::set_seqs(const std::vector<int> & nucs, 
    const SingleSequenceSpec & nuc_ids) {
  const auto & nuc_id_vec = nuc_ids.get_nuc_ids();
  int n_nucs = nuc_id_vec.size();
  std::vector<int> constraints(n_nucs, BASE_N);

  this->nuc_ids = nuc_id_vec;

  for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    int c_nuc = nuc_id_vec[i_nuc];
    NUPACK_CHECK(c_nuc < nucs.size(), 
        "Invalid nucleotide index specified "
        + to_string(c_nuc) + " >= " + to_string(nucs.size()));

    int c_val = nucs[c_nuc];
    constraints[i_nuc] = c_val;
  }
  this->nucs = constraints;
}

// can avoid RNG when only one nucleotide type is possible. TODO fix after 
// refactoring as results will otherwise be different due to different number
// of calls to RNG.
int SequenceState::random_nucleotide(int con, int cur) {
  NUPACK_CHECK(con >= 0 && con <= 15, "Invalid constraint " + to_string(con));

  std::array<double, 5> weight {{0.0, 0.0, 0.0, 0.0, 0.0}};  
  constexpr std::array<double, 5> def_weight {{0.0, 1.0, 1.0, 1.0, 1.0}};
  const std::vector<std::vector<int>> allowed_seqs {
    {1, 2, 3, 4}, // N
    {1}, // A
    {2}, // C
    {3}, // G
    {4}, // U / T
    {1, 3}, // AG
    {1, 2}, // AC
    {2, 3}, // CG
    {1, 4}, // AU
    {3, 4}, // GU
    {2, 4}, // CU
    {1, 2, 3}, // ACG
    {1, 2, 4}, // ACU
    {1, 3, 4}, // AGU
    {2, 3, 4}, // CGU
  };

  int n_choices = allowed_seqs[con].size();
  for (auto i_con = 0; i_con < n_choices; i_con++) {
    weight[allowed_seqs[con][i_con]] = def_weight[allowed_seqs[con][i_con]];
  }
  weight[cur] = 0;

  auto tot_weight = 0.0;
  for (auto i = 1; i < weight.size(); ++i) tot_weight += weight[i];
  
  auto nuc = 1;
  auto r_num = genrand_real1() * tot_weight;
  auto cur_weight = 0.0;
  for ( ; nuc < 4; nuc++) {
    cur_weight += weight[nuc];
    if (cur_weight > r_num) break;
  }
  return nuc;
}

std::vector<int> SequenceState::random_nucleotides(const std::vector<int> & cons) {
  std::vector<int> rval;
  for (auto & con : cons) rval.push_back(random_nucleotide(con, 0));
  return rval;
}

void SequenceState::fill_in_sequences(const SequenceSpec & spec) {
  domains.clear();
  strands.clear();

  for (auto & dom : spec.get_domains()) domains.emplace_back(nucleotides, dom);
  for (auto & strand : spec.get_strands()) strands.emplace_back(nucleotides, strand);
}

void SequenceState::set_variables(const std::vector<int> & vars, const SequenceSpec & spec) {
  const auto & var_map = spec.get_variables();
  NUPACK_CHECK(var_map.size() <= vars.size(), "Invalid variable map " + 
      to_string(var_map.size()) + " out of " + to_string(vars.size()));
  
  variables = vars;
  nucleotides.resize(var_map.size(), 0);

  for (auto i = 0; i < var_map.size(); i++) {
    nucleotides[i] = vars[var_map[i]] + 1;
  }
  fill_in_sequences(spec);
}

std::vector<int> SequenceState::get_sequence(const std::vector<int> & nuc_ids) const {
  std::vector<int> nucs(nuc_ids.size(), 0);
  for (auto i_nuc = 0; i_nuc < nuc_ids.size(); i_nuc++) {
    NUPACK_CHECK(nuc_ids[i_nuc] >= 0 && nuc_ids[i_nuc] < nucleotides.size(), 
      "Invalid nucleotide id[" + to_string(i_nuc) + "]: " + to_string(nuc_ids[i_nuc]));
    nucs[i_nuc] = nucleotides[nuc_ids[i_nuc]];
  }
  return nucs;
}

void SequenceState::serialize(const SequenceSpec & spec, const NupackInvariants & invars,
    std::ostream & out, int indent, std::string prefix) const {
  std::string namestring(indent, ' ');
  out << namestring << "sequences:" << std::endl;
  
  std::string ind_str(indent + 4, ' ');
  ind_str = prefix + ind_str;

  const auto & domains = spec.get_domains();
  out << ind_str << "domains:" << std::endl;
  for (auto i = 0; i < this->domains.size(); i++) {
    out << ind_str << "    ";
    NUPACK_EXC_CHECK(this->domains[i].serialize(domains[i], invars, out),
        "Error serializing strand");
    out << std::endl;
  }

  const auto & strands = spec.get_strands();
  out << ind_str << "strands:" << std::endl;
  for (auto i = 0; i < this->strands.size(); i++) {
    out << ind_str << "    ";
    NUPACK_EXC_CHECK(this->strands[i].serialize(strands[i], invars, out),
        "Error serializing strand");
    out << std::endl;
  }
}
  
}
