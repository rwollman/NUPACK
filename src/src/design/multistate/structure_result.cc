#include "structure_result.h"

#include "complex_result.h"
#include "node_result.h"
#include "sequence_state.h"

#include "structure_utils.h"
#include "sequence_utils.h"
#include "pair_probabilities.h"

#include "design_debug.h"
#include "pathway_utils.h"

#include <fstream>

namespace nupack {

std::string StructureResult::get_structure_string(int i_struc) const {
  return StructureUtils::pairs_to_dpp(this->structures[i_struc], this->breaks);
}

std::string StructureResult::get_sequence_string(
    const NupackInvariants & invars) const {
  std::vector<int> breaks = this->breaks;
  breaks.push_back(this->sequence.size());

  std::vector<int> seqnucs;
  int i_nuc = 0;
  for (auto br : breaks) {
    seqnucs.insert(seqnucs.end(), this->sequence.begin() + i_nuc, 
        this->sequence.begin() + br);
    seqnucs.push_back(STRAND_PLUS);
    i_nuc = br;
  }
  seqnucs.pop_back();

  return SequenceUtils::nuc_to_str(seqnucs, invars.material);
}

void StructureResult::serialize(const StructureSpec & spec, 
    const PhysicalParams & params, const NupackInvariants & invars,
    std::ostream & out, 
    int indent, std::string prefix) const {
  std::string ind_str(indent, ' ');
  std::string name_ind(indent, ' ');
  ind_str = prefix + ind_str;

  if (indent > 2) name_ind = prefix + std::string(indent - 2, ' ') + "- "; 

  std::string seq = this->get_sequence_string(invars);

  for (auto i_tar = 0; i_tar < this->defects.size(); i_tar++) {
    out << name_ind << "name: " << spec.get_names()[i_tar] << std::endl;
    std::string struc = this->get_structure_string(i_tar);
    out << ind_str << "structure: " << struc << std::endl;
    out << ind_str << "sequence:  " << seq << std::endl;
    if (this->tree.get_max_depth() > 0) {
      int ind = 0;
      this->tree.serialize(out, indent + 11, prefix, ind);
    }
    out << ind_str << "defect[nt]: " << flt_format << this->defects[i_tar] << std::endl;
    out << ind_str << "normalized defect: " 
      << flt_format << this->get_normalized_defect(i_tar) << std::endl;
    // if (this->defects.size() == 1) {
      out << ind_str << "target free energy[kcal/(mol K)]: " << flt_format 
          << this->struc_energies[i_tar] << std::endl;
      DBL_TYPE complex_ene = -kB * params.temperature * LOG_FUNC(this->pfunc);
      out << ind_str << "target probability: " << flt_format 
          << EXP_FUNC((complex_ene - this->struc_energies[i_tar]) 
          / (kB * params.temperature)) << std::endl;
    // }
    out << ind_str << "complex free energy[kcal/(mol K)]: " << flt_format 
        << -kB * params.temperature * LOG_FUNC(this->pfunc) << std::endl;
  }

  if (this->defects.size() == 0) {
    out << name_ind << "name: " << spec.get_name() << std::endl;
    out << ind_str << "complex free energy[kcal/(mol K)]: " << flt_format << 
        -kB * params.temperature * LOG_FUNC(this->pfunc) << std::endl;
  }
}

void StructureResult::print_strucfiles(const StructureSpec & spec, 
    const PhysicalParams & params, const NupackInvariants & invars, 
    std::string file_prefix) const {
  const std::vector<std::string> & names = spec.get_names();
  auto it = names.begin();
  for (auto i = 0; it != names.end(); ++i, ++it) {
    std::string ppair_filename = file_prefix + "_" + *it + ".ppairs";
    std::string json_filename = file_prefix + "_" + *it + ".json";

    if (invars.print_ppairs) {
      std::fstream ppair_file(ppair_filename, std::fstream::out);

      ppair_file << COMMENT_STRING << " " << "Version: " << NUPACK_VERSION << std::endl;
      ppair_file << COMMENT_STRING << " " << "Program: multistatedesign" << std::endl;
      ppair_file << COMMENT_STRING << " " << "Start time: " << invars.start_timestamp << std::endl;
      ppair_file << COMMENT_STRING << " " << "Parameters: " << invars.mat_str() << std::endl;
      ppair_file << COMMENT_STRING << " " << "Dangles: " << invars.dangle_str() << std::endl;
      ppair_file << COMMENT_STRING << " " << "Temperature (C): " << flt_format << params.temperature << std::endl;
      ppair_file << COMMENT_STRING << " " << "Sodium concentration (M): " << invars.sodium << std::endl;
      ppair_file << COMMENT_STRING << " " << "Magnesium concentration (M): " << invars.magnesium << std::endl;
      ppair_file << COMMENT_STRING << " " << "Name: " << *it << std::endl;
      ppair_file << COMMENT_STRING << " " << "Sequence:  " << this->get_sequence_string(invars) << std::endl;
      ppair_file << COMMENT_STRING << " " << "Structure: " << this->get_structure_string(i) << std::endl;
      ppair_file << COMMENT_STRING << " " << "v(pi): " << this->symmetry << std::endl;
      ppair_file << COMMENT_STRING << " " << "Structure free energy: " << this->struc_energies[i] << std::endl;
      ppair_file << COMMENT_STRING << " " << "Ensemble free energy: " << (- kB * params.temperature * LOG_FUNC(this->pfunc)) << std::endl;
      ppair_file << COMMENT_STRING << " " << "Partition function: " << exp_format << this->pfunc << std::endl;
      ppair_file << COMMENT_STRING << " " << "Ensemble defect: " << flt_format << this->defects[i] << std::endl;
      ppair_file << COMMENT_STRING << " " << "Normalized ensemble defect: " << flt_format << (this->defects[i] / this->size()) << std::endl;
      ppair_file << std::endl;
      ppair_file << this->size() << std::endl;
      
      this->ppairs.serialize(ppair_file, this->size());
      ppair_file.close();
    }

    if (invars.print_json) {
      std::fstream json_file(json_filename, std::fstream::out);
      
      int indent_size = 4;
      int indent_level = 0;
      std::string material = "rna";
      if (invars.mat_str() == "dna" || invars.mat_str() == "dna1998") material = "dna"; 
      
      json_file << "{" << std::endl;
      indent_level++;
      
      std::string cur_indent(indent_level * indent_size, ' ');
 
      std::vector<std::string> snames = spec.get_strand_names();
      std::stringstream sname_ss;
      sname_ss << "[ \"" << snames[0] << "\"";
      for (auto sid = snames.begin() + 1; sid != snames.end(); ++sid) {
        sname_ss << ", \"" << *sid << "\"";
      }
      sname_ss << "]";
      json_file << cur_indent << "\"nupack version\": \"" << NUPACK_VERSION << "\"" << "," << std::endl;
      json_file << cur_indent << "\"program\": \"cpathway\"" << "," << std::endl;
      json_file << cur_indent << "\"start time\": \"" << invars.start_timestamp << "\"," << std::endl;
      json_file << cur_indent << "\"parameters\": \"" << invars.mat_str() << "\"" << "," << std::endl;
      json_file << cur_indent << "\"material\": \"" << material << "\"" << "," << std::endl;
      json_file << cur_indent << "\"dangles\": \"" << invars.dangle_str() << "\"" << "," << std::endl;
      json_file << cur_indent << "\"temperature (C)\": " << flt_format 
          << params.temperature - ZERO_C_IN_KELVIN << "," << std::endl;
      json_file << cur_indent << "\"name\": \"" << *it << "\"," << std::endl;
      json_file << cur_indent << "\"strand names\": " << sname_ss.str() << "," << std::endl;
      json_file << cur_indent << "\"sodium concentration (M)\": " << invars.sodium << "," << std::endl;
      json_file << cur_indent << "\"magnesium concentration (M)\": " << invars.magnesium << "," << std::endl;
      json_file << cur_indent << "\"sequence\": \"" << this->get_sequence_string(invars) << "\"" << "," << std::endl;
      json_file << cur_indent << "\"structure\": \"" << this->get_structure_string(i) << "\"" << "," << std::endl;
      json_file << cur_indent << "\"v(pi)\":" << this->symmetry << "," << std::endl;
      json_file << cur_indent << "\"energypreamble\": \"Free energy of secondary structure:\"" << "," << std::endl;
      json_file << cur_indent << "\"target free energy\": " << 
        this->struc_energies[i] << "," << std::endl;
      json_file << cur_indent << "\"ensemble free energy\": " <<
        (- kB * params.temperature * LOG_FUNC(this->pfunc) ) << "," << std::endl;
      json_file << cur_indent << "\"partition function\": " << exp_format << this->pfunc << "," << std::endl;
      json_file << cur_indent << "\"ensemble defect\": " << flt_format << this->defects[i] << "," << std::endl;
      json_file << cur_indent << "\"normalized ensemble defect\": " << flt_format << 
        (this->defects[i] / this->size()) << "," << std::endl;
      json_file << cur_indent << "\"nucleotide probabilities\":" << std::endl;
      json_file << cur_indent << "[" << std::endl;
      int j;
      indent_level++;
      cur_indent.assign(indent_level * indent_size, ' ');
      for (j = 0; j < this->size() - 1; j++) {
        json_file << cur_indent << (1 - this->pos_defects[i][j]) << "," << std::endl;
      }
      json_file << cur_indent << (1 - this->pos_defects[i][j]) << std::endl;
      indent_level--;

      cur_indent.assign(indent_level * indent_size, ' ');
      json_file << cur_indent << "]" << std::endl;
      json_file << "}" << std::endl;
      json_file.close();
    }
  }
}

void StructureResult::evaluate(const SequenceState & seqs, 
    const StructureSpec & spec,
    const PhysicalParams & params, const NupackInvariants & invars) {
  
  const auto & strands = seqs.get_strands();
  const auto & strand_ids = spec.get_strands();

  std::vector<int> fullseq;
  this->breaks.clear();
  for (auto c_str : strand_ids) {
    NUPACK_CHECK(c_str >= 0 && c_str < strands.size(), "Invalid strand id " + to_string(c_str) + ", negative or greater than " + to_string(strands.size()));
    
    const auto & curseq = strands[c_str].get_nucs();
    const auto & cur_nuc_ids = strands[c_str].get_nuc_ids();
    append(fullseq, curseq);
    append(nuc_ids, cur_nuc_ids);
    this->breaks.push_back(fullseq.size());
  }
  this->breaks.pop_back();

  int n_nucs = spec.size();
  NUPACK_CHECK(fullseq.size() == n_nucs, "structure and sequence lengths don't agree. " + to_string(fullseq.size()) + " != " + to_string(n_nucs));

  // Compare old and new sequences. Mark any changed nodes
  // If the sequence length changed, mark all nodes as changed
  this->sequence = fullseq;
  this->f_sequence = fullseq;
  this->structures = spec.get_structures();
  this->target = spec.get_target();
  this->strands = spec.get_strands();
  this->domain_map = spec.get_domain_map();
  this->strand_map = spec.get_strand_map();
  this->struc_ids = spec.get_struc_ids();
  this->symmetry = spec.get_symmetry();

  this->pfunc = 0;
  this->ppairs.clear();
  this->eval_time = 0;
  this->nuc_defects.clear();

  // if (nullptr == &this->tree) this->tree = NodeResult();
  this->tree.evaluate(spec.get_tree(), seqs, *this, spec, params, invars);

  this->pfunc = this->get_pfunc(params);
  this->ppairs = this->tree.collect_pair_probs(params, invars);
  this->eval_time = this->tree.collect_eval_times();
  this->update_defects();

  this->struc_energies.resize(this->structures.size());
  for (auto i_tar = 0; i_tar < this->structures.size(); i_tar++) {
    this->struc_energies[i_tar] = this->get_energy(this->structures[i_tar], this->breaks, this->sequence, this->symmetry, params, invars);
  }
}

const Map & StructureResult::get_nuc_defects(int i_target) const {
  NUPACK_DEBUG_CHECK(i_target < this->nuc_defects.size(),
      "Invalid target specified in StructureResult::get_nuc_defects");
  return this->nuc_defects[i_target];
}

DBL_TYPE StructureResult::get_normalized_defect(int i_target) const {
  return this->get_structural_defect(i_target) / this->size();
}

DBL_TYPE StructureResult::get_structural_defect(int i_target) const {
  NUPACK_CHECK(i_target < this->defects.size(), 
      "Invalid target " + to_string(i_target) + " out of " 
      + to_string(this->defects.size()));
  return this->defects[i_target];
}

DBL_TYPE StructureResult::get_pfunc(const PhysicalParams & params) const {
  DBL_TYPE full_pfunc = this->tree.get_pfunc();
  DBL_TYPE bimolec_ene = (BIMOLECULAR + SALT_CORRECTION) * (this->breaks.size());
  
  full_pfunc *= EXP_FUNC(-(bimolec_ene) / (kB * params.temperature));
  full_pfunc /= this->symmetry;

  return full_pfunc;
}

DBL_TYPE StructureResult::get_energy(const std::vector<int> & struc, 
    const std::vector<int> & breaks, const std::vector<int> & seq, int sym, 
    const PhysicalParams & params, const NupackInvariants & invars) {
  std::vector<int> struc_vec(struc.begin(), struc.end());

  std::vector<int> seq_vec(seq.size() + 1 + breaks.size(), 0);
  int j = 0;
  for (auto i = 0, i_br = 0; i < seq.size(); i++) {
    if (i_br < breaks.size() && i == breaks[i_br]) {
      seq_vec[j] = STRAND_PLUS;
      j++;
      i_br++;
    }
    seq_vec[j] = seq[i];
    j++;
  }
  seq_vec[j] = -1;
  
  return naEnergyPairsOrParensFullWithSym(struc_vec.data(), NULL, seq_vec.data(), 
      invars.material, invars.dangle_type, params.temperature - ZERO_C_IN_KELVIN, 
      sym, invars.sodium, invars.magnesium, invars.use_long_helix);
}

void StructureResult::update_defects() {
  this->pos_defects.resize(this->structures.size());
  this->nuc_defects.resize(this->structures.size());
  this->defects.resize(this->structures.size());

  for (auto i_tar = 0; i_tar < this->structures.size(); i_tar++) {
    auto struc_id = this->struc_ids[i_tar];

    this->defects[i_tar] = 0;
    this->pos_defects[i_tar] = this->ppairs.get_nuc_defects(target);

    int n_nucs = this->structures[0].size();
    for (auto i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      std::vector<int> ind{struc_id, this->strand_map[i_nuc], this->domain_map[i_nuc], i_nuc, this->nuc_ids[i_nuc]};

      // this->nuc_defects[i_tar][ind] = 0;
      this->defects[i_tar] += this->pos_defects[i_tar][i_nuc];
      this->nuc_defects[i_tar][ind] = this->pos_defects[i_tar][i_nuc];
    }
  }
}

void StructureResult::replace_node(const StructureResult & other, int k,
    const PhysicalParams & params, const NupackInvariants & invars) {
  this->tree.replace_node(other.tree, k, params, invars);
  this->pfunc = this->get_pfunc(params);
  this->ppairs = this->tree.collect_pair_probs(params, invars);
  this->update_defects();
}
}
