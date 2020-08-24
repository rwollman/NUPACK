#include "nupack_invariants.h"
#include "pathway_utils.h"
#include "physical_spec.h"

#include <sys/time.h>

namespace nupack {
NupackInvariants::NupackInvariants() {
  time_t curtime;
  time(&curtime);
  std::string timestring(ctime(&curtime));
  timeval starttime;
  gettimeofday(&starttime, NULL);

  timestring.resize(timestring.length() - 1);
  this->start_timestamp = timestring;

  this->start_time = starttime.tv_sec + 1e-6 * starttime.tv_usec;
}

void NupackInvariants::deduce_h_split() {
  if (H_split == -1) {
    H_split = (material == RNA || material == RNA37) 
        ? NUPACK_DEF_H_SPLIT_RNA
        : NUPACK_DEF_H_SPLIT_OTHER;
  }
  
}


bool NupackInvariants::opt_time_elapsed() const {
  timeval curtimestruct;
  gettimeofday(&curtimestruct, NULL);

  DBL_TYPE curtime = curtimestruct.tv_sec + 1e-6 * curtimestruct.tv_usec;

  return (curtime - this->start_time > this->allowed_opt_time);
}

std::string NupackInvariants::mat_str() const {
  if (this->material == DNA) {
    return std::string("dna");
  } else if (this->material == RNA) {
    return std::string("rna1995");
  }  else if (this->material == RNA37) {
    return std::string("rna1999");
  }

  return this->material_string;
}

std::string NupackInvariants::dangle_str() const {
  if (this->dangle_type == 0) {
    return std::string("none");
  } else if (this->dangle_type == 1) {
    return std::string("some");
  } 

  return std::string("full");
}

#ifdef JSONCPP_FOUND
Json::Value NupackInvariants::make_json_value() const { 
  Json::Value val;
  Json::Value tmpval;

  tmpval["start time"] = this->start_timestamp;
  tmpval["elapsed[sec]"] = (double) this->elapsed_time;
  tmpval["input prefix"] = this->file_prefix;

  val["general"] = tmpval;

  tmpval.clear();
  tmpval["temperature[K]"] = (double) this->temperature;
  tmpval["material"] = this->mat_str();
  tmpval["sodium[M]"] = (double) this->sodium;
  tmpval["magnesium[M]"] = (double) this->magnesium;
  tmpval["dangles"] = this->dangle_str();
  tmpval["min_ppair"] = (double) this->min_ppair;
  tmpval["use long helix"] = this->use_long_helix;

  val["physical"] = tmpval;

  tmpval.clear();
  tmpval["seed"] = this->seed;
  tmpval["population"] = this->N_population;
  tmpval["N_trials"] = this->N_trials;
  tmpval["H_split"] = this->H_split;
  tmpval["N_split"] = this->N_split;
  tmpval["print_steps"] = this->print_steps;
  tmpval["print_leaves"] = this->print_leaves;
  tmpval["M_bad"] = (double) this->M_bad;
  tmpval["M_reopt"] = (double) this->M_reopt;
  tmpval["f_split"] = (double) this->f_split;
  tmpval["f_passive"] = (double) this->f_passive;
  tmpval["f_stringent"] = (double) this->f_stringent;
  tmpval["f_redecomp"] = (double) this->f_redecomp;
  tmpval["f_refocus"] = (double) this->f_refocus;
  tmpval["gc_init_prob"] = (double) this->gc_init_prob;
  tmpval["cutoff[sec]"] = (double) this->allowed_opt_time;

  tmpval["allow mismatch"] = this->allow_mismatch;
  tmpval["allow wobble"] = this->allow_wobble;
  tmpval["disable weight"] = this->disable_defect_weights;
  tmpval["disable focus"] = this->disable_focus;
  tmpval["forbid split"] = this->forbid_splits;
  tmpval["redecompose"] = this->redecompose;
  tmpval["add default stops"] = this->add_default_stops;

  val["design"] = tmpval;

  return val;
}
#endif // JSONCPP_FOUND

void NupackInvariants::serialize(std::ostream &out, int indent, std::string prefix) const{

  std::string pref_str(indent, ' ');

  pref_str += prefix;

  out << pref_str << "general : "                                               << std::endl;
  out << pref_str << "    version  : "      << NUPACK_VERSION                   << std::endl; 
  out << pref_str << "    start time    : " << start_timestamp                  << std::endl;
  out << pref_str << "    elapsed[sec]  : " << flt_format   << elapsed_time     << std::endl;
  out << pref_str << "    input prefix  : " << param_format << file_prefix      << std::endl; 

  out << pref_str << "physical : " << std::endl;
  out << pref_str << "    temperature[K]: " << flt_format   << temperature      << std::endl;
  out << pref_str << "    material      : " << param_format << mat_str()        << std::endl;
  out << pref_str << "    sodium[M]     : " << flt_format   << sodium           << std::endl;
  out << pref_str << "    magnesium[M]  : " << flt_format   << magnesium        << std::endl;
  out << pref_str << "    dangles       : " << param_format << dangle_str()     << std::endl;
  out << pref_str << "    min_ppair     : " << flt_format   << min_ppair        << std::endl;
  // out << pref_str << "    use long helix: " << bool_format  << use_long_helix   << std::endl;

  out << pref_str << "design : " << std::endl;
  out << pref_str << "    seed          : " << param_format << seed             << std::endl;
  out << pref_str << "    population    : " << param_format << N_population     << std::endl;
  // out << pref_str << "    N_trials      : " << param_format << N_trials         << std::endl;
  out << pref_str << "    H_split       : " << param_format << H_split          << std::endl;
  out << pref_str << "    N_split       : " << param_format << N_split          << std::endl; 
  // out << pref_str << "    print steps   : " << param_format << print_steps      << std::endl;
  //out << pref_str << "    print leaves  : " << bool_format << print_leaves << std::endl;
  out << pref_str << "    M_bad         : " << param_format << M_bad            << std::endl;
  out << pref_str << "    M_reopt       : " << param_format << M_reopt          << std::endl;
  out << pref_str << "    M_reseed      : " << param_format << M_reseed         << std::endl;
  out << pref_str << "    f_split       : " << flt_format   << f_split          << std::endl;
  out << pref_str << "    f_passive     : " << flt_format   << f_passive        << std::endl;
  out << pref_str << "    f_stringent   : " << flt_format   << f_stringent      << std::endl;
  out << pref_str << "    f_redecomp    : " << flt_format   << f_redecomp       << std::endl;
  out << pref_str << "    f_refocus     : " << flt_format   << f_refocus        << std::endl;
  out << pref_str << "    dg_clamp      : " << flt_format   << dG_clamp         << std::endl;
  // out << pref_str << "    gc_init_prob  : " << flt_format   << gc_init_prob     << std::endl;
  out << pref_str << "    cutoff[sec]   : " << flt_format   << allowed_opt_time << std::endl;

  // out << pref_str << "    allow mismatch: " << bool_format << allow_mismatch  << std::endl;
  out << pref_str << "    allow wobble  : " << bool_format  << allow_wobble     << std::endl; 
  //out << pref_str << "    disable weight: " << bool_format << disable_defect_weights << std::endl;
  //out << pref_str << "    disable focus : " << bool_format << disable_focus   << std::endl; 
  //out << pref_str << "    forbid split  : " << bool_format << forbid_splits   << std::endl;
  //out << pref_str << "    redecompose   : " << bool_format << redecompose     << std::endl;
  out << pref_str << "    add default stops   : " << bool_format << add_default_stops   << std::endl;

}
}
