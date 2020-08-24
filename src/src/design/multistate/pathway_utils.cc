
#include "sys/time.h"

#include "design_spec.h"
#include "design_debug.h"
#include "pathway_utils.h"
#include "physical_spec.h"
#include "algorithms.h"

#include <iomanip>

namespace nupack {
DBL_TYPE get_current_time() {
  DBL_TYPE time = 0;
  timeval timev;
  gettimeofday(&timev, NULL);
  time = timev.tv_sec + 1e-6 * timev.tv_usec;
  return time;
}

int sample_weighted_int(const std::vector<DBL_TYPE> & weights) {
  DBL_TYPE tot_weight = 0;
  for (auto weight : weights) tot_weight += weight;
  DBL_TYPE stop = genrand_real1() * tot_weight;
  tot_weight = 0;
  int i = 0;
  for (i = 0; i < weights.size() - 1; i++) {
    tot_weight += weights[i];
    if (tot_weight > stop) break;
  }
  return i;
}

WindowSpec::WindowSpec(const std::string & name, const std::vector<SingleSequenceSpec> & domains) :
    NamedSpec(name, "window") {
  for (auto & dom : domains) {
    const auto & nuc_ids = dom.get_nuc_ids();
    append(this->nuc_ids, nuc_ids);
    this->domains.push_back(dom.get_name());
  }
}

void WindowSpec::set_source(std::vector<SourceSpec *> sources) {
  // need to disallow exclusions in the case of multiple sources
  if (sources.size() > 1) disallow_exclusion();
  
  // vector of false on last window_size - 1 potential substring start points
  // for EACH source generalizes having a vector of allowed of length source
  // size - window size + 1 in single source case
  auto window_size = nuc_ids.size();
  std::vector<bool> disallowed(window_size - 1, false);
  allowed.clear();
  
  for (auto s : sources) {
    const auto & source = *s;
    source_names.push_back(source.get_name());
    const std::vector<int> & nucs = source.get_nucs();
    int n_nucs = nucs.size();
    NUPACK_CHECK(n_nucs >= window_size, 
        "Number of nucleotides in window \"" + get_name() + 
        "\" cannot exceed legnth of source \"" + source_names.back() + "\"");
    std::vector<bool> current_allowed(n_nucs - window_size + 1, true);
    append(allowed, current_allowed);
    append(allowed, disallowed);
    append(source_nucs, nucs);
  }
}

void WindowSpec::allow_similar(const SourceSpec & osource, 
    DBL_TYPE min_sim, DBL_TYPE max_sim) {
  int w_len = this->nuc_ids.size();
  const auto & tnucs = source_nucs;
  const auto & onucs = osource.get_nucs();
  NUPACK_CHECK(tnucs.size() - w_len + 1 == this->allowed.size(),
      "Invalid tsource when setting allow_similar");

  NUPACK_DEBUG("Min sim: " << min_sim << "  Max sim: " << max_sim);
  for (auto i = 0; i < tnucs.size() - w_len + 1; i++) {
    DBL_TYPE sim = 0;
    for (auto j = 0; j < onucs.size() - w_len + 1; j++) {
      int n_same = 0;
      for (auto k = 0; k < w_len; k++) {
        if (tnucs[i + k] == onucs[j + k]) n_same++;
      }
      DBL_TYPE cur_frac = ((DBL_TYPE) n_same) / w_len;
      sim = std::max(cur_frac, sim);
    }
    if (sim < min_sim || sim > max_sim) this->allowed[i] = false;
  }
}

void WindowSpec::exclude_range(int min_exclude, int max_exclude) {
  NUPACK_CHECK(is_exclusion_allowed(), 
    "Exclusion is not allowed because multiple sources are used.")
  
  int w_len = this->nuc_ids.size();
  int start = min_exclude - w_len;
  if (start < 0) start = 0;

  for (auto i = min_exclude; i < max_exclude; i++) {
    this->allowed[i] = false;
  }
}

std::vector<std::vector<int> > WindowSpec::get_constraints() const {
  std::vector<std::vector<int> > rval;
  int w_len = nuc_ids.size();
  
  for (auto i = 0; i < allowed.size(); i++) {
    std::vector<int> cval;
    if (allowed[i]) {
      auto it = source_nucs.begin() + i;
      cval.assign(it, it + w_len);
      rval.push_back(cval);
    }
  }
  return rval;
}

std::vector<int> DomainListSpec::get_nuc_ids_helper(const std::vector<std::string> & dnames, 
    const SequenceSpec & seqspec) const {
  std::vector<int> rval;
  for (auto & dname : dnames) {
    std::vector<int> cur_nuc_ids = seqspec.get_element(dname).get_nuc_ids();
    append(rval, cur_nuc_ids);
  }
  return rval;
}

std::ostream & param_format(std::ostream & out) {
  out << std::setw(7);
  return out;
};

std::ostream & bool_format(std::ostream & out) {
  out << std::boolalpha;
  return out;
}

std::ostream & longflt_format(std::ostream & out) {
  out << std::setw(14) << std::setprecision(10) << std::fixed;
  return out;
}

std::ostream & flt_format(std::ostream & out) {
  out << std::setw(10); out << std::setprecision(6); out << std::fixed;
  return out;
};

std::ostream & exp_format(std::ostream & out) {
  out << std::setw(10) << std::setprecision(6) << std::scientific;
  return out;
};

/*
 * This format is used for partition functions and concentrations,
 * anything that needs to be a scale-free parameter
 */
std::ostream & scale_format(std::ostream & out) {
  out << std::setw(6) << std::setprecision(3) << std::fixed;
  return out;
};

}
