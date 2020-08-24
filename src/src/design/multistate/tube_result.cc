#include "tube_result.h"

#include "structure_result.h"
#include "complex_result.h"

#include "design_debug.h"
#include "pathway_utils.h"

#include "equilibrium_concentrations.h"

#include <fstream>

namespace nupack {  
void TubeResult::evaluate(const std::vector<StructureResult> & strucs, 
    const std::vector<OrderResult> & orders, const std::vector<int> & ord_struc_map,
    const TubeSpec & tube, const PhysicalParams & params,
    const NupackInvariants & opts) {
  using sp = ConcSolverParams;
  
  clear_state();
  
  const auto & complexes = tube.get_complexes();
  const int tube_id = tube.get_id();
  
  std::vector<DBL_TYPE> pfuncs;
  std::vector<int> included_ids;
  std::vector<int> included_local_inds;
  std::vector<DBL_TYPE> cur_target_x;
  std::vector<std::vector<int> > strands;
  std::vector<bool> include_strand;
  std::vector<int> strand_map;
  std::vector<int> rev_strand_map;

  std::vector<double> A;
  std::vector<double> x0;
  std::vector<double> dG;
  std::vector<double> x;

  double deflate = 1.0;
  
  int n_strucs = complexes.size();
  for (auto i_ord = 0; i_ord < n_strucs; i_ord++) {
    const auto & comp = complexes[i_ord];
    auto c_ord = comp.order_ind;
    NUPACK_CHECK(c_ord >= 0 && c_ord < ord_struc_map.size(),
        to_string(c_ord) + " is an invalid structure index");
    
    if (ord_struc_map[c_ord] >= 0) {
      DBL_TYPE cur_pfunc = strucs[ord_struc_map[c_ord]].get_pfunc(params);
      if (cur_pfunc > sp::min_delta) {
        pfuncs.push_back(strucs[ord_struc_map[c_ord]].get_pfunc(params));
        included_ids.push_back(c_ord);
        included_local_inds.push_back(i_ord);
        
        DBL_TYPE cur_x = 0.0;
        for (auto conc : comp.target_concs) cur_x += conc;
        cur_target_x.push_back(cur_x);
        
        strands.push_back(strucs[ord_struc_map[c_ord]].get_strands());
      }
    } else if (orders[c_ord].get_evaluated()) {
      DBL_TYPE cur_pfunc = orders[c_ord].get_pfunc();
      if (cur_pfunc > sp::min_delta) {
        pfuncs.push_back(cur_pfunc);
        included_ids.push_back(c_ord);
        included_local_inds.push_back(i_ord);

        DBL_TYPE cur_x = 0.0;
        for (auto conc : comp.target_concs) cur_x += conc;
        cur_target_x.push_back(cur_x);

        strands.push_back(orders[c_ord].get_strands());
      }
    } else {
      deflate = 1.0 - tube.get_passive();
    }
  }

  int m_strucs = pfuncs.size();
  for (auto i_struc = 0; i_struc < m_strucs; i_struc++) {
    for (auto i_strand = 0; i_strand < strands[i_struc].size(); i_strand++) {
      if (strands[i_struc][i_strand] >= (int)include_strand.size()) {
        include_strand.resize(strands[i_struc][i_strand] + 1, false);
      }
      include_strand[strands[i_struc][i_strand]] = true;
    }
  }

  int j_strand = 0;
  int n_strands = include_strand.size();
  strand_map.resize(include_strand.size());
  for (auto i_strand = 0; i_strand < n_strands; i_strand++) {
    if (include_strand[i_strand]) {
      strand_map[i_strand] = j_strand;
      rev_strand_map.push_back(i_strand);
      j_strand++;
    } else {
      strand_map[i_strand] = -1;
    }
  }

  int m_strands = rev_strand_map.size();

  // Fill out the evaluation arrays
  A.resize(rev_strand_map.size() * m_strucs, 0);

  for (auto i_struc = 0; i_struc < m_strucs; i_struc++) {
    dG.push_back((double) -LOG_FUNC(pfuncs[i_struc]));
    x0.push_back((double) cur_target_x[i_struc] * deflate);
    for (auto i_strand = 0; i_strand < strands[i_struc].size(); i_strand++) {
      j_strand = strand_map[strands[i_struc][i_strand]];
      A[j_strand * m_strucs + i_struc] += 1;
    }
  }

  x.resize(m_strucs);

  calc_conc_from_free_energies(x.data(), A.data(), dG.data(), x0.data(), 
        m_strands, m_strucs, sp::n_points, sp::max_iters, sp::tol, sp::delta_bar,
        sp::eta, sp::min_delta, sp::max_trial, sp::perturb_scale, sp::quiet,
        sp::write_log_file, sp::log_file, sp::seed, NULL);
  
  // Remap to full concentration std::vector 
  this->x.assign(n_strucs, 0);
  for (auto i_struc = 0; i_struc < m_strucs; i_struc++) {
    this->x[included_local_inds[i_struc]] = x[i_struc];
  }

  std::vector<std::vector<DBL_TYPE>> target_x;
  for (auto & c : complexes) target_x.push_back(c.target_concs);
  this->target_x = target_x;
  

  // Calculate nodal and nucleotide defects
  for (auto i_ord = 0; i_ord < n_strucs; i_ord++) {
    auto & comp = complexes[i_ord];
    for (auto i_tar = 0; i_tar < comp.target_concs.size(); i_tar++) {
      auto c_ord = comp.order_ind;
      auto c_tar = comp.target_inds[i_tar];
      auto c_struc = ord_struc_map[c_ord];
      
      NUPACK_CHECK(c_struc >= 0, "on-target marked as off-target");

      DBL_TYPE tx = target_x[i_ord][i_tar];
      DBL_TYPE ax = this->x[i_ord];
      DBL_TYPE x_defect = tx - ax;
      DBL_TYPE x_act = ax;
      if (x_defect < 0) {
        x_defect = 0;
        x_act = tx;
      }
      
      DBL_TYPE struc_defect = strucs[c_struc].get_structural_defect(c_tar);
      this->nuc_conc += strucs[c_struc].size() * tx;
      this->defect += (x_act * struc_defect + (strucs[c_struc].size() * x_defect));
      
      const auto & nuc_defects = strucs[c_struc].get_nuc_defects(c_tar);
      for (auto & nd : nuc_defects) {
        std::vector<int> ind(1, tube_id);
        append(ind, nd.first);
        this->nucleotide_defects[ind] += (x_act * nd.second + x_defect);
      }
    }
  }
}

void TubeResult::clear_state() {
  defect = 0.0;
  nuc_conc = 0.0;
  nucleotide_defects.clear();
}

void TubeResult::serialize(const TubeSpec & tube, const std::vector<StructureSpec> & strucspecs,
    const std::vector<int> & ord_struc_map, const PhysicalParams & params,
    const NupackInvariants & invars, std::ostream & out, int indent, std::string prefix) const {
  const auto & complexes = tube.get_complexes();
  std::string ind_str(indent, ' ');
  std::string name_prefix(indent, ' ');

  ind_str = prefix + ind_str;
  if (indent >= 2) name_prefix = prefix + std::string(indent-2, ' ') + "- ";
  DBL_TYPE water_conc = water_density(invars.temperature - ZERO_C_IN_KELVIN);

  out << name_prefix << "name: " << tube.get_name() << std::endl;
  out << ind_str << "defect[M nt]: " << exp_format 
      << this->get_defect() * water_conc << std::endl;
  out << ind_str << "nucleotide conc[M nt]: " << exp_format 
      << this->nuc_conc * water_conc << std::endl;
  out << ind_str << "normalized defect: " << flt_format
      << this->get_defect() / this->nuc_conc << std::endl;
  out << ind_str << "complexes: " << std::endl;
  
  for (auto i_ord = 0; i_ord < complexes.size(); i_ord++) {
    auto i_struc = ord_struc_map[complexes[i_ord].order_ind];
    if (i_struc >= 0) {
      const std::string & name = strucspecs[i_struc].get_name();

      out << ind_str << "  - name: " << name << std::endl;
      out << ind_str << "    concentration: " << exp_format
          << (water_conc * this->x[i_ord]) << std::endl;

      DBL_TYPE target_x = 0;
      for (auto i = 0; i < this->target_x[i_ord].size(); i++) {
        target_x += this->target_x[i_ord][i];
      }
      out << ind_str << "    target concentration: " << exp_format << 
        (water_conc * target_x) << std::endl;
    } 
  }
}
  
}
