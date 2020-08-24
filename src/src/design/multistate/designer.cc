#include "designer.h"
#include "algorithms.h"

#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <iomanip>

#include <cfloat>
#include <cmath>

namespace nupack {
// Code copied from optimize() and optimize_tubes(). Need to refactor.
void Designer::evaluate_defect() {
  results.clear();
  for (auto i = 0; i < spec.eval.options.N_population; i++) {
    results.emplace_back();
    results.back().init_random(spec);
  }
  
  PhysicalSpec tubespec = spec.eval.physical.get_depth(0, true);
  
  evaluate(tubespec, results);
  pick_best_par(results);
}
  
void Designer::optimize() {
  // Initialize population of results
  results.clear();
  for (auto i = 0; i < spec.eval.options.N_population; i++) 
    results.emplace_back();
  
  optimize_tubes();
}

void Designer::optimize_tubes() {
  auto temp_results = results;
  auto tmp = spec.constraints.init_random(); // unused, but advances RNG
  
  PhysicalSpec tubespec = spec.eval.physical.get_depth(0, true);
  int ref_ind = 0;

  // init random
  for (auto & tr : temp_results) tr.init_random(spec);

  spec.eval.physical.decompose(spec.eval.options);

  // optimize trees
  optimize_forest(spec.eval.physical, temp_results);

  // Evaluate at tube level
  auto tube_res = temp_results;
  evaluate(tubespec, tube_res);
  results = tube_res;

  while (!tubes_satisfied(tubespec, tube_res, temp_results) && !spec.eval.options.opt_time_elapsed()) {
    // Add worst offenders
    NUPACK_DEBUG("Refocusing");
    if (spec.eval.options.print_steps >= PRINT_REFOCUS) {
      print_leaf_root(tubespec, tube_res, "refocus", ref_ind);
      ref_ind++;
    }
    add_off_targets(spec.eval.physical, tube_res, temp_results);

    optimize_forest(spec.eval.physical, temp_results);
    
    tube_res = temp_results;
    evaluate(tubespec, tube_res);
    append(results, tube_res);

    pick_best_par(results);
  }
}

void Designer::copy_sequences(Results & a, Results & b) {
  NUPACK_DEBUG("a.size() " + to_string(a.size()) + " b.size() " + to_string(b.size()));
  for (auto it_a = a.begin(), it_b = b.begin(); it_a != a.end(); ++it_a, ++it_b)
      it_a->eval.sequences = it_b->eval.sequences;
}

void Designer::optimize_forest(PhysicalSpec & phys_spec, Results & res) {
  std::vector<PhysicalSpec> specs;
  std::vector<Results> results;

  int n_levels = phys_spec.get_max_depth() + 1;
  for (auto i = 0; i < n_levels; i++) {
    specs.push_back(phys_spec.get_depth(i));
    results.push_back(Results());
  }

  DesignResult tmpres;
  for (auto & r : res) {
    tmpres.eval.sequences = r.eval.sequences;
    results[n_levels - 1].push_back(tmpres);
  }

  bool root_accepted = false;
  static int red_ind = 0;
  while (!root_accepted) {
    root_accepted = true;

    optimize_leaves(specs[n_levels-1], results[n_levels-1]);

    for (auto i = n_levels - 2; i >= 0 && root_accepted; i--) {
      // Merge sequences
      Results tmp_result = results[i];
      tmp_result.resize(results[i + 1].size());      
      copy_sequences(tmp_result, results[i + 1]);

      evaluate(specs[i], tmp_result);

      append(results[i], tmp_result);

      pick_best_par(results[i]);
      copy_sequences(results[i + 1], results[i]);

      evaluate(specs[i + 1], results[i + 1]);

#ifndef NDEBUG
      std::stringstream ss;
      ss << "PE " << std::setw(4) << i << ":";
      serialize_defects(std::cout, specs[i], results[i], 0, ss.str());
#endif // NDEBUG

      if (!parent_satisfied(specs[i], tmp_result, results[i + 1])
          && !spec.eval.options.opt_time_elapsed()) {
        NUPACK_DEBUG("Redecomposing");
        root_accepted = false;
        
        if (spec.eval.options.print_steps >= PRINT_REDECOMPOSE) {
          print_leaf_root(specs[i], tmp_result, "redecompose", red_ind);
          red_ind++;
        }
        redecompose(phys_spec, results[i], results[i + 1]);
        
        for (auto j = i + 1; j < n_levels; j++) {
          PhysicalSpec tmp = phys_spec.get_depth(j);
          if (j < specs.size()) {
            specs[j] = tmp;
            results[j] = Results();
          } else {
            specs.push_back(tmp);
            results.emplace_back();
          }

          results[n_levels - 1].resize(results[i].size());
          copy_sequences(results[n_levels - 1], results[i]);
        }
      }
    }
  }
  res = results[0];
}

void Designer::print_leaf_root(PhysicalSpec & phys_spec, Results & res, std::string type, int index) {
  std::stringstream ss;
  ss << spec.eval.options.file_prefix << "_" << type << "_" << index << ".npo";
  std::string new_filename = ss.str();

  NUPACK_DEBUG("Printing leaves: " << type << " " << index);

  std::fstream npofile(new_filename, std::fstream::out | std::fstream::trunc);
  spec.eval.options.serialize(npofile, 0, "");
  res[0].eval.sequences.serialize(spec.eval.sequences, spec.eval.options, npofile, 0, "");
  res[0].eval.physical.serialize(phys_spec, spec.eval.options, npofile, 0, "");
  spec.objectives.serialize(spec.eval, res[0].eval, npofile, 0, "");
}

void Designer::optimize_leaves(PhysicalSpec & phys_spec, Results & res) {
  auto & opts = spec.eval.options;
  mutate_leaves(phys_spec, res);

  bool all_sat = all_satisfied(phys_spec, res);
  int m = 0;
  static int i_opt = 0;
  int cur_reopt = 0;

  // std::stringstream sufss;
  // sufss << "_" << i_opt;
  // std::string suffix(sufss.str());
  std::string suffix = "_" + to_string(i_opt);
  
// #ifndef NDEBUG
//   std::stringstream ss;
//   ss << "[DEBUG] L1O" << std::setw(4) << m << " / " << std::setw(4) << 
//       opts.M_reopt << ": ";
//   serialize_defects(std::cout, phys_spec, res, 0, ss.str());
// #endif // NDEBUG

  if (opts.print_steps >= PRINT_REOPTIMIZE) {
    print_leaf_root(phys_spec, res, "leafropt" + suffix, cur_reopt);
    cur_reopt++;
  }

  while (!all_sat && m < opts.M_reopt && !opts.opt_time_elapsed()) {
    NUPACK_DEBUG("Reoptimizing leaves");
    Results new_gen = res;
    make_offspring(phys_spec, new_gen, opts.M_reseed);

    if (opts.print_steps >= PRINT_RESEED) {
      evaluate(phys_spec, new_gen);
      print_leaf_root(phys_spec, new_gen, "leafreseed" + suffix, cur_reopt);
    }

    mutate_leaves(phys_spec, new_gen);
    evaluate(phys_spec, new_gen);

    if (opts.print_steps >= PRINT_REOPTIMIZE) {
      evaluate(phys_spec, new_gen);
      print_leaf_root(phys_spec, new_gen, "leafreopt" + suffix, cur_reopt);
      cur_reopt++;
    }

    append(res, new_gen);

    if (pick_best_par(res)) {
// #ifndef NDEBUG
//       std::stringstream ss;
//       ss << "[DEBUG] LRO" << std::setw(4) << m << " / " << std::setw(4) << 
//           opts.M_reopt << ": ";
//       serialize_defects(std::cout, phys_spec, res, 0, ss.str());
// #endif // NDEBUG
      m = 0;

      all_sat = all_satisfied(phys_spec, res);
    } else {
      m++;
// #ifndef NDEBUG
//       std::stringstream ss;
//       ss << "[DEBUG] REO" << std::setw(4) << m << " / " << std::setw(4) << 
//           opts.M_reopt << ": ";
//       serialize_defects(std::cout, phys_spec, res, 0, ss.str());
// #endif // NDEBUG
    }
  }
  i_opt++;
}

bool Designer::parent_satisfied(PhysicalSpec & curspec, Results & parent, Results & child) {
  NUPACK_CHECK(parent.size() == child.size(), "Parent and child have different lengths");

  for (auto par_it = parent.begin(), chi_it = child.begin(); 
      par_it != parent.end(); par_it++, chi_it++) {
    auto parent_ob = spec.objectives.get_defects(spec.eval, par_it->eval);
    auto child_ob = spec.objectives.get_defects(spec.eval, chi_it->eval);
    auto parent_sat = spec.objectives.satisfied(spec.eval, par_it->eval);
    for (auto i_ob = 0; i_ob < parent_ob.size(); i_ob++) {
      if (!parent_sat[i_ob] && 
          parent_ob[i_ob] > (child_ob[i_ob] / spec.eval.options.f_stringent)) {
        return false;
      }
    }
  }
  return true;
}

bool Designer::tubes_satisfied(PhysicalSpec & curspec, Results & tube_res, Results & tree_res) {
  NUPACK_CHECK(tube_res.size() == tree_res.size(), "Tube and tree have different lengths");

  for(auto tube_it = tube_res.begin(), tree_it = tree_res.begin();
     tube_it != tube_res.end(); tube_it++, tree_it++) {
    auto tube_ob = spec.objectives.get_defects(spec.eval, tube_it->eval);
    auto tree_ob = spec.objectives.get_defects(spec.eval, tree_it->eval);
    auto tube_sat = spec.objectives.satisfied(spec.eval, tube_it->eval);
    for (auto i_ob = 0; i_ob < tube_ob.size(); i_ob++) {
      if (!tube_sat[i_ob] && 
          tube_ob[i_ob] > (tree_ob[i_ob] / spec.eval.options.f_stringent)) {
        return false;
      }
    }
  }
  return true;
}

void Designer::add_off_targets(PhysicalSpec & curspec, Results & tube_res, Results & tree_res) {
  NUPACK_CHECK(tube_res.size() == tree_res.size(), "Tube and tree have different lengths");

  std::vector<std::vector<DBL_TYPE> > tree_stop_ob(tree_res.size());

  auto i_res = 0;
  for (auto tube_it = tube_res.begin(), tree_it = tree_res.begin(); 
      tube_it != tube_res.end(); ++tube_it, ++tree_it, ++i_res) {
    auto tube_ob = spec.objectives.get_defects(spec.eval, tube_it->eval);
    auto tree_ob = spec.objectives.get_defects(spec.eval, tree_it->eval);
    auto tube_sat = spec.objectives.satisfied(spec.eval, tube_it->eval);
    for (auto i_ob = 0; i_ob < tube_ob.size(); i_ob++) {
      if (!tube_sat[i_ob] && tree_ob[i_ob] < tube_ob[i_ob]) {
        tree_stop_ob[i_res].push_back(tube_ob[i_ob] - spec.eval.options.f_refocus
            * (tube_ob[i_ob] - tree_ob[i_ob]));
      } else {
        tree_stop_ob[i_res].push_back(0);
      }
    }
  }

  bool satisfied = false;
  while (!satisfied) {
    satisfied = true;
    
    auto fail_tube = tube_res.end();
    auto fail_tree = tree_res.end();
    int fail_ob = -1;
    
    auto i_res = 0;
    for (auto tube_it = tube_res.begin(), tree_it = tree_res.begin();
        tube_it != tube_res.end(); ++tube_it, ++tree_it, ++i_res) {
      auto tree_ob = spec.objectives.get_defects(spec.eval, tree_it->eval);
      auto tube_ob = spec.objectives.get_defects(spec.eval, tube_it->eval);
      for (auto i_ob = 0; i_ob < tree_ob.size() && satisfied; i_ob++) {
        if (tree_ob[i_ob] < tree_stop_ob[i_res][i_ob]) {
          NUPACK_DEBUG("Tube failed " << i_res << " " << i_ob << " tube: " << 
              tube_ob[i_ob] << "   stop: " << tree_stop_ob[i_res][i_ob] << "   tree:" << tree_ob[i_ob]);
          satisfied = false;
          fail_tube = tube_it;
          fail_tree = tree_it;
          fail_ob = i_ob;
        }
      }
    }

    if (!satisfied) {
      int phys_id = spec.objectives.get_physical_id(fail_ob);
      auto tube_id = spec.objectives.get_tube_id(fail_ob);
      NUPACK_CHECK(!tube_id.empty(),
          "Invalid objective type failed during focus check. "
          "Must be a tube or combined objective");
      int max_struc = -1;
      DBL_TYPE max_conc = 0.0;
      
      for (auto & c_tube : tube_id) {
        NUPACK_CHECK(c_tube >= 0, "Invalid tube when adding off-targets. "
            "Objective type mismatch");
        auto & cur_tube = fail_tube->eval.physical.get_results()[phys_id]
            .get_tubes()[c_tube];
        auto tube_conc = cur_tube.get_concentrations();
        auto tot_nuc_conc = cur_tube.get_nuc_conc();
        
        auto child_conc = fail_tree->eval.physical.get_results()[phys_id]
            .get_tubes()[c_tube].get_concentrations();
        
        auto & cur_param_spec = curspec.get_param_specs()[phys_id];
        auto complexes = cur_param_spec.get_tubes()[c_tube].get_complexes();
        auto order_to_struc = cur_param_spec.get_ord_struc_map();
        
        NUPACK_CHECK(tube_conc.size() == child_conc.size() && 
            tube_conc.size() == complexes.size(), "Structure size mismatch");
        
        for (auto i_struc = 0; i_struc < tube_conc.size(); i_struc++) {
          auto c_struc = complexes[i_struc].order_ind;
          if (order_to_struc[c_struc] < 0) {
            auto struc_nuc_conc = fail_tube->eval.physical.get_results()[phys_id]
              .get_orders()[c_struc].size() * tube_conc[i_struc];
            if (struc_nuc_conc / tot_nuc_conc > max_conc) {
              max_conc = struc_nuc_conc / tot_nuc_conc;
              max_struc = c_struc;
            }
          }
        }
      }

      curspec.add_off_target(phys_id, max_struc, spec.eval.sequences);
      max_struc = curspec.get_param_specs()[phys_id].get_strucs().size() - 1;

      evaluate(curspec, tree_res);
      curspec.decompose_ppair(phys_id, max_struc, 0, fail_tree->eval.physical,
          fail_tree->eval.sequences, spec.eval.options);

      PhysicalSpec treespec = curspec.get_depth(0, false);
      evaluate(treespec, tree_res);      
    }
  }
}

void Designer::redecompose(PhysicalSpec & phys_spec, Results & parent, 
    Results & child) {
  NUPACK_CHECK(parent.size() == child.size(), "Parent child size mismatch");
  
  int child_depth = child.front().eval.physical.get_max_depth();
  std::vector<std::vector<DBL_TYPE> > child_stop_ob(parent.size());
  
  auto & opts = spec.eval.options;
  EvalSpec tmpspec = spec.eval;
  tmpspec.physical = phys_spec;

  auto par_it = parent.begin();
  auto chi_it = child.begin();
  int i_res = 0;
  for ( ; par_it != parent.end(); ++par_it, ++chi_it, ++i_res) {
    auto par_ob = spec.objectives.get_defects(tmpspec, par_it->eval);
    auto chi_ob = spec.objectives.get_defects(tmpspec, chi_it->eval);
    auto tube_sat = spec.objectives.satisfied(tmpspec, par_it->eval);

    for (auto i_ob = 0; i_ob < par_ob.size(); i_ob++) {
      if (!tube_sat[i_ob] && chi_ob[i_ob] / opts.f_stringent < par_ob[i_ob]) {
        // Solve ~C_d - ~C^*_{d+1} / f_stringent <= (f_redecomp * (~C_d - ~C_{d+1} / f_stringent))
        // for ~C^*_{d+1}
        // ~C^*_{d+1} >= (~C_d - (f_redecomp * (~C_d - ~C_{d+1}))) * f_stringent
        child_stop_ob[i_res].push_back(
            (par_ob[i_ob] - (opts.f_redecomp * (par_ob[i_ob] - chi_ob[i_ob] /
            opts.f_stringent))) * opts.f_stringent);
      } else {
        child_stop_ob[i_res].push_back(0);
      }
    }
  }

  bool satisfied = false;
  while (!satisfied) {
    satisfied = true;
    tmpspec.physical = phys_spec;
    auto fail_parent = parent.end();
    auto fail_child = child.end();
    int fail_ob = -1;
    for (par_it = parent.begin(), chi_it = child.begin(), i_res = 0;
        par_it != parent.end(); par_it++, chi_it++, i_res++) {
      std::vector<DBL_TYPE> tree_ob = this->spec.objectives.get_defects(
        tmpspec, chi_it->eval);
      for (auto i_ob = 0; i_ob < tree_ob.size() && satisfied; i_ob++) {
        if (tree_ob[i_ob] < child_stop_ob[i_res][i_ob]) {
          NUPACK_DEBUG("Tube failed " << i_res << " " << i_ob << " " << 
              child_stop_ob[i_res][i_ob] << " / " << tree_ob[i_ob]);
          satisfied = false;
          fail_parent = par_it;
          fail_child = chi_it;
          fail_ob = i_ob;
        }
      }
    }

    if (!satisfied) {
      NUPACK_DEBUG("Redecomposing based ob: " << fail_ob);
      int min_spec = -1;
      int min_struc = -1;
      int min_node = -1;
      DBL_TYPE min_defect = 1e100;

      const DesignResult & temp_parent = *fail_parent;
      const DesignResult & temp_child = *fail_child;
      std::vector<DBL_TYPE> merged_ob;

      for (auto i_spec = 0; i_spec < temp_parent.eval.physical.get_results().size(); i_spec++) {
        const std::vector<StructureResult> & strucs = 
          temp_parent.eval.physical.get_results()[i_spec].get_strucs();
        for (auto i_struc = 0; i_struc < strucs.size(); i_struc++) {
          int n_leaves = strucs[i_struc].get_n_leaves();
          NUPACK_DEBUG("N leaves: " << n_leaves);
          for (auto i_node = 0; i_node < n_leaves; i_node++) {
            DesignResult merged = temp_parent;
            merged.eval.physical.replace_node(temp_child.eval.physical, phys_spec, i_spec, 
                i_struc, i_node,
                this->spec.eval.options);
            merged_ob = this->spec.objectives.get_defects(
              this->spec.eval, merged.eval);
            NUPACK_DEBUG("Objective: " << merged_ob[fail_ob]);
            if (merged_ob[fail_ob] < min_defect) {
              min_defect = merged_ob[fail_ob];
              min_spec = i_spec;
              min_struc = i_struc;
              min_node = i_node;
            }
          }
        }
      }

      phys_spec.decompose_ppair(min_spec, min_struc, min_node, temp_parent.eval.physical,
          temp_parent.eval.sequences, this->spec.eval.options);

      NUPACK_DEBUG("Redecompose at: " << min_spec << " " << min_struc << " " << min_node);
      PhysicalSpec childspec = phys_spec.get_depth(child_depth);
      this->evaluate(childspec, child);
    }
  }
}

void Designer::mutate_leaves(PhysicalSpec & phys_spec, Results & res) {
  NUPACK_DEBUG("Mutating leaves");
  
  auto & opts = spec.eval.options;
  int m = 0;
  evaluate(phys_spec, res);
  bool all_sat = all_satisfied(phys_spec, res);
  for (auto & r : res) r.clear_tabu();

  while (!all_sat && m < opts.M_bad && !opts.opt_time_elapsed()) {
    auto new_gen = res;
    make_offspring(phys_spec, new_gen);
    update_tabu(res, new_gen);
    evaluate(phys_spec, new_gen);
    append(res, new_gen);
    
    if (pick_best_par(res)) {
#ifndef NDEBUG
      std::stringstream ss;
      ss << "[DEBUG] LM " << std::setw(4) << m << " / " << std::setw(4) << 
          opts.M_bad << ": ";
      serialize_defects(std::cout, phys_spec, res, 0, ss.str());
#endif
      
      all_sat = all_satisfied(phys_spec, res);
      m = 0;
    } else {
      m++;
    }
  }
  for (auto & r : res) r.clear_tabu();
}

void Designer::update_tabu(Results & old_gen, Results & new_gen) {
  NUPACK_CHECK(old_gen.size() == new_gen.size(), "Invalid new gen size in update_tabu");

  std::vector<Results::iterator> to_clear;

  auto new_it = new_gen.begin();
  auto old_it = old_gen.begin();
  for ( ; new_it != new_gen.end(); ++new_it, ++old_it) {
    auto & new_vars = new_it->eval.sequences.get_variables();
    if (old_it->is_tabu(new_vars)) {
      to_clear.push_back(new_it);
    } else {
      old_it->add_tabu(new_vars);
    }
  }

  for (auto & cl : to_clear) new_gen.erase(cl);
}

bool Designer::all_satisfied(PhysicalSpec & phys_spec, Results & res) {
  EvalSpec tmp = spec.eval;
  tmp.physical = phys_spec;
  using el = typename std::remove_reference<decltype(res)>::type::value_type;
  return all_of(res, [&](const el & r) {return spec.objectives.all_satisfied(tmp, r.eval); });
}

// Can also make offspring by recombination?
// Can bias offspring to fix correct objectives?
void Designer::make_offspring(PhysicalSpec & phys_spec, Results & res, int n) {
  std::vector<int> mut_vars;
  EvalSpec tmp = spec.eval;
  tmp.physical = phys_spec;

  for (auto & r : res) {
    std::vector<DBL_TYPE> cur_weights(spec.eval.sequences.get_variables().size(), 0);
    spec.objectives.get_variable_mut_weights(tmp, r.eval, cur_weights);

    mut_vars.clear();
    for (auto i = 0; i < n; i++) {
      int samp = sample_weighted_int(cur_weights);
      if (samp >= 0 && cur_weights[samp] > 0) {
        mut_vars.push_back(samp);
        cur_weights[samp] = 0;
      }
    }
    const std::vector<int> & cur_vars = r.eval.sequences.get_variables();
    std::vector<int> new_vars = spec.constraints.make_mutation(mut_vars, cur_vars);
    if (new_vars[0] != -1) 
      r.eval.sequences.set_variables(new_vars, spec.eval.sequences);
  }
}

bool Designer::pick_best_par(Results & res) {
  bool ret = false;

  auto stops = spec.objectives.get_stops(spec.eval, res.begin()->eval);
  std::vector<DBL_TYPE> summaries;

  for (auto & r : res) {
    const auto & def = r.objectives;
    const auto & sat = r.satisfied;
    DBL_TYPE summary = 0;
    for (auto i = 0; i < def.size(); i++) {
      if (sat[i]) {
        summary += 1.0;
      } else {
        NUPACK_DEBUG_CHECK(!(def[i] < stops[i]), 
            "def / stop < 1.0 in unsatisfied objective");
        summary += def[i] / stops[i];
      }
    }
    summaries.push_back(summary);
  }
  auto unsorted_summaries = summaries;
  std::sort(summaries.begin(), summaries.end());

  int n_to_include = spec.eval.options.N_population;
  DBL_TYPE cutoff = summaries.back();
  if (n_to_include < summaries.size()) 
    cutoff = (summaries[n_to_include] + summaries[n_to_include - 1]) / 2;
  cutoff *= (1.0 + 1e-10);

  std::vector<Results::iterator> it_vec;

  auto it = res.begin();
  for ( ; it != res.end(); ++it) {
    auto i = it - res.begin();

    int n_left = res.size() - i;
    int n_to_delete = res.size() - n_to_include - it_vec.size();

    if (n_to_delete > 0 && 
        (unsorted_summaries[i] > cutoff || n_left == n_to_delete)) {
      it_vec.push_back(it);
      if (i < n_to_include) ret = true; 
    } else {
      if (i >= n_to_include) it->clear_tabu();
    }
  }
  
  for (auto & it : it_vec) res.erase(it);

  return ret;
}


void Designer::serialize_defects(std::ostream & out, PhysicalSpec & phys_spec, 
    Results & res, int indent, std::string prefix) {
  EvalSpec tmp = spec.eval;
  tmp.physical = phys_spec;
  
  for (auto it = res.begin(); it != res.end(); ++it) {
    auto i_res = it - res.begin();
    std::stringstream ss;
    ss << prefix << std::setw(4) << i_res << ": ";
    spec.objectives.short_serial(tmp, it->eval, std::cout, 4, ss.str());
  }
}

void Designer::evaluate(PhysicalSpec & phys_spec, Results & res) {
  EvalSpec tmp = spec.eval;
  tmp.physical = phys_spec;
  
  for (auto & r : res) {
    r.eval.physical.evaluate(r.eval.sequences, phys_spec, tmp.options);
    r.eval.symmetries.evaluate(r.eval.sequences, tmp.symmetries);

    r.objectives = spec.objectives.get_defects(tmp, r.eval);
    r.satisfied = spec.objectives.satisfied(tmp, r.eval);
  }
}
}
