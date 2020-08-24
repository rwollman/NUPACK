#include "pathway_input.h"

#include "parsestruc.h"
#include "pathway_parser.h"
#include "pathway_lexer.h"

#include "symmetry_objective.h"
#include "objective_handler.h"
#include "constraint_handler.h"
#include "physical_spec.h"
#include "design_spec.h"
#include "weight_spec.h"
#include "sequence_utils.h"

#include "design_debug.h"
#include "utils.h"
#include "algorithms.h"
#include "types.h"

#include "sys/time.h"
#include <typeinfo>
#include <array>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <list>
#include <utility>
#include <stdio.h>

namespace nupack {

DBL_TYPE ScriptProcessor::get_frac_units(const std::string & unit_name) {
  // Default for fractional types is percent  
  const std::vector<std::pair<std::string, DBL_TYPE>> units = {{"frac", 1.0}, {"%", 0.01}};

  DBL_TYPE multiplier = 1.0;
  for (auto & u : units) {
    if (unit_name == u.first) multiplier = u.second;
  }
  return multiplier;
}

DBL_TYPE ScriptProcessor::process_frac_units(value_struc_t * cur, DBL_TYPE baseval) {
  const std::vector<std::string> units = {"frac", "%"};
  
  if (cur) {
    DBL_TYPE multiplier = get_frac_units(cur->strval);
    for (auto & u : units) {
      if (cur->strval == u) return baseval * multiplier;
    }
  }
  else {
    return baseval;
  }
  NUPACK_ERROR("Invalid unit " + to_string(cur->strval));
}

DBL_TYPE ScriptProcessor::process_conc_units(value_struc_t * cur, DBL_TYPE baseval) {
  const std::vector<std::pair<std::string, DBL_TYPE>> multiplier_map = {
    {"frac", 1.0},
    {"M", 1.0},
    {"dM", 1.0e-1},
    {"cM", 1.0e-2},
    {"mM", 1.0e-3},
    {"uM", 1.0e-6},
    {"nM", 1.0e-9},
    {"pM", 1.0e-12},
    {"fM", 1.0e-15},
    {"aM", 1.0e-18},
    {"zM", 1.0e-21},
    {"yM", 1.0e-24},
  };

  DBL_TYPE water_conc = water_density(invars.temperature - ZERO_C_IN_KELVIN);
  DBL_TYPE multiplier = 1.0 / water_conc;

  if (cur) {
    std::string name;
    for (auto & m : multiplier_map) {
      if (cur->strval == m.first) {
        multiplier = m.second;
        name = m.first;
        break;
      }
    }

    if (name != "frac") multiplier /= water_conc;
  }

  return baseval * multiplier;
}

std::vector<std::string> ScriptProcessor::collect_string_list(value_struc_t * linkedlist) {
  std::vector<std::string> strlist;
  for (auto cur = linkedlist; cur; cur = cur->next) {
    if (cur->strval) {
      strlist.push_back(cur->strval);
    } else {
      break;
    }
  }
  return strlist;
}

NameList ScriptProcessor::collect_complexes(value_struc_t * list) {
  NameList complexes;
  for (auto cur = list; cur; cur = cur->next) 
    complexes.push_back(collect_string_list(cur->list));
  return complexes;
}

void ScriptProcessor::resolve_all(WeightSpec & weights) {
  SequenceSpec & seqs = spec.eval.sequences;

  std::map<std::string, int> strandmap = seqs.get_strand_map();

  resolve_SSM();
  
  for (auto & item : structure_map) {
    auto & cur_struc = item.second;
    if (cur_struc.get_strand_names().size() == 0) {
      /* No domains or strands were added to this structure */
      std::vector<int> breaks = cur_struc.get_breaks();
      std::vector<std::string> names;
      breaks.push_back(cur_struc.size());
      int i_nuc = 0;
      int i_str = 0;
      for (auto br : breaks) {
        std::string cons(br - i_nuc, 'N');
        std::string dname = "-" + cur_struc.get_name() + "_domain_" + to_string(i_str + 1);
        std::string sname = "-" + cur_struc.get_name() + "_strand_" + to_string(i_str + 1);
        
        seqs.add_domain_by_constraints(dname, cons, invars);
        seqs.add_strand_by_domains(sname, std::vector<std::string>(1, dname));
        names.push_back(sname);
        NUPACK_DEBUG("Adding strand " << sname);
        strandmap = seqs.get_strand_map();
        i_nuc = br;
        i_str++;
      }
      cur_struc.set_strand_names(names);
    } else {
      const std::vector<std::string> & strandnames = cur_struc.get_strand_names();
      std::vector<std::string> new_strands;
      bool is_set = false;
      int i_nuc = 0;
      int i_br = 0;
      int n_nucs = cur_struc.size();
      auto breaks = cur_struc.get_breaks();
      std::vector<std::string> cur_domains;
      breaks.push_back(n_nucs);
      int curlen = 0;
      
      for (auto & s_name : strandnames) {
        NUPACK_CHECK(i_nuc < n_nucs,
            "Strand / domain length mismatch in structure: " + cur_struc.get_name());
        
        if (!contains(s_name, strandmap)) {
          // Must be a domain
          is_set = true;
          
          NUPACK_CHECK(i_nuc < n_nucs && i_br < breaks.size() &&
              i_nuc < breaks[i_br] && (cur_domains.size() > 0
              || i_nuc == 0 || (i_br > 0 && breaks[i_br - 1] == i_nuc)),
              "Break in structure " + cur_struc.get_name() +
              " must align with strand/domain break ");
          cur_domains.push_back(s_name);
          
          NUPACK_EXC_CHECK(curlen = seqs.get_domain(s_name).size(),
              "Invalid domain or strand name " + s_name);
          
          i_nuc += curlen;
          
          if (i_nuc == breaks[i_br]) {
            i_br++;
            std::string newname;
            for (auto & dom : cur_domains) newname += dom + "-";
            newname.pop_back();
            newname += " strand";
            
            if (!contains(newname, strandmap)) {
              seqs.add_strand_by_domains(newname, cur_domains);
              strandmap = seqs.get_strand_map();
            }
            new_strands.push_back(newname);
            cur_domains.clear();
          }
        } else {
          NUPACK_EXC_CHECK(curlen = seqs.get_strand(s_name).size(), 
              "Error retrieving strand");
          new_strands.push_back(s_name);
          i_nuc += curlen;
          NUPACK_CHECK(i_br < breaks.size() && i_nuc == breaks[i_br],
              "structure break in " + cur_struc.get_name() +
              " must align with strand or domain breaks "
              + to_string(i_br)
              + " " + to_string(breaks.size()) + " " + to_string(i_nuc)
              + " " + to_string(breaks[i_br]));
          i_br++;
        }
      }
      if (is_set) {
        cur_struc.set_strand_names(new_strands);
      }
    }
  }

  resolve_weights(weights);
  generate_default_parameter_spec();
}

void ScriptProcessor::resolve_SSM() {
  for (auto & item : SSM_map) {
    auto & cur_sym = item.second;
    cur_sym.resolve(spec.eval.sequences);
    spec.eval.symmetries.add_objective(cur_sym);
  }
}

void ScriptProcessor::resolve_weights(WeightSpec & weights) {
  SequenceSpec & seqs = spec.eval.sequences;
  
  std::vector<std::map<std::string, int> > idmap;
  std::map<std::string, int> pmap;
  std::map<std::string, int> tubemap;
  std::map<std::string, int> strucmap;
  std::map<std::string, int> new_strandmap = seqs.get_strand_map();
  std::map<std::string, int> new_domainmap = seqs.get_domain_map();

  int i_tube = 0, i_struc = 0;
  for (auto & item : tube_map) {
    item.second.set_id(i_tube);
    tubemap[item.first] = i_tube;
    i_tube++;
  }
  for (auto & item : structure_map) {
    auto & cur_struc = item.second;
    cur_struc.resolve_strand_names(new_strandmap);
    cur_struc.resolve_nuc_ids(seqs);
    cur_struc.set_id(i_struc);
    strucmap[item.first] = i_struc;
    i_struc++;
  }
  
  idmap.push_back(pmap);
  idmap.push_back(tubemap);
  idmap.push_back(strucmap);
  idmap.push_back(new_strandmap);
  idmap.push_back(new_domainmap);
  weights.resolve_names(idmap);
}

void ScriptProcessor::generate_default_parameter_spec() {
  std::set<std::string> tubes_included;
  std::set<std::string> strucs_included;
  std::map<std::string, std::pair<int, int> > struc_inds_cur;
  
  PhysicalParams params;
  params.temperature = invars.temperature;

  SingleParamSpec default_spec;
  default_spec.set_id(0); // currently only one parameter set per design
  default_spec.set_parameters(params);

  struc_inds_cur.clear();

  int n_strucs_cur = 0;
  for (auto & item : tube_map) {
      TubeSpec & cur_tube = item.second;
      std::string cur_tube_name = cur_tube.get_name();
      if (!contains(cur_tube_name, tubes_included)) {
        // This tube wasn't included, add it to global spec
        auto strucs_in_tube = cur_tube.get_struc_names();
        for (auto & struc : strucs_in_tube) {
          if (!contains(struc, struc_inds_cur)) {
            NUPACK_CHECK(contains(struc, structure_map),
                "Invalid structure: " + struc + " specified in tube: " + cur_tube_name);
            auto & cur_struc = structure_map.find(struc)->second;
            auto tmp_id = default_spec.add_structure(cur_struc);
            struc_inds_cur[struc] = tmp_id;
            strucs_included.insert(struc);
            n_strucs_cur++;
          }
        }

        cur_tube.resolve_structures(struc_inds_cur);
        default_spec.add_tube(cur_tube);
        tubes_included.insert(cur_tube_name);
      }
  }

  for (auto & item : structure_map) {
    auto & cur_struc = item.second;
    std::string cur_struc_name = cur_struc.get_name();
    
    if (!contains(cur_struc_name, strucs_included)) {
      // This structure wasn't included, need to add correctly
      default_spec.add_structure(cur_struc);
      n_strucs_cur++;
    }
  }

  if (n_strucs_cur > 0) {
    default_spec.enumerate_off_targets(spec.eval.sequences.get_strand_map());
    default_spec.compute_target_matrices();
    spec.eval.physical.add_param_spec(default_spec);
  }
}

void ScriptProcessor::add_pattern_constraints(std::vector<Prevent> & prevents) {
  const auto & strands = spec.eval.sequences.get_strands();
  std::vector<int> poss_nucs;
  NUPACK_EXC_CHECK(poss_nucs = spec.constraints.get_possible_nucleotides(),
      "Preliminary constraint satisfaction failed");
  
  for (const auto & p : prevents) {
    std::vector<std::string> cur_targets;
    if (!p.has_specific_targets()) {
      for (const auto & s : strands) cur_targets.push_back(s.get_name());
    } else {
      cur_targets = p.strands;
    } 
    
    for (const auto & pat : p.prevented_sequences)
      for (auto & cur_target : cur_targets)
        spec.constraints.add_constraint(PatternConstraint(cur_target, pat, spec.eval.sequences, poss_nucs));
  }
  
}

void ScriptProcessor::add_library_constraints(LibraryMap lib_defs, LibraryAssignmentMap lib_assignments) {
  for (auto & l : lib_assignments) {
    const auto & cur_name = l.first;
    const auto & cur_lib_name = l.second.first;
    const auto & cur_definitions = l.second.second;

    const auto & curspec = spec.eval.sequences.get_element(cur_name);
    const auto & nuc_ids = curspec.get_nuc_ids();
    const auto & nuc_to_vars = spec.eval.sequences.get_variables();

    NUPACK_CHECK(contains(cur_lib_name, lib_defs), cur_lib_name + " is not a defined library");
    const auto & curlib = lib_defs.find(cur_lib_name)->second;

    // add constraints on each subsequence of domain independently
    int i_nuc = 0;
    for (auto & def : cur_definitions) {
      auto cur_def = curlib.find(def);
      NUPACK_CHECK(contains_it(cur_def, curlib),
          def + " does not exist in library " + cur_lib_name);

      const std::vector<std::string> & cur_defs = cur_def->second;
      NUPACK_CHECK(cur_defs.size() > 0, "no possibilities for current library");

      int cur_len = cur_defs[0].size();

      NUPACK_CHECK(i_nuc + cur_len <= nuc_ids.size(), "Library/element length mismatch");

      std::vector<int> cur_vars;
      for (int j_nuc = 0; j_nuc < cur_len; j_nuc++) {
        cur_vars.push_back(nuc_to_vars[nuc_ids[i_nuc + j_nuc]]);
      }
      
      std::vector<trinary> poss(cur_defs.size(), true);
      int supp_var = spec.constraints.add_variable(poss);

      spec.constraints.add_constraint(WordConstraint(cur_vars, cur_defs, supp_var));

      i_nuc += cur_len;
    }
  }
}

void ScriptProcessor::add_match_constraints() {
  for (auto & item : match_spec_map) {
    auto & matchspec = item.second;
    auto matches = matchspec.make_match_constraints(spec.eval.sequences);
    for (auto & m : matches) spec.constraints.add_constraint(m);
  }
}

void ScriptProcessor::add_implied_constraints() {
  SequenceSpec & seqs = spec.eval.sequences;
  PhysicalSpec & physspec = spec.eval.physical;
  ConstraintHandler & constraints = spec.constraints;
  // Create the nucleotide variables
  
  const auto & nucs = seqs.get_nucleotides();
  for (int i_nuc_id = 0; i_nuc_id < nucs.size(); ++i_nuc_id) {
    int nuc_con = nucs[i_nuc_id];
    int i_var = constraints.add_nucleotide_variable(nuc_con);
    seqs.set_variable(i_nuc_id, i_var);
  }
  const auto & vars = seqs.get_variables();
  const auto & doms = seqs.get_domains();

  for (auto & dom : doms) {
    auto name = dom.get_name();
    if (name.size() > 0 && name[name.size() - 1] != '*') {
      std::string c_name = name + "*";
      int o_ind = seqs.get_domain_id(c_name);
      if (o_ind >= 0 && o_ind < doms.size()) {
        auto ids1 = dom.get_nuc_ids();
        auto ids2 = doms[o_ind].get_nuc_ids();

        NUPACK_CHECK(ids1.size() == ids2.size(), "complement size mismatch");        
        for (auto j = 0; j < ids1.size(); j++) {
          int nid1 = ids1[j];
          int nid2 = ids2[ids1.size() - 1 - j];

          constraints.add_constraint(CompConstraint(nid1, nid2, NUPACK_CS_STRONG));
        }
      }
    }
  }

  // Go through the structures and add implied complementarity constraints
  const auto & specs = physspec.get_param_specs();
  for (auto & spec : specs) {
    const auto & strucs = spec.get_strucs();
    for (auto & struc : strucs) {
      auto nuc_ids = struc.get_nuc_ids();
      auto pairings = struc.get_structures();

      for (auto & pairing : pairings) {
        NUPACK_CHECK(nuc_ids.size() == pairing.size(),
            "Nucleotide sequence length (" + to_string(nuc_ids.size()) +
            ") and structure size (" + to_string(pairing.size())
            + ") are not equal");
        for (auto i_nuc = 0; i_nuc < nuc_ids.size(); i_nuc++) {
          int j_nuc = pairing[i_nuc];
          if (j_nuc > i_nuc) {
            int i_nuc_id = vars[nuc_ids[i_nuc]];
            int j_nuc_id = vars[nuc_ids[j_nuc]];

            if (invars.allow_wobble || invars.allow_mismatch) {
              constraints.add_constraint(
                  CompConstraint(i_nuc_id, j_nuc_id, NUPACK_CS_WEAK));
            } else if (!invars.allow_wobble && !invars.allow_mismatch) {
              constraints.add_constraint(
                  CompConstraint(i_nuc_id, j_nuc_id, NUPACK_CS_STRONG));
            }
          }
        }
      }
    }
  }
}

void ScriptProcessor::set_window_source(WindowSpec & wind, value_struc_t * source) {
  auto source_names = collect_string_list(source);
  
  std::vector<SourceSpec *> sources;  
  for (auto & tsource_name : source_names) {
    NUPACK_CHECK(already_named(tsource_name),
        "Source " + tsource_name + " not found in window " + wind.get_name());
    
    auto it = source_map.find(tsource_name);
    NUPACK_CHECK(contains(tsource_name, source_map),
        tsource_name + " not previously defined as source.");
    
    sources.push_back(&it->second);
  }
  
    wind.set_source(sources);
}

void ScriptProcessor::set_window_similarity(WindowSpec* wind, 
    value_struc_t * unitsource, value_struc_t * range) {
  NUPACK_CHECK(range->next, "Unbalanced range when setting window similarity");
  NUPACK_CHECK(range->next->next == NULL,
      "Only one range in similarity assignment");

  DBL_TYPE min_sim = range->doubleval;
  DBL_TYPE max_sim = range->next->doubleval;
  DBL_TYPE scale = 1.0;
  NUPACK_CHECK(!(max_sim < min_sim), "Maximum similarity must be >= minimum similarity");
  
  auto tmpunit = unitsource;
  NUPACK_CHECK(tmpunit, "No source provided");

  std::string tmpunit_str(tmpunit->strval);

  if (tmpunit_str == "frac" || tmpunit_str == "%") {
    if (tmpunit_str == "%") scale = 0.01;
    tmpunit = tmpunit->next;
    NUPACK_CHECK(tmpunit, "No source provided");
  }

  std::string tsource_name(wind->get_source());
  std::string osource_name(tmpunit->strval);
  // check for valid source included here because possible for user to attempt
  // to define similarity before assigning a source, which is not handled.
  NUPACK_CHECK(already_named(tsource_name),
      "Source " + tsource_name + " not found in window " + wind->get_name());
  NUPACK_CHECK(already_named(osource_name),
      "Source " + osource_name + " not found in window " + wind->get_name());
  
  NUPACK_CHECK(contains(tsource_name, source_map),
      tsource_name + " not previously defined as source.");
  NUPACK_CHECK(contains(osource_name, source_map),
      osource_name + " not previously defined as source.");

  auto & osource = source_map.find(osource_name)->second;
  wind->allow_similar(osource, scale * min_sim, scale * max_sim);
}

void ScriptProcessor::add_exclusions(WindowSpec * wind, value_struc_t * vals) {
  auto cur = vals;
  while (cur) {
    int start = cur->intval;
    cur = cur->next;
    NUPACK_CHECK(cur, "Invalid exclusion list, it must be even");

    int end = cur->intval;
    wind->exclude_range(start, end);
    cur = cur->next;
  }
}

void ScriptProcessor::set_ssm_wordsize(value_struc_t * cur, int line) {
  std::string name(cur->id->strval);
  NUPACK_CHECK(already_named(name), "Setting wordsize for nonexistent objective");
  
  auto it = SSM_map.find(name);
  NUPACK_CHECK(contains_it(it, SSM_map), "Setting wordsize for non SSM entity");
  auto & cur_SSM = it->second;

  std::vector<int> sizes;
  auto currange = cur->def;
  while (currange) {
    int start = round(currange->doubleval);
    currange = currange->next;
    NUPACK_CHECK(currange, "Internal error, single element found in range");
    
    int end = round(currange->doubleval);
    NUPACK_CHECK(start > 0, "Word size must be greater than 0, " +
        to_string(start) + " found on line " + to_string(line));
    NUPACK_CHECK(start <= end, "Start must be less than or equal to end"
        " on line " + to_string(line));

    for (auto i = start; i <= end; i++) sizes.push_back(i);
    currange = currange->next;
  }

  cur_SSM.set_word_len(sizes);
}

void ScriptProcessor::add_window(std::string name, std::vector<std::string> domains, const SequenceSpec & seqs) {
  std::vector<SingleSequenceSpec> seqlist;
  for (auto & dom : domains) {
    NUPACK_EXC_CHECK(
        seqlist.push_back(SingleSequenceSpec(seqs.get_element(dom))),
        "cannot find domain or strand " + dom + " in "
        "definition for window " + name);
  }
  window_map.emplace(name, WindowSpec(name, seqlist));
  names.insert(name);
}

void ScriptProcessor::set_item_weight(value_struc_t * cur, int line, WeightSpec & weights) {
  std::string name(cur->id->strval);
  NUPACK_CHECK(already_named(name), 
      "Weight being set for an undefined entity " + 
      name + " on line " + to_string(line));

  std::vector<std::string> idlist = collect_string_list(cur->id);

  DBL_TYPE curdbl = (DBL_TYPE)cur->def->doubleval;
  curdbl = process_frac_units(cur->stype->next, curdbl);
  weights.add_weight(idlist, curdbl);
}

void ScriptProcessor::add_match_domain(value_struc_t * cur, int line) {
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(already_named(namestring),
      "Setting similarity range on nonexistent similarity sequence " + namestring
      + " on line " + to_string(line));
  
  NUPACK_CHECK(cur->stype && cur->stype->next &&
      cur->stype->next->strval,
      "Invalid units provided");
  NUPACK_CHECK(cur->stype && cur->stype->next
      && cur->stype->next->next && cur->stype->next->next->strval,
      "Invalid name provided on line " + to_string(line));

  std::string dom_name(cur->stype->next->next->strval);
  NUPACK_CHECK(already_named(dom_name),
      "Domain \"" + dom_name + "\" not found on line: " + to_string(line));
  
  std::string unit_type(cur->stype->next->strval);
  DBL_TYPE multiplier = get_frac_units(unit_type);

  std::vector<double> mmin;
  std::vector<double> mmax;
  value_struc_t * ranges = cur->def;
  while (ranges) {
    NUPACK_CHECK(ranges->type == NP_TT_FLOAT || ranges->type == NP_TT_INTEGER,
        "internal error, invalid range type");
    
    double min = ranges->doubleval * multiplier;
    ranges = ranges->next;
    NUPACK_CHECK(ranges, "internal error, ranges must come in pairs");
    
    double max = ranges->doubleval * multiplier;
    NUPACK_CHECK(min > -1e-100 && max < (1.0 + 1e-8),
        "fractional ranges must lie between 0 and 1 (0% and 100%)");
    NUPACK_CHECK(!(min > max), "similarity range specified out of order on line " + to_string(line));

    mmin.push_back(min);
    mmax.push_back(max);
    ranges = ranges->next;
  }

  auto & cur_match = match_spec_map.find(namestring)->second;
  cur_match.add_match(dom_name, mmin, mmax);
}

void ScriptProcessor::gather_independent_entities(std::set<std::string> & tubes_to_add, std::set<std::string> & strucs_to_add) {
  const std::vector<SingleParamSpec> & sets = spec.eval.physical.get_param_specs();
  for (auto & set : sets) {
    const auto & strucs = set.get_strucs();
    const auto & tubes = set.get_tubes();
    for (auto & struc : strucs) {
      std::vector<std::string> names = struc.get_names();
      for (auto & name : names) strucs_to_add.insert(name);
    }
    for (auto & tube : tubes) {
      tubes_to_add.insert(tube.get_name());
      auto & names = tube.get_struc_names();
      for (auto & name : names) {
        strucs_to_add.erase(name);
      }
    }
  }
}

void ScriptProcessor::add_global_objective(DBL_TYPE stop) {
  std::set<std::string> tubes_to_add;
  std::set<std::string> strucs_to_add;
  gather_independent_entities(tubes_to_add, strucs_to_add);

  for (auto & ob : spec.objectives.get_names()) {
    strucs_to_add.erase(ob);
    tubes_to_add.erase(ob);
  }

  std::vector<std::string> element_names(strucs_to_add.begin(), strucs_to_add.end());
  append(element_names, tubes_to_add);
  NUPACK_CHECK(!element_names.empty(),
      "Global stop condition specified, but no elements are unused in other objectives");

  spec.objectives.add_objective(CombinedObjective(element_names, stop));
}

void ScriptProcessor::add_default_stop_conditions() {
  std::set<std::string> tubes_to_add;
  std::set<std::string> strucs_to_add;
  gather_independent_entities(tubes_to_add, strucs_to_add);

  for (auto & struc : strucs_to_add) {
    spec.objectives.add_objective(StructureObjective(struc, NUPACK_DEF_STOP));
    NUPACK_DEBUG("Adding structure objective: " << struc);
  }
  for (auto & tube : tubes_to_add) {
    spec.objectives.add_objective(TubeObjective(tube, NUPACK_DEF_STOP));
    NUPACK_DEBUG("Adding tube objective: " << tube);
  }
}

void ScriptProcessor::process_temperature() { 
  std::string unitstring("C");
  if (cur->stype->next) unitstring.assign(cur->stype->next->strval);
  if (unitstring == "C") {
    invars.temperature = (DBL_TYPE) cur->def->doubleval + ZERO_C_IN_KELVIN;
  } else if (unitstring == "K") {
    invars.temperature = (DBL_TYPE) cur->def->doubleval;
  } else {
    NUPACK_ERROR("Invalid units for temperature.");
  }
}

void ScriptProcessor::process_material() { 
  std::string matstring(cur->def->strval);
  convert_to_lower(matstring);
  if (matstring == "rna" || matstring == "rna1995") {
    invars.material = RNA;
  } else if (matstring == "rna37" || matstring == "rna1999") {
    invars.material = RNA37;
  } else if (matstring == "dna" || matstring == "dna1998") {
    invars.material = DNA;
  } else {
    invars.material = USE_SPECIFIED_PARAMETERS_FILE;
    invars.material_string = matstring;
    strncpy(PARAM_FILE, cur->def->strval, 99);
    PARAM_FILE[99] = '\0';
  }
}

void ScriptProcessor::process_sodium() { 
  DBL_TYPE curdbl = cur->def->doubleval;
  invars.sodium = process_conc_units(cur->stype->next, curdbl);
}

void ScriptProcessor::process_magnesium() { 
  DBL_TYPE curdbl = cur->def->doubleval;
  invars.magnesium = process_conc_units(cur->stype->next, curdbl);
}

DBL_TYPE ScriptProcessor::process_global_stopdef() { 
  DBL_TYPE curdbl = cur->def->doubleval;
  auto global_stop = process_frac_units(cur->stype->next, curdbl);
  NUPACK_DEBUG("Global stop set: " + to_string(global_stop));
  NUPACK_CHECK(global_stop < 1.0, "Invalid global fraction specified, must be < 1.0 or 100 %");
  return global_stop;
}

void ScriptProcessor::process_dangles() { 
  std::string cur_str(cur->def->strval);
  if (cur_str == "all") {
    invars.dangle_type = 2;
  } else if (cur_str == "some") {
    invars.dangle_type = 1;
  } else if (cur_str == "none") {
    invars.dangle_type = 0;
  }
}

void ScriptProcessor::process_domain(int line) { 
  auto & seqs = spec.eval.sequences;
  std::string dom_name(cur->id->strval);
  NUPACK_CHECK(!already_named(dom_name),
      "Name " + dom_name + " already defined by line "+ to_string(line));
  NUPACK_EXC_CHECK(seqs.add_domain_by_constraints(dom_name, cur->def->strval, invars), 
      "Error adding domain " + dom_name);

  names.insert(dom_name);
}

void ScriptProcessor::process_strand(int line) { 
  std::string str_name(cur->id->strval);
  auto & seqs = spec.eval.sequences;
  auto strlist = collect_string_list(cur->def);
  
  NUPACK_CHECK(!already_named(str_name),
      "Name " + str_name + " already defined by line "+ to_string(line));
  NUPACK_EXC_CHECK(seqs.add_strand_by_domains(str_name, strlist),
      "Error adding strand " + str_name);

  names.insert(str_name);
}

void ScriptProcessor::process_structure(int line) {
  std::string name(cur->id->strval);
  NUPACK_CHECK(already_named(name), "Structure being set for undefined complex:"
      + std::string(name));
  
  auto it = complex_map.find(name);
  NUPACK_CHECK(contains_it(it, complex_map), 
      name + " not previously defined as a complex.");
  
  StructureSpec temp(name, cur->def->strval);
  temp.set_strand_names(it->second);
  
  structure_map.emplace(name, std::move(temp));
  
}

void ScriptProcessor::process_complex(int line) {
  std::string name(cur->id->strval);
  NUPACK_CHECK(!already_named(name),
      "Name " + name + " already defined by line" + to_string(line));
  
  auto temp = collect_string_list(cur->def);
  
  complex_map.emplace(name, temp);
  names.insert(name);
}

void ScriptProcessor::process_tube(int line) { 
  auto strlist = collect_string_list(cur->def);
  std::string namestring(cur->id->strval);
  
  NUPACK_CHECK(!already_named(namestring),
      "Name " + namestring + " already defined before line "
      + to_string(line));
  
  auto water_conc = water_density(invars.temperature);
  
  auto temp_tube = TubeSpec();
  temp_tube.set_name(namestring);
  for (auto & str : strlist) 
    temp_tube.add_target(str, DEFAULT_CONCENTRATION / water_conc);
  
  tube_map.emplace(namestring, temp_tube);
  names.insert(namestring);
}

void ScriptProcessor::process_concdef(int line) {
  std::string tube_name(cur->id->strval);
  NUPACK_CHECK(already_named(tube_name),
      "Concentration being set for undefined tube " + tube_name);
  
  auto it = tube_map.find(tube_name);
  NUPACK_CHECK(contains_it(it, tube_map), tube_name + " not previously defined as tube.");
   
  DBL_TYPE curdbl = cur->def->doubleval;
  curdbl = process_conc_units(cur->stype->next, curdbl); 
  NUPACK_CHECK(curdbl >= -1e-50, "Invalid concentration/units specified on "
     "line " + to_string(line));
  
  std::string struc_name(cur->id->next->strval);
  it->second.set_target_conc(struc_name, curdbl);
}

void ScriptProcessor::process_stopdef(int line) { 
  DBL_TYPE curdbl = cur->def->doubleval;
  curdbl = process_frac_units(cur->stype->next, curdbl);
  NUPACK_CHECK(curdbl < 1.0,
      "Invalid stop condition specified, must be < 1.0 or 100 % "
      "on line " + to_string(line));
  NUPACK_CHECK(curdbl >= -1e-50, "Invalid stop condition specified on line " +
      to_string(line) + " stop condition must be > 0.0");
  
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(already_named(namestring),
      namestring + " not found on line " + to_string(line));
  
  auto it = tube_map.find(namestring);
  if (contains_it(it, tube_map)) {
    it->second.set_passive(curdbl * invars.f_passive);
    spec.objectives.add_objective(TubeObjective(namestring, curdbl));
  } else if (contains(namestring, SSM_map)) {
    spec.objectives.add_objective(SymmetryObjective(namestring, curdbl));
  } else if (contains(namestring, structure_map)) {
    spec.objectives.add_objective(StructureObjective(namestring, curdbl));
  } else {
    NUPACK_ERROR("Default objective not defined for " + namestring +
        " on line " + to_string(line));
  }
}

void ScriptProcessor::process_maxsize() { 
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(already_named(namestring), "tube " + namestring + " not found");
  
  auto it = tube_map.find(namestring);
  NUPACK_CHECK(contains_it(it, tube_map), namestring + " not previously defined as tube.")
  
  it->second.set_maxsize(cur->def->intval);
}

void ScriptProcessor::process_prevent(std::vector<Prevent> & prevents) {
  Prevent current;
  current.strands = collect_string_list(cur->id);
  current.prevented_sequences = collect_string_list(cur->def);
  prevents.push_back(std::move(current));
}

void ScriptProcessor::process_library(LibraryMap & libraries) { 
  std::string unitstring("");
  if (cur->stype->next) unitstring.assign(cur->stype->next->strval);
  
  std::string namestring(cur->id->strval);
  auto strlist = collect_string_list(cur->def);
  if (!contains(unitstring, libraries)) {
    libraries[unitstring] = std::map<std::string, std::vector<std::string> >();
  }
  NUPACK_CHECK(!contains(namestring, libraries[unitstring]),
      "Invalid library definition for " + namestring + " already exists");
  
  libraries[unitstring][namestring] = strlist;
}

void ScriptProcessor::process_libseq(LibraryAssignmentMap & library_assignments) { 
  std::string namestring(cur->id->strval);
  std::string unitstring("");
  if (cur->stype->next) unitstring.assign(cur->stype->next->strval);
  
  auto strlist = collect_string_list(cur->def);
  
  auto tmp_lib_val = std::make_pair(unitstring, strlist);
  auto tmp_lib_pair = std::make_pair(namestring, tmp_lib_val);
  library_assignments.insert(tmp_lib_pair);
}

void ScriptProcessor::process_symmetry_min(int line) { 
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(!already_named(namestring),
      "Name " + namestring + " already defined before line "
      + to_string(line));
  
  auto strlist = collect_string_list(cur->def);
  SSM_map.emplace(namestring, SingleSymmetrySpec(namestring, strlist));
  names.insert(namestring);
}

void ScriptProcessor::process_source(int line) { 
  std::string namestring(cur->id->strval);
  
  if (!already_named(namestring)) {
    std::vector<int> nucs = SequenceUtils::str_to_nuc(cur->def->strval);
    source_map.emplace(namestring, SourceSpec(namestring, nucs));
    names.insert(namestring);
  } else {
    auto it = window_map.find(namestring);
    NUPACK_DEBUG("Assigning window source");
    NUPACK_CHECK(contains_it(it, window_map),
        "window source assignment failure on line " + to_string(line));
    
    set_window_source(it->second, cur->def);
  }
}

void ScriptProcessor::process_similarity(int line) { 
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(already_named(namestring),
      "window definition failure on line " + to_string(line));
  
  auto it = window_map.find(namestring);
  NUPACK_CHECK(contains_it(it, window_map),
      "window source assignment failure on line " + to_string(line));
  set_window_similarity(&it->second, cur->stype->next, cur->def);
}

void ScriptProcessor::process_window() { 
  auto & seqs = spec.eval.sequences;
  std::string namestring(cur->id->strval);
  auto strlist = collect_string_list(cur->def);
  add_window(namestring, strlist, seqs);
}

void ScriptProcessor::process_exclude() { 
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(already_named(namestring), "Setting exclude for nonexistent window");
  
  auto it = window_map.find(namestring);
  NUPACK_CHECK(contains_it(it, window_map), "Setting exclude for non-window entity.");
  add_exclusions(&it->second, cur->def);
}

void ScriptProcessor::process_complementary(int & n_comp, int line) { 
  auto strlist = collect_string_list(cur->def);
  auto strlist2 = collect_string_list(cur->id);
  
  std::string complement_type("complement");
  if (cur->stype && cur->stype->next) {
    if (cur->stype->next->strval) {
      complement_type.assign(cur->stype->next->strval);
    } else {
      complement_type.clear();
    }
  }
  NUPACK_CHECK(complement_type == "complement" || complement_type == "wobble",
      "Invalid complementarity type specified on line " + to_string(line));
  
  std::string comp_name = "-comp" + to_string(n_comp);
  complement_map.emplace(comp_name, DomainListSpec(complement_type, strlist2, strlist, line));
  names.insert(comp_name);
  n_comp++;
}

void ScriptProcessor::process_identical(int & n_ident, int line) { 
  auto strlist = collect_string_list(cur->def);
  auto strlist2 = collect_string_list(cur->id);
  
  std::string ident_name = "-ident" + to_string(n_ident);
  identical_map.emplace(ident_name, DomainListSpec("identical", strlist2, strlist, line));
  names.insert(ident_name);
  n_ident++;
}

void ScriptProcessor::process_match(int line) { 
  std::string namestring(cur->id->strval);
  NUPACK_CHECK(!already_named(namestring),
      "Name " + namestring + " already defined by line "
      + to_string(line));
  
  match_spec_map.emplace(namestring, MatchSpec(namestring, std::string(cur->def->strval)));
  names.insert(namestring);
}

void ScriptProcessor::process_off_targets(int line) {  
  auto info = cur->def;
  std::string namestring(cur->id->strval);
  
  NUPACK_CHECK(already_named(namestring), "tube " + namestring + " not found");
  auto it = tube_map.find(namestring);
  NUPACK_CHECK(contains_it(it, tube_map), namestring + " not previously defined as tube.")
  auto & cur_tube = it->second;
  
  // set maxsize
  NUPACK_CHECK(NP_TT_INTEGER == info->type || NP_TT_FLOAT == info->type, 
      "max size must be integer valued")
  cur_tube.set_maxsize(info->intval);
  
  
  // collect whitelist
  if (info->next != NULL) {
    info = info->next;
    auto comps = collect_complexes(info->list);
    make_off_target_list(comps, cur_tube.get_whitelist(), line);
    
    // collect blacklist    
    if (info->next != NULL) {
      info = info->next;
      comps = collect_complexes(info->list);
      make_off_target_list(comps, cur_tube.get_blacklist(), line);
    }
  }
}

void ScriptProcessor::make_off_target_list(NameList & names, NameList & list, int line) {
  const auto & strands = spec.eval.sequences.get_strand_map();
  for (auto & name : names) {
    for (auto & n : name) 
      NUPACK_CHECK(already_named(n), "off-target component \"" + n 
          + "\" on line " + to_string(line) + " not previously defined");
    
    if (name.size() == 1 && contains(name[0], structure_map)) {
      auto it = structure_map.find(name[0]);
      list.push_back(it->second.get_strand_names());
    }
    else {
      for (auto & n : name) 
        NUPACK_CHECK(contains(n, strands), "off-target component \"" + n 
            + "\" on line " + to_string(line) + " is not a previously defined strand");
      
      list.push_back(name);
    }
  }
}

void ScriptProcessor::process_design() {
  NUPACK_DEBUG("start processing");
  auto & seqs = spec.eval.sequences;

  int line = 1;
  cur = root;
  
  std::vector<Prevent> prevents;

  LibraryMap libraries;
  LibraryAssignmentMap library_assignments;
  LibraryValue tmp_lib_val;
  LibraryPair tmp_lib_pair;

  WeightSpec weights;
  DBL_TYPE global_stop = -1;

  int n_comp = 0, n_ident = 0;

  while (cur) {
    if (cur->type == NP_TT_DEFINITION) {
      switch(cur->stype->intval) {
        /* Physical parameters */
        case TOK_TEMPERATURE:
          process_temperature();
          break;
        case TOK_MATERIAL:
          process_material();
          break;
        case TOK_SODIUM:
          process_sodium();
          break;
        case TOK_MAGNESIUM:
          process_magnesium();
          break;
        case TOK_GLOBAL_STOPDEF:
          global_stop = process_global_stopdef();
          break;
        case TOK_DANGLES:
          process_dangles();
          break;
        
        /* Optimization parameters */
        case TOK_SEED:
          invars.seed = (unsigned int) cur->def->intval;
          break;
        case TOK_HSPLIT:
          invars.H_split = cur->def->intval;
          break;
        case TOK_NSPLIT:
          invars.N_split = cur->def->intval;
          break;
        case TOK_MBAD:
          invars.M_bad = (DBL_TYPE) cur->def->intval;
          break;
        case TOK_MREOPT:
          invars.M_reopt = (DBL_TYPE) cur->def->intval;
          break;
        case TOK_MRESEED:
          invars.M_reseed = (DBL_TYPE) cur->def->intval;
          break;
        case TOK_FSPLIT:
          invars.f_split = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FSTRINGENT:
          invars.f_stringent = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FREDECOMP:
          invars.f_redecomp = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FREFOCUS:
          invars.f_refocus = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FPASSIVE:
          invars.f_passive = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_GC_INIT:
          invars.gc_init_prob = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_DGCLAMP:
          invars.dG_clamp = (DBL_TYPE) cur->def->doubleval;
          break;
        
        /* Optimization modes */
        case TOK_ALLOWWOBBLE:
          invars.allow_wobble = cur->def->intval;
          break;
        case TOK_ALLOWMISMATCH:
          invars.allow_mismatch = cur->def->intval;
          break;
        case TOK_DISABLEMUTWEIGHTS:
          invars.disable_defect_weights = cur->def->intval;
          break;
        case TOK_REDECOMPOSE:
          invars.redecompose = cur->def->intval;
          break;
        case TOK_POPULATION:
          invars.N_population = cur->def->intval;
          break;
        case TOK_TRIALS:
          invars.N_trials = cur->def->intval;
          break;
        case TOK_OPTTIME:
          invars.allowed_opt_time = (DBL_TYPE) cur->def->doubleval;
          break;
        
        /* Output parameters */
        case TOK_PRINTLEAVES:
          invars.print_leaves = cur->def->intval;
          break;
        case TOK_PRINTSTEPS:
          invars.print_steps = cur->def->intval;
          break;
        case TOK_MINPAIR:
          invars.min_ppair = (DBL_TYPE) cur->def->doubleval;
          break;

        /* Physical components */
        case TOK_DOMAIN:
          process_domain(line);
          break;
        case TOK_STRAND:
          process_strand(line);
          break;
        case TOK_STRUCTURE:
          process_structure(line);
          break;
        case TOK_WEIGHT:
          set_item_weight(cur, line, weights);
          break;
        case TOK_TUBE:
          process_tube(line);
          break;
        
        case TOK_CONCDEF: // Segfaults if trying to set concentration for tube without specific structure.
          process_concdef(line);
          break;
        case TOK_STOPDEF:
          process_stopdef(line);
          break;
        case TOK_MAXSIZE:
          process_maxsize();
          break;
        case TOK_OFFTARGETS:
          process_off_targets(line);
          break;
        
        /* constraints */
        case TOK_PREVENT:
          process_prevent(prevents);
          break;
        case TOK_LIBRARY:
          process_library(libraries);
          break;
        case TOK_LIBSEQ:
          process_libseq(library_assignments);
          break;
        case TOK_SYMMETRY_MIN:
          process_symmetry_min(line);
          break;
        case TOK_WORDSIZE:
          set_ssm_wordsize(cur, line);
          break;
        
        /* External sequence constraints */
        case TOK_SOURCE:
          process_source(line);
          break;
        case TOK_SIMILARITY:
          process_similarity(line);
          break;
        case TOK_WINDOW:
          process_window();
          break;
        case TOK_ACCESSIBLE:
          NUPACK_LOG_WARN("Accessibility constraints are not enabled; they will have no effect.");
          break;
        case TOK_EXCLUDE:
          process_exclude();
          break;
        case TOK_COMPLEMENTARY:
          process_complementary(n_comp, line);
          break;
        case TOK_IDENTICAL:
          process_identical(n_ident, line);
          break;
        case TOK_MATCH:
          process_match(line);
          break;
        case TOK_MATCHRANGE:
          add_match_domain(cur, line);
          break;
        case TOK_COMPLEX:
          process_complex(line);
          break;
        
        default:
          NUPACK_ERROR("Invalid definition code " + to_string(cur->stype->intval) +
              " on line " + to_string(line));
      }
    }
    line += 1;
    cur = cur->next;
  }
  invars.deduce_h_split();

  // Most absolute properties are now filled out. Names are currently in
  // place instead of reference indices. This fills out reference indices
  resolve_all(weights);
  spec.objectives.set_weights(weights);
  add_implied_constraints();
  add_match_constraints();
  add_pattern_constraints(prevents);
  add_library_constraints(libraries, library_assignments);
  
  // post-process window constraint
  for (auto & item : window_map) {
    const auto & nuc_to_vars = spec.eval.sequences.get_variables();
    
    auto & cur_wind = item.second;
    std::vector<std::vector<int> > cur_defs = cur_wind.get_constraints();
    std::vector<int> nuc_ids = cur_wind.get_nuc_ids();
    
    int n_poss = cur_defs.size();
    std::vector<trinary> poss(n_poss, true);
    int supp_var = spec.constraints.add_variable(poss);
    
    std::vector<int> cur_vars;
    for (auto j_nuc = 0; j_nuc < nuc_ids.size(); j_nuc++) {
      cur_vars.push_back(nuc_to_vars[nuc_ids[j_nuc]]);
    }
    
    std::vector<std::string> cur_def_strings;
    for (const auto & cur_def : cur_defs) {
      std::string tmp = SequenceUtils::nuc_to_str(cur_def, invars.material);
      cur_def_strings.push_back(tmp);
    }
    
    spec.constraints.add_constraint(WordConstraint(cur_vars, cur_def_strings, supp_var));
  }
  
  for (auto & item : complement_map) {
    auto & cur_dlist = item.second;
    std::vector<int> nuc_ids1 = cur_dlist.get_nuc_ids1(seqs);
    std::vector<int> nuc_ids2 = cur_dlist.get_nuc_ids2(seqs);
    
    NUPACK_CHECK(nuc_ids1.size() == nuc_ids2.size(),
        "Length mismatch in constraint on line "
        + to_string(cur_dlist.get_line()));
    
      auto comp_type = NUPACK_CS_STRONG;
      for (auto j_nuc = 0; j_nuc < nuc_ids1.size(); j_nuc++) {
        spec.constraints.add_constraint(
            CompConstraint(nuc_ids1[j_nuc], 
            nuc_ids2[nuc_ids2.size() - j_nuc - 1], comp_type));
      }
  }
  
  for (auto & item : wobble_map) {
    auto & cur_dlist = item.second;
    std::vector<int> nuc_ids1 = cur_dlist.get_nuc_ids1(seqs);
    std::vector<int> nuc_ids2 = cur_dlist.get_nuc_ids2(seqs);
    
    NUPACK_CHECK(nuc_ids1.size() == nuc_ids2.size(),
        "Length mismatch in constraint on line "
        + to_string(cur_dlist.get_line()));
    
      auto comp_type = NUPACK_CS_WEAK;
      for (auto j_nuc = 0; j_nuc < nuc_ids1.size(); j_nuc++) {
        spec.constraints.add_constraint(
            CompConstraint(nuc_ids1[j_nuc], 
            nuc_ids2[nuc_ids2.size() - j_nuc - 1], comp_type));
      }
  }
  
  for (auto & item : identical_map) {
    auto & cur_dlist = item.second;
    std::vector<int> nuc_ids1 = cur_dlist.get_nuc_ids1(seqs);
    std::vector<int> nuc_ids2 = cur_dlist.get_nuc_ids2(seqs);
    
    NUPACK_CHECK(nuc_ids1.size() == nuc_ids2.size(),
        "Length mismatch in constraint on line "
        + to_string(cur_dlist.get_line()));
    
    for (auto j_nuc = 0; j_nuc < nuc_ids2.size(); j_nuc++) {
      spec.constraints.add_constraint(
          IdentConstraint(nuc_ids1[j_nuc], nuc_ids2[j_nuc]));
    }
  }
  
  if (invars.add_global_stop || global_stop > -1) {
    if (global_stop < 0.0) global_stop = NUPACK_DEF_STOP; 
    add_global_objective(global_stop);
  }

  if (invars.add_default_stops) add_default_stop_conditions(); 

  NUPACK_CHECK(spec.objectives.get_names().size() > 0, 
      "No objectives provided! No design can be performed!");
}

int ScriptProcessor::parse_design() {
  std::ifstream npfile(filename, std::fstream::in);
  std::stringstream buffer;
  NUPACK_DEBUG("Before parsing");
  NUPACK_CHECK(npfile.is_open(), "Error reading file");

  if (npfile.good()) buffer << npfile.rdbuf();

  root = get_ast(buffer.str().c_str());
  NUPACK_CHECK(root != NULL, "Error parsing file");
  process_design();

  if (invars.seed == 0) {
    timeval curtime;
    gettimeofday(&curtime, NULL);
    invars.seed = (unsigned int) curtime.tv_usec;
  }

  return ERR_OK;
}

std::pair<std::string, std::string> ScriptProcessor::split_filename(std::string filename) {
  std::string prefix = filename;
  int filename_len = filename.size();
  if (filename_len > 3) {
    std::string suffix(filename.begin() + filename_len - 3, filename.end());
    if (suffix == ".np") {
      prefix.assign(filename.begin(), filename.begin() + (filename_len - 3));
    } else {
      filename += ".np";
    }
  } else {
    filename += ".np";
  }
  return {filename, prefix};
}
}
