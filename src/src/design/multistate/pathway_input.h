#pragma once

#include "parsestruc.h"
#include "constraint_handler.h"
#include "physical_spec.h"
#include "weight_spec.h"
#include "design_spec.h"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <tuple>

namespace nupack {
  using Complex = std::vector<std::string>;
  
  struct Prevent {
    std::vector<std::string> strands;
    std::vector<std::string> prevented_sequences;
    
    bool has_specific_targets() const { return !strands.empty(); }
  };
  
  class ScriptProcessor {
    using Names = std::set<std::string>;
    
    using LibraryValue = std::pair<std::string, std::vector<std::string> >;
    using LibraryPair = std::pair<std::string, LibraryValue>;
    using LibraryMap = std::map<std::string, std::map<std::string, std::vector<std::string> > >;
    using LibraryAssignmentMap = std::multimap<std::string, LibraryValue>;
    
    using PreventTargets = std::vector<std::string>;
    using Prevents = std::vector<std::vector<std::string>>;

    public:
      ScriptProcessor(std::string & script_name, DesignSpec & spec) : 
          spec(spec), invars(spec.eval.options) {
        std::tie(filename, invars.file_prefix) = split_filename(script_name);
      }
      ~ScriptProcessor() { np_tt_destroy_value_struc(root); }

      bool already_named(std::string & name) { return !(names.find(name) == names.end()); }

      
      void add_library_constraints(LibraryMap lib_defs, LibraryAssignmentMap lib_assignments);
      void add_match_constraints();
      void add_pattern_constraints(std::vector<Prevent> &);
      void add_implied_constraints();
      
      void set_window_source(WindowSpec & wind, value_struc_t * source);
      void set_window_similarity(WindowSpec * wind, value_struc_t * unitsource, value_struc_t * range);
      void add_window(std::string name, std::vector<std::string> domains, const SequenceSpec & seqs);
      
      void add_exclusions(WindowSpec * wind, value_struc_t * vals);
      void set_ssm_wordsize(value_struc_t * cur, int line);
      void set_item_weight(value_struc_t * cur, int line, WeightSpec & weights);
      void add_match_domain(value_struc_t * cur, int line);
      
      void gather_independent_entities(Names & tubes_to_add, Names & strucs_to_add);
      
      void add_global_objective(DBL_TYPE stop);
      void add_default_stop_conditions();
      
      DBL_TYPE process_conc_units(value_struc_t * cur, DBL_TYPE baseval);
      DBL_TYPE process_frac_units(value_struc_t * cur, DBL_TYPE baseval);
      DBL_TYPE get_frac_units(const std::string & unit_name);

      int parse_design();
      void process_design();
    
      static void convert_to_lower(std::string & str) { for (auto & s : str) s = tolower(s); }
      static std::vector<std::string> collect_string_list(value_struc_t * linkedlist);
      static NameList collect_complexes(value_struc_t * list);
      static std::pair<std::string, std::string> split_filename(std::string filename);
    private:
      DesignSpec & spec;
      NupackInvariants & invars;
      Names names;
      std::string filename;
      
      std::map<std::string, SingleSymmetrySpec> SSM_map;
      std::map<std::string, StructureSpec> structure_map;
      std::map<std::string, SingleParamSpec> specset_map;
      std::map<std::string, TubeSpec> tube_map;
      std::map<std::string, MatchSpec> match_spec_map;
      std::map<std::string, SourceSpec> source_map;
      std::map<std::string, WindowSpec> window_map;
      std::map<std::string, DomainListSpec> complement_map;
      std::map<std::string, DomainListSpec> identical_map;
      std::map<std::string, DomainListSpec> wobble_map;
      std::map<std::string, Complex> complex_map;

      value_struc_t * root;
      value_struc_t * cur;
      
      void resolve_all(WeightSpec & weights);
      void resolve_SSM();
      void resolve_weights(WeightSpec & weights);
      
      void generate_default_parameter_spec();
      
      // used during AST traversal for the different token types
      void process_temperature();
      void process_material();
      void process_sodium();
      void process_magnesium();
      DBL_TYPE process_global_stopdef();
      void process_dangles();
      void process_domain(int line);
      void process_strand(int line);
      void process_structure(int line);
      void process_tube(int line);
      void process_concdef(int line);
      void process_stopdef(int line);
      void process_maxsize();
      void process_prevent(std::vector<Prevent> &);
      void process_library(LibraryMap & libraries);
      void process_libseq(LibraryAssignmentMap & library_assignments);
      void process_symmetry_min(int line);
      void process_source(int line);
      void process_similarity(int line);
      void process_window();
      void process_exclude();
      void process_complementary(int & n_comp, int line);
      void process_identical(int & n_ident, int line);
      void process_match(int line);
      void process_complex(int line);
      
      void process_off_targets(int line);
      void make_off_target_list(NameList & names, NameList & list, int line);
  };    
}
