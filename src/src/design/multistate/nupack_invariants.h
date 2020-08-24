#pragma once

#include <thermo.h>
#include <json/json.h>

#include <string>

namespace nupack {
  enum print_steps_e {
    PRINT_NONE = 0,
    PRINT_REFOCUS = 1,
    PRINT_REDECOMPOSE = 2,
    PRINT_REOPTIMIZE = 3,
    PRINT_RESEED = 4,
    PRINT_ALL = 5
  };
      
  class NupackInvariants {
    public:
      DBL_TYPE temperature {NUPACK_DEF_TEMPERATURE};                    // temperature (K)
      DBL_TYPE sodium {NUPACK_DEF_SODIUM};                              // sodium concentration (M)
      DBL_TYPE magnesium {NUPACK_DEF_MAGNESIUM};                        // magnesium concentration (M)
      DBL_TYPE min_ppair {NUPACK_DEF_MIN_PPAIR};                        // minimum pair probability saved
      int M_bad {NUPACK_DEF_M_BAD};                             // scaled number of unfavorable leaf mutations to allow
      int M_reopt {NUPACK_DEF_M_REOPT};                               // number of failed leaf reoptimizations to allow
      int M_reseed {NUPACK_DEF_M_RESEED};                               // number of nucleotides to reseed at beginning of leaf reoptimization
      DBL_TYPE f_split {NUPACK_DEF_F_SPLIT};                            // minimum pair prob of helix for ppair decomp
      DBL_TYPE f_passive {NUPACK_DEF_F_PASSIVE};                        // fraction that concentrations are deflated by
      DBL_TYPE f_stringent {NUPACK_DEF_F_STRINGENT};                    // Margin for relaxation during tree optimization
      DBL_TYPE f_redecomp {NUPACK_DEF_F_REDECOMP};                      // Margin of correction during decomposition correction
      DBL_TYPE f_refocus {NUPACK_DEF_F_REFOCUS};                        // Margin of correction during decomposition correction
      DBL_TYPE gc_init_prob {NUPACK_DEF_GC_INIT_PROB};                  // Probability of choosing GC vs AU on initialization (UNUSED)
      DBL_TYPE dG_clamp {-20.0};                                        // clamp bonus for enforcing pairs (different from tubedesign)
      parameter_set material {NUPACK_DEF_MATERIAL};                     // DNA1998/RNA1995/RNA1999 type
      int dangle_type {NUPACK_DEF_DANGLE_TYPE};                         // method of accounting for dangle energies 
      unsigned int seed {NUPACK_DEF_SEED};                              // RNG seed
      int H_split {-1};                                                 // minimum number of base pairs on either side of a split point, -1 == unset sentinel value
      int N_split {NUPACK_DEF_N_SPLIT};                                 // minimum number of bases in a child ensemble
      int N_trials {1};                                                 // number of separate seeds to run the design with (NOT IMPLEMENTED)
      int N_population {NUPACK_DEF_POPULATION};                         // populations used for Pareto dominance work
      bool print_leaves {false};                                        // UNUSED
      int print_steps {PRINT_NONE};                                     // print intermediate design evaluation
      bool allow_wobble {NUPACK_DEF_ALLOW_WOBBLE};      
      bool allow_mismatch {NUPACK_DEF_ALLOW_MISMATCH};
      bool use_long_helix {NUPACK_DEF_USE_LONG_HELIX};    
      bool disable_defect_weights {NUPACK_DEF_DISABLE_DEFECT_WEIGHTS};  // UNUSED
      bool disable_focus {NUPACK_DEF_DISABLE_FOCUS};                    // UNUSED
      bool forbid_splits {NUPACK_DEF_FORBID_SPLITS};                    // UNUSED
      bool redecompose {NUPACK_DEF_REDECOMPOSE};                        // UNUSED
      bool include_dummies {false};                                     // UNUSED
      bool add_default_stops {false};
      bool add_global_stop {false};
      bool print_json {false};
      bool print_ppairs {false};
      std::string file_prefix {""};
      DBL_TYPE elapsed_time {0};
      DBL_TYPE allowed_opt_time {86000000};

      std::string material_string;
      std::string start_timestamp;
      DBL_TYPE start_time;

      NupackInvariants();
      
      void deduce_h_split();
      
      std::string mat_str() const;
      std::string dangle_str() const;
      void serialize(std::ostream & out, int indent = 0, std::string prefix = "") const;

      bool opt_time_elapsed() const;
#ifdef JSONCPP_FOUND
      Json::Value make_json_value() const;
#endif
  };
}
