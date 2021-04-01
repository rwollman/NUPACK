#pragma once

#include "adapter.h"
// #include <json/json.h>

#include "../../model/ModelVariants.h"

namespace nupack {
namespace custom_csp {

enum class print_level {
    NONE = 0,
    REFOCUS = 1,
    REDECOMPOSE = 2,
    REOPTIMIZE = 3,
    RESEED = 4,
    ALL = 5
};

struct NupackInvariants {
    real temperature {NUPACK_DEF_TEMPERATURE};                          // temperature (K)
    real sodium {NUPACK_DEF_SODIUM};                                    // sodium concentration (M)
    real magnesium {NUPACK_DEF_MAGNESIUM};                              // magnesium concentration (M)
    real f_sparse {NUPACK_DEF_MIN_PPAIR};                              // minimum pair probability saved
    int M_bad {NUPACK_DEF_M_BAD};                                       // scaled number of unfavorable leaf mutations to allow
    int M_reopt {NUPACK_DEF_M_REOPT};                                   // number of failed leaf reoptimizations to allow
    int M_reseed {NUPACK_DEF_M_RESEED};                                 // number of nucleotides to reseed at beginning of leaf reoptimization
    real f_split {NUPACK_DEF_F_SPLIT};                                  // minimum pair prob of helix for ppair decomp
    real f_passive {NUPACK_DEF_F_PASSIVE};                              // fraction that concentrations are deflated by
    real f_stringent {NUPACK_DEF_F_STRINGENT};                          // Margin for relaxation during tree optimization
    real f_redecomp {NUPACK_DEF_F_REDECOMP};                            // Margin of correction during decomposition correction
    real f_refocus {NUPACK_DEF_F_REFOCUS};                              // Margin of correction during decomposition correction
    real gc_init_prob {NUPACK_DEF_GC_INIT_PROB};                        // Probability of choosing GC vs AU on initialization (UNUSED)
    real dG_clamp { -20.0};                                             // clamp bonus for enforcing pairs (different from tubedesign)
    parameter_set material {NUPACK_DEF_MATERIAL};                       // DNA1998/RNA1995/RNA1999 type
    Ensemble ensemble = Ensemble::min;                                  // method of accounting for dangle energies
    unsigned int seed {NUPACK_DEF_SEED};                                // RNG seed
    int H_split { -1};                                                  // minimum number of base pairs on either side of a split point, -1 == unset sentinel value
    int N_split {NUPACK_DEF_N_SPLIT};                                   // minimum number of bases in a child ensemble
    int N_trials {1};                                                   // number of separate seeds to run the design with (NOT IMPLEMENTED)
    bool print_leaves {false};                                          // UNUSED
    print_level print_steps {print_level::NONE};                                       // print intermediate design evaluation
    bool allow_wobble {NUPACK_DEF_ALLOW_WOBBLE};
    bool allow_mismatch {NUPACK_DEF_ALLOW_MISMATCH};
    bool use_long_helix {NUPACK_DEF_USE_LONG_HELIX};
    bool disable_defect_weights {NUPACK_DEF_DISABLE_DEFECT_WEIGHTS};    // UNUSED
    bool disable_focus {NUPACK_DEF_DISABLE_FOCUS};                      // UNUSED
    bool forbid_splits {NUPACK_DEF_FORBID_SPLITS};                      // UNUSED
    bool redecompose {NUPACK_DEF_REDECOMPOSE};                          // UNUSED
    bool include_dummies {false};                                       // UNUSED
    bool add_default_stops {false};
    bool add_global_stop {false};
    bool print_json {false};
    bool print_ppairs {false};
    string file_prefix {""};
    real elapsed_time {0};
    real allowed_opt_time {86000000};

    string material_string;
    string start_timestamp;
    real start_time;

    NupackInvariants();

    void deduce_h_split();

    string mat_str() const;
    string dangle_str() const;
    void serialize(std::ostream & out, int indent = 0, string prefix = "") const;

    bool opt_time_elapsed() const;
#ifdef JSONCPP_FOUND
    Json::Value make_json_value() const;
#endif

    static int num_mutations;             // to get a speed-invariant metric to compare with newdesign
    static int num_redecompositions;
};

}
}
