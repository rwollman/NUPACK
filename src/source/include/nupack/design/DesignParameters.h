#pragma once
#include "../common/Config.h"
#include "../common/Random.h"

namespace nupack { namespace newdesign {

struct DesignParameters {
    /** seed for the StaticRNG; 0 is sentinel value indicating that Static_RD should be used */
    uint rng_seed {0};
    /** stop condition (attempt to find sequence with defect lower than this) */
    real f_stop {0.02};
    /** the fraction of the stop condition to allocate to off-targets in passive
     * set when using ensemble focusing */
    real f_passive {0.01};
    /** the number of flanking base pairs required to be in a structure on
     * either side of a split point */
    int H_split {2};
    /** the minimum number of nucleotides that must be in a leaf node in the
     * decomposition tree */
    int N_split {12};
    /** the fraction of the partition function that must be captured by a set of
     * exclusive split points during decomposition */
    real f_split {0.99};
    /** for the estimate of the defect at depth d, d factors of f_stringent are
     * applied to f_stop to yield the stop condition for the estimate at depth d
     * */
    real f_stringent {0.99};
    /** the bonus free energy applied to each enforced pair in a decomposed
     * structure (to ensure that it forms with probability near 1) */
    real dG_clamp {-20};
    /** the maximum number of times non-improving sequences can be reencountered
     * during leaf mutation without finding a single improving sequence before
     * leaf mutation exits in failure */
    int M_bad {300};
    /** the number of nucleotide variables that are mutated sequentially without
     * intermediate evaluation during leaf reseeding */
    int M_reseed {50};
    /** the maximum number of times leaf reseeding and subsequent leaf
     * reoptimization can occur without finding a sequence with a better defect
     * before leaf optimization exits in failure */
    int M_reopt {3};
    /** Not currently used because all active complexes are redecomposed when
     * redecomposition is called for. When redecomposing individual complexes,
     * this is the fraction of the initial difference between estimates at depth
     * d and d+1 (after dividing out f_stringent) that is allowed to remain
     * after redecomposing. */
    real f_redecomp {0.03};
    /** When refocusing, this is the fraction of the initial difference between
     * the full root ensemble defect and the focused ensemble defect estimate
     * that is allowed to remain after refocusing. */
    real f_refocus {0.03};
    /** number of bytes of RAM split evenly amongst the models in the design for their caches
     */
    std::size_t cache_bytes_of_RAM {0};
    /**
     * The maximum archive size at the leaves and hence an upper bound on the number of returned solutions
     */
    // uint archive_size {10};

    /**
     * whether to run final analysis calculation at the end (default to false for release)
     */
    bool time_analysis {false};

    /**
     * log file paths
     */
    string log {};
    string decomposition_log {};
    string thermo_log {};

    std::map<string, string> log_file_paths() const {
        std::map<string, string> paths;
        if (!log.empty()) paths.emplace("basic", log);
        if (!decomposition_log.empty()) paths.emplace("decomposition", decomposition_log);
        if (!thermo_log.empty()) paths.emplace("thermo", thermo_log);
        return paths;
    }

    /**
     * The cutoff for ppairs that make it from the dense pair probability matrix
     * into the sparse matrices used during design.
     */
    real f_sparse {0.00001};

    /**
     * For profiling; allows running all thermodynamics code a multiple of times
     * to disentangle design from thermo time contributions
     */
    uint slowdown {0};

    // Set the global seed based on the one held here. Do not change the one held here.
    void init_rng() const;

    NUPACK_REFLECT(DesignParameters, rng_seed, f_stop, f_passive, H_split, N_split, f_split,
            f_stringent, dG_clamp, M_bad, M_reseed, M_reopt, f_redecomp, f_refocus,
            cache_bytes_of_RAM, f_sparse, slowdown, log, decomposition_log, thermo_log,
            time_analysis);
};


}}
