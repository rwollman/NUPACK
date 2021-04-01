#pragma once
#include "Defect.h"
#include "Granularity.h"
#include "Logging.h"
#include <chrono>

namespace nupack { namespace newdesign {

vec<real> ord_lin_lsq(vec<real> const &x, vec<real> const &y);

struct Timer {
    using simple_type = True;
    using clock = std::chrono::high_resolution_clock;
    using time = decltype(std::chrono::high_resolution_clock::now());
    using duration = std::chrono::duration<real>;

    time _start;
    time _stop;

    Timer & start();
    real elapsed() const;
    real stop();

    auto save_repr() const {
        return std::pair<real, real>(
                duration(_start.time_since_epoch()).count(),
                duration(_stop.time_since_epoch()).count());
    }

    void load_repr(std::pair<real, real> x) {
        _start = time(std::chrono::duration_cast<clock::duration>(std::chrono::duration<real>(x.first)));
        _stop = time(std::chrono::duration_cast<clock::duration>(std::chrono::duration<real>(x.second)));
    }
};



struct DesignStats {
    uint num_leaf_evaluations {0};
    uint num_reseeds {0};
    vec<uint> num_redecompositions;
    vec<uint> offtargets_added_per_refocus;
    real design_time {0};
    real analysis_time {0};
    EnsemblePartition final_Psi;

    NUPACK_REFLECT(DesignStats, num_leaf_evaluations, num_reseeds, num_redecompositions, offtargets_added_per_refocus, design_time, analysis_time, final_Psi);
};



}}
