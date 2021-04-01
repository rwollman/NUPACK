/** @brief contains classes for holding all the conceptual objects in a design
  (Design) and an object representing the design logic*/
#pragma once
#include "Design.h"
#include "DesignParameters.h"
#include "DesignComponents.h"
#include "Objectives.h"
#include "Result.h"
#include "Weights.h"

namespace nupack {
namespace newdesign {

struct Designer {
    Design design;
    vec<Objective> objectives;

    DesignParameters parameters;
    Weights weights;
    EnsemblePartition Psi;

    int max_depth;

    DesignStats stats;
    Timer timer;
    Logs logs;
    EngineObserver obs;

    ResultState best {inf_result};
    // ArchiveState archive;

    std::set<Sequence> known_bads;

    std::function<void(Designer &, bool)> checkpoint {NoOp()}; // checks if should checkpoint and does

    Designer() = default;
    Designer(Design d, vec<Objective> objs, Weights weights, DesignParameters params={}) :
            design(std::move(d)), objectives(std::move(objs)),
            parameters(std::move(params)),
            weights(std::move(weights)),
            Psi(design.complexes, parameters.f_passive * parameters.f_stop),
            logs(parameters.log_file_paths()),
            obs{parameters.slowdown, &logs} {}

    NUPACK_REFLECT(Designer, design, objectives, parameters, Psi, max_depth, weights, stats, best, known_bads)

    void initialize(bool decompose=true);
    void subset_decompose(vec<uint> subset, uint depth=0);

    void redecompose_active(Local const &env, uint depth);
    bool redecompose(uint depth, Sequence const &sequence);
    void refocus(Local const &env, Sequence const &sequence);


    /* alternate stuff */
    // Result alternate_optimize_tubes(Local const &env);
    // bool length_extrapolation_refocus(Local const &env);
    // bool sum_pf_refocus(Local const &env);
    /* end stuff */


    Result optimize_tubes(Local const &env);
    Result optimize_tubes_impl(Local const &env);
    Result optimize_forest(Local const &env, Sequence seq);
    Result optimize_leaves(Local const &env, Sequence seq);
    Result mutate_leaves(Local const &env, Sequence seq);

    /* multiobjective */
    Result evaluate_objectives(Local const &env, uint depth, EnsemblePartition const &part, Weights const &weights);
    Result reevaluate_objectives(Local const &env, Result const &res, uint depth, EnsemblePartition const &part, Weights const &weights);

    Sequence best_sequence(Local const &env);

    bool improvement_slowing(vec<uint> const &x, vec<real> const &y);

    auto time_elapsed() const { return stats.design_time + timer.elapsed(); }

    bool success() const {return best.full.weighted_total() <= parameters.f_stop;}
    void time_analysis(Local const &env);

    static constexpr auto repr_names() {return make_names("design", "parameters", "weights", "Psi", "stats", "timer", "best", "max_depth", "known_bads");}

    auto save_repr() const {return make_members(design, parameters, weights, Psi, stats,
            timer, best, max_depth, known_bads);}

    void load_repr(Design design, DesignParameters parameters, Weights weights, EnsemblePartition Psi, DesignStats stats,
            Timer timer, ResultState best, int max_depth, std::set<Sequence> known_bads) {
        *this = Designer();
        this->design = std::move(design);
        this->parameters = std::move(parameters);
        this->weights = std::move(weights);
        this->Psi = std::move(Psi);
        this->stats = std::move(stats);
        this->timer = std::move(timer);
        this->logs = Logs(parameters.log_file_paths());
        this->obs = {parameters.slowdown, &logs};
        this->best = std::move(best);
        this->max_depth = std::move(max_depth);
        this->known_bads = std::move(known_bads);
    }
};




}}
