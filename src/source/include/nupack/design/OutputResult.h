#pragma once
#include "Design.h"
#include "DesignComponents.h"
#include "DesignParameters.h"
#include "Objectives.h"
#include "Weights.h"

namespace nupack { namespace newdesign {

struct Designer;
struct Result;

struct ComplexResult {
    string name;
    ::nupack::Complex sequence;
    Structure structure;
    ProbabilityMatrix pair_probabilities;
    real log_partition_function;
    real defect;
    real normalized_defect;

    ComplexResult() = default;

    NUPACK_REFLECT(ComplexResult, name, sequence, structure, pair_probabilities, log_partition_function, defect, normalized_defect);
};


struct TubeComplex {
    string name;
    real concentration;
    real target_concentration;
    real defect;
    real structural_defect;
    real concentration_defect;
    real normalized_defect_contribution;

    TubeComplex() = default;

    NUPACK_REFLECT(TubeComplex, name, concentration, target_concentration, defect, structural_defect, concentration_defect, normalized_defect_contribution);
};

struct TubeResult {
    string name;
    real nucleotide_concentration;
    real defect;
    real normalized_defect;
    vec<TubeComplex> complexes;

    TubeResult() = default;

    NUPACK_REFLECT(TubeResult, name, nucleotide_concentration, defect, normalized_defect, complexes);
};

struct SingleResult {
    std::map<string, Sequence> domains;
    std::map<string, Sequence> strands;
    vec<ComplexResult> complexes;
    vec<TubeResult> tubes;

    vec<real> defects;
    vec<real> weighted_defects;

    SingleResult() = default;
    SingleResult(Designer const &, Result const &);

    NUPACK_REFLECT(SingleResult, domains, strands, complexes, tubes, defects, weighted_defects);
};

struct DesignResult {
    Model<real> model;
    DesignParameters parameters;
    DesignStats stats;
    vec<Objective> objectives;

    vec<SingleResult> results;
    Weights weights;

    bool success;

    DesignResult() = default;
    // DesignResult(Specification const &, Designer const &, Result const &);
    // DesignResult(Specification const &spec, Designer const &designer) : DesignResult(spec, designer, designer.best.full) {}
    DesignResult(Designer const &);

    NUPACK_REFLECT(DesignResult, model, parameters, stats, objectives, results, weights, success);
};


}}
