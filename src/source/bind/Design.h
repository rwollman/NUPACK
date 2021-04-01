#pragma once
#include "Bind.h"
#include <nupack/design/Objectives.h>
#include <nupack/design/OutputResult.h>
#include <nupack/design/Specification.h>
#include <nupack/design/Design.h>
#include <nupack/design/DesignComponents.h>

namespace nupack::newdesign {

/******************************************************************************/

void render(Document &doc, Type<ComplexObjective> t);
void render(Document &doc, Type<ComplexResult> t);
void render(Document &doc, Type<ComplexSpec> t);
void render(Document &doc, Type<ConstraintSpec> t);
void render(Document &doc, Type<DesignParameters> t);
void render(Document &doc, Type<DesignResult> t);
void render(Document &doc, Type<DesignStats> t);
void render(Document &doc, Type<DiversitySpec> t);
void render(Document &doc, Type<DomainSpec> t);
void render(Document &doc, Type<DualListSpec> t);
void render(Document &doc, Type<EnergyEqualizationObjective> t);
void render(Document &doc, Type<EnsemblePartition> t);
void render(Document &doc, Type<MultitubeObjective> t);
void render(Document &doc, Type<Objective> t);
void render(Document &doc, Type<PatternObjective> t);
void render(Document &doc, Type<PatternSpec> t);
void render(Document &doc, Type<ReversedComplex> t);
void render(Document &doc, Type<SSMObjective> t);
void render(Document &doc, Type<Sequence> t);
void render(Document &doc, Type<SimilarityObjective> t);
void render(Document &doc, Type<SimilaritySpec> t);
void render(Document &doc, Type<SingleResult> t);
void render(Document &doc, Type<Specification> t);
void render(Document &doc, Type<Strand> t);
void render(Document &doc, Type<StrandSpec> t);
void render(Document &doc, Type<Structure> t);
void render(Document &doc, Type<Timer> t);
void render(Document &doc, Type<TubeComplex> t);
void render(Document &doc, Type<TubeObjective> t);
void render(Document &doc, Type<TubeResult> t);
void render(Document &doc, Type<TubeSpec> t);
void render(Document &doc, Type<Weight> t);
void render(Document &doc, Type<Weights> t);
void render(Document &doc, Type<WordSpec> t);

/******************************************************************************/

inline Optional<Objective> request(Type<Objective> t, rebind::Variable const &v, rebind::Dispatch &msg) {
    if (auto p = v.request<Objective::Var>()) return Objective(std::move(*p));
    return msg.error("Cannot convert to Objective", t);
}

/******************************************************************************/

}
