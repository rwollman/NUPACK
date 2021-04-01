#include <nupack/design/OutputResult.h>
#include <nupack/design/Designer.h>
#include <nupack/reflect/Serialize.h>
#include <limits>
#include <cmath>

namespace nupack { namespace newdesign {

// DesignResult::DesignResult(Specification const &spec, Designer const &designer, Result const &result) :
DesignResult::DesignResult(Designer const &designer) :
        model(at(designer.design.complexes, 0).target.model), parameters(designer.parameters),
        stats(designer.stats), objectives(designer.objectives),
        success(designer.success()),
        results({designer.best.full.full_evaluation(designer)}), weights(designer.weights) {}
        // results(indirect_view(designer.archive.full.results, [&](auto const &res) {
        //     return res.evaluated;
        // })) {}


SingleResult::SingleResult(Designer const &designer, Result const &res) :
        defects(res.totals()), weighted_defects(res.weighted_totals()) {
    auto const &design = designer.design;
    auto const &seqs = design.sequences;
    auto const &sequence = res.sequence;
    auto env = Local();

    /* domains */
    for (auto const &el : seqs.domains) { domains.emplace(el.first, el.second.to_sequence(sequence)); }
    /* strands */
    for (auto const &el : seqs.strands) { strands.emplace(el.first, el.second.to_sequence(sequence)); }

    auto const &models = design.models;

    auto & engobs = const_cast<EngineObserver &>(designer.obs);
    /* complexes */
    for (auto const &d : design.complexes) {
        ComplexResult comp;
        comp.name = d.name;
        comp.sequence = to_nick_sequence(d.strands, sequence);
        comp.structure = d.target.structure;
        auto x = d.log_pfunc(env, models, sequence, 0, {}, engobs);
        x = std::isfinite(x) ? x : std::numeric_limits<decltype(x)>::lowest();
        comp.log_partition_function = x;
        if (d.is_on_target()) {
            comp.pair_probabilities = d.pair_probabilities(env, models, sequence, 0, {}, engobs);
        }
        comp.defect = d.defect(env, models, sequence, 0, {}, engobs).total();
        comp.normalized_defect = comp.defect / len(d);

        complexes.emplace_back(comp);
    }

    auto log_pfuncs = design.log_pfuncs(env, 0, {}, {}, engobs);
    auto complex_defects = design.complex_defects(env, 0, {}, {}, engobs);

    /* tubes */
    for_each(design.tubes, [&](auto const &t) {
        TubeResult tube;
        tube.name = t.name;
        tube.nucleotide_concentration = t.nucleotide_concentration;
        tube.defect = t.defect(log_pfuncs, complex_defects).total();
        tube.normalized_defect = tube.defect / tube.nucleotide_concentration;

        auto concentrations = t.concentrations(log_pfuncs);

        zip(t.targets, concentrations, [&](auto const &c, auto conc) {
            auto const &ref = at(complexes, c.complex_index);
            TubeComplex comp;
            comp.name = ref.name;
            comp.concentration = conc;
            comp.target_concentration = c.target_conc;

            auto comp_defect = at(complex_defects, c.complex_index);
            comp.structural_defect = structural_defect(c, comp_defect, comp.concentration).total();
            comp.concentration_defect = concentration_defect(c, comp.concentration).total();
            comp.defect = comp.structural_defect + comp.concentration_defect;
            comp.normalized_defect_contribution = comp.defect / tube.nucleotide_concentration;

            tube.complexes.emplace_back(comp);
        });

        tubes.emplace_back(tube);
    });
}





}}
