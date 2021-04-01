#include "Design.h"

namespace nupack::newdesign {

void render(Document &doc, Type<Timer> t) {
    doc.type(t, "design.components.Timer");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "start", &Timer::start);
    doc.method(t, "elapsed", &Timer::elapsed);
    doc.method(t, "stop", &Timer::stop);
}

void render(Document &doc, Type<DesignStats> t) {
    doc.type(t, "design.results.Stats");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<EnsemblePartition> t) {
    doc.type(t, "design.results.Partition");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<MultitubeObjective> t) {
    doc.type(t, "design.objectives.MultitubeObjective");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<TubeObjective> t) {
    doc.type(t, "design.objectives.TubeObjective");
    doc.method(t, "new", rebind::construct<string>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ComplexObjective> t) {
    doc.type(t, "design.objectives.ComplexObjective");
    doc.method(t, "new", rebind::construct<string>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<SSMObjective> t) {
    doc.type(t, "design.objectives.SSMObjective");
    doc.method(t, "new", rebind::construct<vec<string>, uint>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<PatternObjective> t) {
    doc.type(t, "design.objectives.PatternObjective");
    doc.method(t, "new", rebind::construct<vec<string>, vec<Sequence>>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<SimilarityObjective> t) {
    doc.type(t, "design.objectives.SimilarityObjective");
    doc.method(t, "new", rebind::construct<vec<string>, vec<Sequence>, vec<std::pair<real, real>>>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<EnergyEqualizationObjective> t) {
    doc.type(t, "design.objectives.EnergyEqualizationObjective");
    doc.method(t, "new", rebind::construct<vec<string>, Optional<real>>(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<Objective> t) {
    doc.type(t, "design.objectives.Objective");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}



void render(Document &doc, Type<SingleResult> t) {
    doc.type(t, "design.results.Single");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ComplexResult> t) {
    doc.type(t, "design.results.Complex");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<TubeComplex> t) {
    doc.type(t, "design.results.TubeComplex");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<TubeResult> t) {
    doc.type(t, "design.results.Tube");

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<DesignResult> t) {
    doc.type(t, "design.results.Result");
    doc.method(t, "new", rebind::construct(t));
    render_json(doc, t);
    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<DomainSpec> t) {
    doc.type(t, "design.components.Domain");
    doc.method(t, "new", rebind::construct<string, string>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<StrandSpec> t) {
    doc.type(t, "design.components.Strand");
    doc.method(t, "new", rebind::construct<string, vec<string>>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ComplexSpec> t) {
    doc.type(t, "design.components.Complex");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<TubeSpec> t) {
    doc.type(t, "design.components.Tube");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

// Start constraint specs

void render(Document &doc, Type<DualListSpec> t) {
    doc.type(t, "design.components.DualList");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<PatternSpec> t) {
    doc.type(t, "design.components.Pattern");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<DiversitySpec> t) {
    doc.type(t, "design.components.Diversity");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<WordSpec> t) {
    doc.type(t, "design.components.Word");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<SimilaritySpec> t) {
    doc.type(t, "design.components.Similarity");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<ConstraintSpec> t) {
    doc.type(t, "design.components.Constraints");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


using Condition = rebind::Callback<rebind::Integer>; // This has to be so for now, bindings are a bit broken if using a non-primitive type.
using Handler = rebind::Callback<void>;

struct DesignRunner {
    auto operator()(Specification const &spec, Local const &env, Condition condition, Handler handler, std::optional<DesignResult> restart_={}) const {
        SignalRuntime l;

        if (handler.function && !condition.function) {
            NUPACK_ERROR("If using checkpointing with designer, you must supply a checkpoint condition");
        }

        Designer d(spec);
        d.initialize();

        if (bool(restart_)) {
            try {
                auto const &restart = restart_.value();
                auto const &res = restart.results[0];
                auto mapping = Specification::ensure_compatibility(spec, res); // throw useful errors where there are mismatches
                auto &seqs = d.design.sequences;

                for (auto const &domain: res.domains) seqs.set_domain(domain.first, domain.second);

                d.Psi.mask = vmap<decltype(d.Psi.mask)>(mapping, [&](auto i) {return restart.stats.final_Psi.active(i);});
                d.stats = restart.stats;
                d.redecompose_active(env, 0);
            } catch (...) {
                print("nupack: Failure in loading Design specification from intermediate result. Does it correspond to the correct design?");
                throw;
            }
        }
        if (condition.function) {
            auto real_condition = [c=std::move(condition)] (Designer const &des, bool done) {
                switch (c(des.stats, des.timer, done)) {
                    case +1: return true;
                    case -1: throw SignalError::sigint();
                    default: return false;
                }
            };
            if (handler.function) {
                auto real_handler = [&, h=std::move(handler)] (Designer &des) {
                    des.stats.design_time += des.timer.stop();
                    des.stats.final_Psi = des.Psi;
                    auto sequence = des.best_sequence(env);
                    auto result = DesignResult(des);
                    h(std::move(result));
                    des.timer.start();
                };
                d.checkpoint = [handler=std::move(real_handler),
                    condition=std::move(real_condition)] (Designer &des, bool done) {if (condition(des, done)) handler(des);};
            } else {
                d.checkpoint = [condition=std::move(real_condition)] (Designer &des, bool done) {condition(des, done);};
            }
        }

        d.optimize_tubes(env);

        // d.design.set_sequence(d.best.full.sequence);
        return DesignResult(d);
    }
};



void render(Document &doc, Type<Specification> t) {
    doc.type(t, "design.core.Specification");
    doc.method(t, "new", rebind::construct<Model<real> const &, bool>(t));
    render_json(doc, t);

    doc.method<4>(t, "()", DesignRunner());
    doc.method(t, "evaluate", [](Specification const &spec, Local const &env) {
        Designer d(spec);

        auto & seqs = d.design.sequences;

        auto condition = seqs.all_nucleotides_fixed();
        if (!condition.first) NUPACK_ERROR("there are variable nucleotides in the design in domain: " + condition.second);

        d.initialize();
        d.time_analysis(env);
        d.best.full = d.evaluate_objectives(env, 0, {}, d.weights);
        d.best.full.full_evaluation(d);

        return DesignResult(d);
    });

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<Weight> t) {
    doc.type(t, "design.weights.Weight");
    doc.method(t, "new", rebind::construct<
            Optional<string>, Optional<string>, Optional<string>, Optional<string>,
            real>(t));

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}


void render(Document &doc, Type<ReversedComplex> t) {
    doc.type(t, "design.weights.ReversedComplex");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "domains", &ReversedComplex::domains);
    doc.method(t, "strands", &ReversedComplex::strands);

    doc.method(t, "{}", dumpable(t));
}


void render(Document &doc, Type<Weights> t) {
    doc.type(t, "design.weights.Weights");
    doc.method(t, "new", rebind::construct(t));
    doc.method(t, "add", &Weights::add);
    doc.method(t, "add_objective_weight", &Weights::add_objective_weight);

    doc.method(t, "{}", dumpable(t));
    render_public(doc, t);
}

void render(Document &doc, Type<DesignParameters> t) {
    using base_type = DesignParameters;
    doc.type(t, "design.components.Parameters");
    doc.method(t, "new", rebind::construct(t));

    doc.method(t, "{}", dumpable(t));
    NUPACK_PUBLIC(DesignParameters, rng_seed, f_stop, f_passive, H_split, N_split, f_split,
            f_stringent, dG_clamp, M_bad, M_reseed, M_reopt, f_redecomp, f_refocus,
            cache_bytes_of_RAM, f_sparse, slowdown, log, decomposition_log, thermo_log,
            time_analysis);
}

}

namespace nupack {

void render_design(Document &doc) {
    doc.render<newdesign::MultitubeObjective>();
    doc.render<newdesign::TubeObjective>();
    doc.render<newdesign::ComplexObjective>();
    doc.render<newdesign::SSMObjective>();
    doc.render<newdesign::SimilarityObjective>();
    doc.render<newdesign::EnergyEqualizationObjective>();
    doc.render<newdesign::PatternObjective>();
    doc.render<newdesign::Objective>();

    doc.render<newdesign::Specification>();
    doc.render<newdesign::Timer>();

    doc.render<newdesign::DesignResult>();

    doc.render<newdesign::EnsemblePartition>();
}

}
