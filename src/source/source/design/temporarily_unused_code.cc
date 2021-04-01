
 Result Designer::alternate_optimize_tubes(Local const &env) {
     timer = Timer().start();
     logger.log("time", "type", "depth", "psi_active", "psi_passive", "sequence", "defect");

     /******************************************************************/
     max_depth = design.max_depth();

     auto estimate = optimize_forest(env, design.sequence());
     design.set_sequence(estimate.sequence);

     bool done = false;
     Result full {inf_result};

     while (!done) {
         checkpoint(*this);
         // if (length_extrapolation_refocus(env)) {
         if (sum_pf_refocus(env)) {
             estimate = evaluate_objectives(env, 0, Psi);
             // {design.sequence(), {design.normalized_defect(env, 0, Psi)}};
             full = evaluate_objectives(env, 0, {});
             // {design.sequence(), {design.normalized_defect(env)}};
             done = full.total() <= max(parameters.f_stop, estimate.total());

             if (full.total() < best.full.total()) {
                 logger.log(timer.elapsed(), "root accepted", 0, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(full.sequence), LAZY(full.total()));
                 best.full = full;
             } else {
                 logger.log(timer.elapsed(), "root rejected", 0, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(full.sequence), LAZY(full.total()));
             }

             if (!done) refocus(env, full, estimate);
         }

         if (!done) {
             estimate = optimize_forest(env, full.sequence);
             design.set_sequence(estimate.sequence);
         }
     }

     // Result full {design.sequence(), design.normalized_defect(env)};
     // if (full.total() < best.full.total()) best.full = full;
     // logger.log(timer.elapsed(), "root accepted", 0, Psi.num_active(), Psi.num_inactive(), design.sequences.json_domains(full.sequence), LAZY(full.total()));

     // while (full.total() > max(parameters.f_stop, estimate.total())) {

     // }
     /******************************************************************/
     stats.design_time += timer.stop(); // if checkpointed and restarted, the += will make the output stats reflect the total design time instead of just the most recent segment
     stats.final_Psi = Psi;

     time_analysis(env);
     BEEP(design.models.ram());
     // std::cout << json(design.sequences.correlation_matrix) << std::endl;

     return best.full;
 }

/**
 * @brief
 *
 * @param env
 * @return true if first off-target added barely nudged
 * @return false if more than one off-target were added before curve levelled out
 */
 bool Designer::length_extrapolation_refocus(Local const &env) {
     /* exit immediately if there are no off-targets remaining to add */
     if (Psi.all_active()) return true;

     /* compute pfuncs for active */
     vec<real> log_pfuncs = vmap(Psi.actives(), [&](auto i) {
         return at(design.complexes, i).log_pfunc(env, design.models, design.sequence(), 0);
     });
     /* least squares to get predictor */
     vec<real> lengths = vmap<vec<real>>(Psi.actives(), [&](auto i) {
         return len(at(design.complexes, i));
     });
     auto coefficients = ord_lin_lsq(lengths, log_pfuncs);
     /* predict passive */
     log_pfuncs = vec<real>(len(design.complexes), 0.0);
     izip(design.complexes, [&](auto i, auto const &c) {
         at(log_pfuncs, i) = Psi.active(i) ?
             c.log_pfunc(env, design.models, design.sequence(), 0) :
             coefficients[0] + len(c) * coefficients[1]; // regression prediction
     });

     /* compute estimated fractions with prediction */
     vec<real> fractions(len(design.complexes), 0.0);
     for (auto const &tube : design.tubes) {
         zip(tube.targets, tube.fractions(log_pfuncs, design.complexes), [&] (auto const &c, auto frac) {
             auto i = c.complex_index;
             if (!Psi.active(i)) at(fractions, i) += frac;
         });
     }

     vec<std::pair<uint, real>> passive;
     izip(Psi.mask, [&](auto i, auto b) {if (!b) passive.emplace_back(i, at(fractions, i));});
     sort(passive, [](auto const &a, auto const &b) {return a.second > b.second;});

     auto order = key_view(passive);
     auto cur = begin_of(order);
     auto part = Psi;
     /* add in order (same as normal refocus); stopping when new vs old is in agreement */

     at(part.mask, *cur) = true;

     Result prev = evaluate_objectives(env, 0, Psi);
     // {design.sequence(), {design.normalized_defect(env, 0, Psi)}};
     Result estimate = evaluate_objectives(env, 0, part);
     // {design.sequence(), {design.normalized_defect(env, 0, part)}};
     logger.log(timer.elapsed(), "refocused", 0, part.num_active(), part.num_inactive(), design.sequences.json_domains(estimate.sequence), LAZY(estimate.total()));

     auto condition = [&]() -> bool { return (estimate.total() - prev.total())/prev.total() < parameters.f_refocus; };
     bool immediate = condition(); /* if most likely off-target didn't change defect much, then evaluate full ensemble */

     while (++cur != end_of(order) && !condition()) {
         // if (cur >= end_of(order)) break;
         at(part.mask, *cur) = true;
         prev = estimate;
         estimate = evaluate_objectives(env, 0, part);
         // {design.sequence(), {design.normalized_defect(env, 0, part)}};
         logger.log(timer.elapsed(), "refocused", 0, part.num_active(), part.num_inactive(), design.sequences.json_domains(estimate.sequence), LAZY(estimate.total()));
     }

     vec<uint> changed;
     izip(part.mask, Psi.mask, [&](auto i, auto n, auto o) {if (n && !o) changed.emplace_back(i);});
     subset_decompose(changed);
     stats.offtargets_added_per_refocus.emplace_back(len(changed));

     Psi = part;
     known_bads.clear();
     return immediate;
 }

 bool Designer::sum_pf_refocus(Local const &env) {
     /* exit immediately if there are no off-targets remaining to add */
     if (Psi.all_active()) return true;

     /* predict passive */
     auto log_pfuncs = vec<real>(len(design.complexes), 0.0);
     izip(design.complexes, [&](auto i, auto const &c) {
         at(log_pfuncs, i) = Psi.active(i) ?
             c.log_pfunc(env, design.models, design.sequence(), 0) :
             c.log_pf_single_strands(env, design.models, design.sequence()); // single stranded prediction
     });

     /* compute estimated fractions with prediction */
     vec<real> fractions(len(design.complexes), 0.0);
     for (auto const &tube : design.tubes) {
         zip(tube.targets, tube.fractions(log_pfuncs, design.complexes), [&] (auto const &c, auto frac) {
             auto i = c.complex_index;
             if (!Psi.active(i)) at(fractions, i) += frac;
         });
     }

     vec<std::pair<uint, real>> passive;
     izip(Psi.mask, [&](auto i, auto b) {if (!b) passive.emplace_back(i, at(fractions, i));});
     sort(passive, [](auto const &a, auto const &b) {return a.second > b.second;});

     auto order = key_view(passive);
     auto cur = begin_of(order);
     auto part = Psi;
     /* add in order (same as normal refocus); stopping when new vs old is in agreement */

     at(part.mask, *cur) = true;

     Result prev = evaluate_objectives(env, 0, Psi);
     // {design.sequence(), {design.normalized_defect(env, 0, Psi)}};
     Result estimate = evaluate_objectives(env, 0, part);
     // {design.sequence(), {design.normalized_defect(env, 0, part)}};
     logger.log(timer.elapsed(), "refocused", 0, part.num_active(), part.num_inactive(), design.sequences.json_domains(estimate.sequence), LAZY(estimate.total()));

     auto condition = [&]() -> bool { return (estimate.total() - prev.total())/prev.total() < parameters.f_refocus; };
     bool immediate = condition(); /* if most likely off-target didn't change defect much, then evaluate full ensemble */

     while (++cur != end_of(order) && !condition()) {
         // if (cur >= end_of(order)) break;
         at(part.mask, *cur) = true;
         prev = estimate;
         estimate = evaluate_objectives(env, 0, part);
         // {design.sequence(), {design.normalized_defect(env, 0, part)}};
         logger.log(timer.elapsed(), "refocused", 0, part.num_active(), part.num_inactive(), design.sequences.json_domains(estimate.sequence), LAZY(estimate.total()));
     }

     vec<uint> changed;
     izip(part.mask, Psi.mask, [&](auto i, auto n, auto o) {if (n && !o) changed.emplace_back(i);});
     subset_decompose(changed);
     stats.offtargets_added_per_refocus.emplace_back(len(changed));

     Psi = part;
     known_bads.clear();
     return immediate;
 }
