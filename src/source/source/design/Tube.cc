#include <nupack/design/Tube.h>
#include <nupack/concentration/Equilibrate.h>

namespace nupack { namespace newdesign {

vec<StrandView> Tube::strand_types(vec<Complex> const &cs) const {
    auto comps = complexes(cs);
    std::set<StrandView> temp_strands;
    for (auto const &c : comps)
        for (auto strand : c.strands)
            temp_strands.emplace(strand);

    return vec<StrandView>(view(temp_strands));
}



/**
 * @brief fill in stoichiometry once complexes is stable
 * @details Uses comparison of the strand views within complexes to generate
 *     the stoichiometric matrix (vector of how many of each type of strand
 *     per complex)
 */
void Tube::compute_stoichiometry(vec<Complex> const &cs) {
    /* make vector of strand sequences in some order */
    auto strands = strand_types(cs);
    auto comps = complexes(cs);

    /* make matrix and vector appropriate sizes */
    real_mat stoich_mat;
    stoich_mat.zeros(len(comps), len(strands));

    auto find_index = [&](auto const &strand) {
        for (auto i : range(len(strands))) if (strands[i] == strand) return i;
        throw std::out_of_range("strand not found");
    };

    /* process complexes into stoich_mat matrix and vector */
    izip(comps, [&](auto i, auto const &c) {
        for (auto const &s : c.strands) {
            auto strand_index = find_index(s);
            stoich_mat(i, strand_index) += 1;
        }
    });

    stoichiometry = stoich_mat;
}


void remove_added_strands(uint num_strands, real_col &x) {
    if (num_strands > 0) x.shed_rows(x.n_rows - num_strands, x.n_rows - 1);
}

/**
 * @brief
 *
 * @param env compute environment allowing possible parallelization
 * @param map design object caching thermodynamic model(s) needed to evaluate
 * complex and tube properties
 * @param cs the set of all complexes in the design
 * @param s the sequence for complexes to evaluating themselves using
 * @param depth the depth in the decomposition tree of each complex in the tube
 * to evaluate estimates of the
 * @param part partition of complexes into active and passive sets
 * @return concentrations or concentration estimates of all complexes in tube in
 * same order as targets
 */
vec<real> Tube::concentrations(vec<real> const &log_pfuncs, EnsemblePartition const &part) const {
    auto water_conc = water_molarity(model.conditions.temperature);

    real_col dG, x0; real_mat A;
    // whether to use deflated mass constraints or not
    bool estimate = len(part) > 0 && any_of(targets, [&](auto const &t) {return !part.active(t.complex_index);});
    if (estimate) {
        std::tie(A, x0, dG) = deflate(log_pfuncs, part);
    } else {
        vec<real> vec_dG(indirect_view(targets, [&](auto const &t) {return -at(log_pfuncs, t.complex_index);}));
        dG = real_col(vec_dG);
        x0 = real_col(vmap<vec<real>>(targets, [&](auto const &c) {return c.target_conc / water_conc;}));
        A = stoichiometry;
    }

    real_col x = newdesign::concentrations(A, x0, dG);

    if (estimate) x = reinflate(x, part);
    return vmap<vec<real>>(x, [&](auto const &c) {return c * water_conc;});
}


/**
 * @brief compute the fraction (or estimated fraction) of total nucleotide
 * concentration represented by each of the complexes in the tube
 *
 * @param env compute environment allowing possible parallelization
 * @param map design object caching thermodynamic model(s) needed to evaluate
 * complex and tube properties
 * @param cs the set of all complexes in the design
 * @param s the sequence for complexes to evaluating themselves using
 * @param depth the depth in the decomposition tree of each complex in the tube
 * to evaluate estimates of the
 * @param part partition of complexes into active and passive sets
 * @return fractions or estimates of all complexes in tube in same order as
 * targets
 */
vec<real> Tube::fractions(vec<real> const &log_pfuncs, EnsemblePartition const &part) const {
    auto x = concentrations(log_pfuncs, part);
    for_each(x, [&](auto &i) {i /= nucleotide_concentration;});
    return x;
}


/**
 * @brief remove rows corresponding to passive complexes from stoichiometric
 * matrix and initial complex concentrations. Don't compute information for
 * passive complexes. Deflate total strand concentrations according to
 * parameters in partition.
 *
 * @param env compute environment allowing possible parallelization
 * @param map design object caching thermodynamic model(s) needed to evaluate
 * complex and tube properties
 * @param cs the set of all complexes in the design
 * @param s the sequence for complexes to evaluating themselves using
 * @param depth the depth in the decomposition tree of each complex in the tube
 * to evaluate estimates of the
 * @param part partition of complexes into active and passive sets
 * @return std::tuple<real_mat, real_col, real_col>
 */
std::tuple<real_mat, real_col, real_col> Tube::deflate(vec<real> const &log_pfuncs, EnsemblePartition const &part) const {
    auto water_conc = water_molarity(model.conditions.temperature);

    real_col x0 = real_col(vmap<vec<real>>(targets, [&](auto const &c) {return c.target_conc / water_conc;}));
    auto A = stoichiometry;

    vec<value_type_of<arma::uvec>> active;
    vec<uint> active_complex_inds;
    izip(targets, [&](auto i, auto const &t) {
        if (part.active(t.complex_index)) {
            active_complex_inds.emplace_back(t.complex_index);
            active.emplace_back(i);
        }
    });
    arma::uvec slice(active);

    vec<real> vec_dG(indirect_view(active_complex_inds, [&](auto i) {return -at(log_pfuncs, i);}));
    real_col ret_dG = real_col(vec_dG);

    real_mat ret_A = A.rows(slice);
    real_col ret_x0 = x0(slice);
    ret_x0 *= (1 - part.deflate);

    return std::make_tuple(std::move(ret_A), std::move(ret_x0), std::move(ret_dG));
}


/**
 * @brief Undo the dimension changes from deflation by adding in zero rows for
 * passive off-targets to the concentrations column
 *
 * @param x the solved concentrations in the deflated ensemble
 * @param part partition of complexes into active and passive sets
 * @return vector with positive concentrations for active complexes and zeroes
 * for passive complexes in the correct order
 */
real_col Tube::reinflate(real_col const &x, EnsemblePartition const &part) const {
    vec<value_type_of<arma::uvec>> active;
    izip(targets, [&](auto i, auto const &t) {if (part.active(t.complex_index)) active.emplace_back(i);});
    arma::uvec slice(active);

    real_col ret_x = la::zeros<real_col>(len(targets));
    ret_x(slice) = x;

    return ret_x;
}


/**
 * @brief compute the nucleotides' contributions to the structural defect
 * component of a complex's contribution to the test tube ensemble defect
 *
 * @param t the target information (including target concentration) for the
 * complex in the tube
 * @param comp_defect the nucleotide contributions to the complex ensemble
 * defect
 * @param concentration the actual concentration (or estimate) of the
 * corresponding complex
 * @return per-nucleotide contribution to structural defect of complex
 */
Defect structural_defect(TubeTarget const &t, Defect const &comp_defect, real concentration) {
    auto min_conc = min(concentration, t.target_conc);
    return Defect(vmap(comp_defect.contributions, [&](auto const &nd) {
        return std::make_pair(nd.first, min_conc * nd.second);
    }));
}


/* TODO: comment on non-collapsing version */
/**
 * @brief compute the nucleotides' contributions to the concentration defect component of a complex's contribution to the test tube ensemble defect
 *
 * @param t the target information (including target concentration) for the
 * complex in the tube
 * @param complexes the set of all complexes in the design
 * @param concentration the actual concentration (or estimate) of the
 * corresponding complex
 * @return per-nucleotide contribution to concentration defect of complex
 */
Defect concentration_defect(TubeTarget const &t, real concentration) {
    auto floor_diff = max(t.target_conc - concentration, 0);
    Defect defect;
    for (auto i : t.indices) defect.contributions.emplace_back(i, floor_diff);
    return defect;
}


/**
 * @brief compute the per-nucleotide (underlying variables) contributions to the
 * test tube ensemble defect (or estimate) on a concentration basis.
 *
 * @param env compute environment allowing possible parallelization
 * @param map design object caching thermodynamic model(s) needed to evaluate
 * complex and tube properties
 * @param cs the set of all complexes in the design
 * @param s the sequence for complexes to evaluating themselves using
 * @param depth the depth in the decomposition tree of each complex in the tube
 * to evaluate estimates of the
 * @param part partition of complexes into active and passive sets
 * @return the per-nucleotide test tube ensemble defect (or estimate)
 */
Defect Tube::defect(vec<real> const &log_pfuncs, vec<Defect> const &comp_defects,
        EnsemblePartition const &part, ComplexWeights const &weights) const {
    auto concs = concentrations(log_pfuncs, part);
    NUPACK_REQUIRE(targets.size(), ==, concs.size(), "Mismatch in number of specified concentrations", log_pfuncs.size(), comp_defects.size());

    std::map<uint, real> defects;
    auto accum = [&](auto const &def) {
        for (auto x : def.contributions) {
            auto it = defects.find(x.first);
            if (it == end_of(defects)) it = defects.emplace(x.first, 0).first;
            it->second += x.second;
        }
    };

    zip(targets, concs, [&](auto const &t, auto const &conc) {
        if (t.is_on_target()) {
            auto const &comp_def = at(comp_defects, t.complex_index);
            /* without weights */
            if (len(weights) == 0) {
                accum(structural_defect(t, comp_def, conc).reduced());
                accum(concentration_defect(t, conc).reduced());
            }
            /* with weights */
            else {
                auto str_def = structural_defect(t, comp_def, conc).contributions;
                auto conc_def = concentration_defect(t, conc).contributions;
                auto const &comp_weights = weights.at(t.complex_index);

                NUPACK_REQUIRE(len(str_def), ==, len(conc_def));
                NUPACK_REQUIRE(len(str_def), ==, len(comp_weights));
                NUPACK_REQUIRE(len(str_def), ==, len(t.indices), "all weighted defects must be non-collapsed for weighting");

                Defect tot_def;
                zip(str_def, conc_def, comp_weights, [&](auto s, auto c, auto w) {
                    NUPACK_REQUIRE(s.first, ==, c.first, "must be same underlying nucleotides in the same order");
                    tot_def.contributions.emplace_back(s.first, w * (s.second + c.second));
                });
                accum(tot_def.reduced());
            }
        }
    });


    // repackage vector of defects into vector of pairs from nucleotides to defects
    defect_vec defs(view(defects));
    // izip(mapped_defects, [&](auto i, auto d) {if (d > 0) defs.emplace_back(i, d);});
    return Defect(defs);
}


/**
 * @brief compute the per-nucleotide (underlying variables) contributions to the
 * test tube ensemble defect (or estimate) scaled by total nucleotide concentration.
 *
 * @param env compute environment allowing possible parallelization
 * @param map design object caching thermodynamic model(s) needed to evaluate
 * complex and tube properties
 * @param cs the set of all complexes in the design
 * @param s the sequence for complexes to evaluating themselves using
 * @param depth the depth in the decomposition tree of each complex in the tube
 * to evaluate estimates of the
 * @param part partition of complexes into active and passive sets
 * @return the normalized per-nucleotide test tube ensemble defect (or estimate)
 */
Defect Tube::normalized_defect(vec<real> const &log_pfuncs, vec<Defect> const &comp_defects,
        EnsemblePartition const &part, ComplexWeights const &weights) const {
    auto def = defect(log_pfuncs, comp_defects, part, weights);
    for (auto & c : def.contributions) c.second /= nucleotide_concentration;
    return def;
}


vec<real> Tube::concentrations(vec<real> const &log_pfuncs) const {
    real_mat A = stoichiometry;

    auto water_conc = water_molarity(model.conditions.temperature);
    real_col x0(vmap<vec<real>>(targets, [&](auto const &c) { return c.target_conc / water_conc; }));
    real_col dG(vmap<vec<real>>(targets, [&](auto const &c) { return -at(log_pfuncs, c.complex_index); }));

    real_col x = newdesign::concentrations(A, x0, dG);

    x *= water_conc;
    return {begin_of(x), end_of(x)};
}


/**
 * @brief compute the concentrations of complexes in a test tube given the their
 * stoichiometry, initial fractions and free energies
 *
 * @param A the stoichiometry matrix
 * @param x0 the initial mole fractions (from target concentrations)
 * @param dG the unscaled free energy (-log(pfunc))
 * @return real_col the mole fractions of each of the complexes
 */
real_col concentrations(real_mat const &A, real_col const &x0, real_col const &dG) {
    real_col logq = -dG;
    concentration::Options options;
    options.method = concentration::Method::cd;

    concentration::Output<real> sol;
    std::string exception;

    try {
        sol = concentration::equilibrate(A, la::log(A.t() * x0), logq, options);
    } catch (std::exception const &e) {
        exception = e.what();
    }

    if (exception.empty() && (sol.converged || sol.error <= 1e-3)) return std::move(sol.solution);

    json info{{"A", A}, {"x0", x0}, {"logq", logq},
              {"options", options}, {"solution", sol}, {"exception", exception}};
    NUPACK_ERROR("nupack::design: equilibrium convergence failure: " + info.dump());
}


/**
 * @brief
 *
 * @param log_pfuncs
 * @return vec<real>
 */
vec<real> Tube::fractions(vec<real> const &log_pfuncs) const {
    auto x = concentrations(log_pfuncs);
    for_each(x, [&](auto &i) {i /= nucleotide_concentration;});
    return x;
}


/**
 * @brief save indices for each complex so that design.complexes doesn't need to
 * be passed in to compute concentration_defect
 *
 * @param cs
 */
void Tube::store_complex_indices(vec<Complex> const &cs) {
    for (auto &t : targets) t.indices = at(cs, t.complex_index).to_indices();
}

}}
