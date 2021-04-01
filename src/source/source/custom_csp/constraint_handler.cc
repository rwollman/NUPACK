#include <nupack/design/custom_csp/constraint_handler.h>
#include <nupack/design/custom_csp/sequence_utils.h>

#include <nupack/design/custom_csp/design_debug.h>
#include <nupack/design/custom_csp/algorithms.h>

#include <vector>
#include <array>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>


namespace nupack {
namespace custom_csp {
VariableNode::VariableNode() :
    depth(0), cost(0), tiebreaker(random_float()), parent(nullptr) {}

VariableNode::VariableNode(std::shared_ptr<VariableNode> parent) :
    depth(parent->get_depth() + 1), cost(0),
    tiebreaker(random_float()), parent(parent) {}

/* modify incoming Domains based on changes in stack */
bool VariableNode::add_implications(AllowTable & allow_table,
                                    const SolveStack & stack) {
    v.reserve(stack.size() + v.size());
    bool success = true;
    for (auto const &c : stack.v) {
        NUPACK_DEBUG_CHECK(allow_table.size() > c.var,
                           "Invalid variable being set: " + to_string(c.var)
                           + " >= " + to_string(allow_table.size()));
        NUPACK_DEBUG_CHECK(allow_table[c.var].size() > c.val,
                           "invalid value being set for " + to_string(c.var) + ": " + to_string(c.val)
                           + " >= " + to_string(allow_table[c.var].size()));

        if (allow_table[c.var][c.val] == NUPACK_VV_UNSET) {
            allow_table[c.var][c.val] = c.trit;
            v.emplace_back(c.var, c.val, c.trit);
        } else if (allow_table[c.var][c.val] != c.trit) {
            success = false;
            break;
        }
    }

    return success;
}



void VariableNode::change_branch(std::shared_ptr<VariableNode> from,
                                 std::shared_ptr<VariableNode> to, AllowTable & allow_table) {
    std::shared_ptr<VariableNode> c_from = from;
    std::shared_ptr<VariableNode> c_to = to;
    vec<std::shared_ptr<VariableNode> > to_stack;

    if (!c_to) {
        while (c_from) {
            c_from->rollback_variables(allow_table);
            c_from = c_from->get_parent();
        }
    } else if (!c_from) {
        while (c_to) {
            c_to->assign_variables(allow_table);
            c_to = c_to->get_parent();
        }
    } else {
        while (!c_to->is_root() && c_to->get_depth() > c_from->get_depth()) {
            to_stack.push_back(c_to);
            c_to = c_to->get_parent();
        }

        while (!c_from->is_root() && c_from->get_depth() > c_to->get_depth()) {
            c_from = c_from->get_parent();
        }

        if (c_from && c_to) {
            NUPACK_CHECK(c_from->get_depth() == c_to->get_depth(),
                         "Invalid depths in change_branch after equal depth. "
                         "this->depth = " + to_string(c_from->get_depth())
                         + " other->depth = " + to_string(c_to->get_depth()));

            while (!c_from->is_root() && !c_to->is_root() && c_from != c_to) {
                NUPACK_DEBUG_CHECK(c_to->get_depth() == c_from->get_depth(),
                                   "Depths out of sync in change_branch "
                                   "this->depth = " + to_string(c_from->get_depth()) +
                                   "  other->depth = " + to_string(c_to->get_depth()));

                to_stack.push_back(c_to);

                c_to = c_to->get_parent();
                c_from = c_from->get_parent();
            }

            NUPACK_CHECK(c_from && c_to && c_from == c_to, "Disjoint tree found");

            // NUPACK_DEBUG_PRINT("Rolling back from " << from->get_depth());
            // NUPACK_DEBUG_PRINT("Rolling back to " << c_from->get_depth());
            c_from = from;
            while (c_from != c_to) {
                c_from->rollback_variables(allow_table);
                c_from = c_from->get_parent();
            }

            // NUPACK_DEBUG_PRINT("Rolling forward from " << c_to->get_depth());
            // NUPACK_DEBUG_PRINT("Rolling forward to " << to->get_depth());
            for (auto const &r : reversed(to_stack)) r->assign_variables(allow_table);
        }
    }
}



void VariableNode::rollback_variables(AllowTable & allow_table) {
    // NUPACK_DEBUG_PRINT("BEFORE ROLLBACK");
    // print_table(allow_table, std::cout);
    for (auto const &c : v) {
        NUPACK_DEBUG_CHECK((int)allow_table[c.var].size() > c.val,
                           "usize-value mismatch at " + to_string(c.var) + ": "
                           + to_string(allow_table[c.var].size())
                           + " <= " + to_string(c.val));

        NUPACK_CHECK((int)allow_table[c.var][c.val] == c.trit,
                     "Invalid rollback at variable: " + to_string(c.var) +
                     " value: " + to_string(c.val) + " allowed: "
                     + to_string((int)c.trit) + "  assigned: " + to_string((int)allow_table[c.var][c.val]));

        allow_table[c.var][c.val] = NUPACK_VV_UNSET;
    }
    // NUPACK_DEBUG_PRINT("AFTER ROLLBACK");
    // print_table(allow_table, std::cout);
}



void VariableNode::assign_variables(AllowTable & allow_table) {
    for (auto const &c : v) {
        NUPACK_DEBUG_CHECK((int)allow_table[c.var].size() > c.val,
                           "usize-value mismatch at " + to_string(c.var) + ": "
                           + to_string(allow_table[c.var].size())
                           + " <= " + to_string(c.val));

        NUPACK_CHECK(allow_table[c.var][c.val] == NUPACK_VV_UNSET,
                     "Assigning doubly to variable: " + to_string(c.var) +
                     " value: " + to_string(c.val) + " allowed: "
                     + to_string(c.trit) + "  old_all: " + to_string((int)allow_table[c.var][c.val]));

        allow_table[c.var][c.val] = c.trit;
    }
}



void VariableNode::update_cost(const vec<int> & orig) {
    int new_cost = 0;

    for (auto const &c : v) {
        NUPACK_DEBUG_CHECK(c.var >= 0 && c.var < orig.size(),
                           "Original size (" + to_string(orig.size()) + ") variable ("
                           + to_string(c.var) + ") mismatch");
        if (orig[c.var] == c.val && c.trit == NUPACK_VV_FALSE) {
            new_cost++;
        }
    }

    n_assigned = v.size();
    if (parent) {
        new_cost += parent->get_cost();
        n_assigned += parent->get_n_assigned();
    }

    cost = new_cost;
}


void SolveStruc::init_weights(AllowTable const &table) {
    NUPACK_CHECK(len(table) == len(value_allowed), "incoming AllowTable is different size than solver AllowTable");
    weight = vmap(table, [](auto const &c) {return double(ConstraintHandler::get_n_unset(c));});
}


bool CompConstraint::propagate_constraint(int modified, SolveStack & sstack,
        const SolveStruc & ss) const {
    constexpr int n_bases = 4;
    constexpr std::array<int, n_bases> base =   {{0, 1, 2, 3}};
    constexpr std::array<int, n_bases> comp =   {{3, 2, 1, 0}};
    constexpr std::array<int, n_bases> w_comp = {{ -1, -1, 3, 2}};

    bool satisfied = true;
    if (strength != NUPACK_CS_NONE) {
        int i_var = modified;
        int j_var = -1;
        if (i_var == i) {
            j_var = j;
        } else if (i_var == j) {
            j_var = i;
        }
        if (j_var >= 0) {
            vec<trinary> poss(n_bases, false);
            for (auto i = 0; i < n_bases; i++) {
                if (ss.value_allowed[i_var][i]) {
                    poss[comp[base[i]]] = true;
                    if (strength == NUPACK_CS_WEAK && w_comp[base[i]] >= 0) {
                        poss[w_comp[base[i]]] = true;
                    }
                }
            }

            for (auto i = 0; i < n_bases; i++) {
                if (!poss[i]) {
                    if (ss.value_allowed[j_var][i] == NUPACK_VV_UNSET) {
                        sstack.push_back(j_var, i, NUPACK_VV_FALSE);
                    } else if (ss.value_allowed[j_var][i] == NUPACK_VV_TRUE) {
                        satisfied = false;
                    }
                }
            }
        }
    }
    return satisfied;
}



bool IdentConstraint::propagate_constraint(int modified, SolveStack & sstack,
        const SolveStruc & ss) const {
    int i_var = modified;
    int j_var = -1;
    bool success = true;
    if (i_var == i) {
        j_var = j;
    } else if (i_var == j) {
        j_var = i;
    }
    if (j_var >= 0) {
        vec<trinary> poss(ss.value_allowed[i_var]);

        constexpr int n_bases = 4;
        for (auto i = 0; i < n_bases; i++) {
            if (ss.value_allowed[j_var][i] == NUPACK_VV_TRUE && !poss[i]) {
                success = false;
            } else if (ss.value_allowed[j_var][i] == NUPACK_VV_UNSET && !poss[i]) {
                sstack.push_back(j_var, i, NUPACK_VV_FALSE);
            }
        }
    }

    return success;
}



PatternConstraint::PatternConstraint(vec<int> const & vars, const string & constraint,
        const vec<int> & poss_nucs) :
        nuc_ids(vars), constraint(constraint),
        pattern(SequenceUtils::nucs_to_bools(SequenceUtils::str_to_nuc(constraint))) {

    for (auto i = 0; i < nuc_ids.size(); i++) nuc_id_map[nuc_ids[i]] = i;

    // If pattern is longer than target, there is effectively no constraint
    if (pattern.size() > nuc_ids.size()) {
        starts = {};
        return;
    }

    for (auto i = 0; i < nuc_ids.size() - ((int)pattern.size() - 1); i++) {
        trinary allowed = false; // is breaking the pattern even allowed?
        for (auto j = 0; j < pattern.size(); j++) {
            int nuc_id = nuc_ids[i + j];
            auto callowed = SequenceUtils::nuc_to_bool(poss_nucs[nuc_id]);
            NUPACK_DEBUG_CHECK(pattern[j].size() == callowed.size(), "invalid callowed size");

            for (auto k = 0; k < pattern[j].size(); k++) {
                if (callowed[k] && !pattern[j][k]) allowed = true;
            }
        }
        starts.push_back(allowed);
    }
}



bool PatternConstraint::propagate_constraint(int modified, SolveStack & sstack,
        const SolveStruc & ss) const {
    // trivial constraint, see constructor
    if (starts.empty()) return true;

    int modified_index = nuc_id_map.find(modified)->second;
    int max_index = nuc_ids.size() - (pattern.size() - 1);
    int min_index = 0;

    int tmp_max = modified_index + pattern.size() - 1;
    int tmp_min = modified_index - pattern.size() + 1;
    max_index = std::min(tmp_max, max_index);
    min_index = std::max(tmp_min, min_index);


    bool success = true;
    if (min_index > max_index) return true;

    for (auto i = min_index; i == min_index || i < max_index; i++) {
        if (starts[i]) {
            // Number of nucleotides that have possibilities to break the match
            auto n_nuc_choices = 0;
            auto nuc_choice = -1;
            for (auto j = 0; j < pattern.size(); j++) {
                auto c_nuc = nuc_ids[i + j];
                auto c_nuc_choices = false;
                for (auto k = 0; k < pattern[j].size(); k++) {
                    if (ss.value_allowed[c_nuc][k] && !pattern[j][k]) {
                        // value allowed, but not in pattern
                        if (!c_nuc_choices) {
                            // mark the current
                            c_nuc_choices = true;
                            nuc_choice = j;
                            n_nuc_choices++;
                            break;
                        }
                    }
                }

                if (n_nuc_choices >= 2) {
                    break;
                }
            }

            if (n_nuc_choices == 0) {
                success = false;
            } else if (n_nuc_choices == 1) {
                auto c_nuc = nuc_ids[i + nuc_choice];
                for (auto k = 0; k < pattern[nuc_choice].size(); k++) {
                    if (ss.value_allowed[c_nuc][k] && pattern[nuc_choice][k]) {
                        NUPACK_DEBUG_CHECK(ss.value_allowed[c_nuc][k] == NUPACK_VV_UNSET,
                                           "Invalid state to change to false "
                                           + to_string((int)ss.value_allowed[c_nuc][k]));
                        sstack.push_back(c_nuc, k, NUPACK_VV_FALSE);
                    }
                }
            }
        }
        else {
            success = false;
        }
    }

    return success;
}



WordConstraint::WordConstraint(const vec<int> & vars,
                               const vec<string> & constraint, int supp_var) :
    supp_var(supp_var), nuc_ids(vars) {
    varval_to_ids.resize(vars.size());
    for (auto & v : varval_to_ids) v.resize(4);

    auto i_val = 0;
    for (auto & con : constraint) {
        NUPACK_CHECK(vars.size() == con.size(), "usize mismatch in WordConstraint");
        AllowTable cur_con;
        allowed_ind.emplace_back();
        vec<int> int_pattern = SequenceUtils::str_to_nuc(con);
        for (auto i = 0; i < con.size(); i++) {
            auto cur_allowed = SequenceUtils::nuc_to_bool(int_pattern[i]);
            cur_con.push_back(cur_allowed);
            int n_allowed = 0;
            int all_ind = 0;

            for (auto j = 0; j < cur_allowed.size(); j++) {
                if (cur_allowed[j]) {
                    varval_to_ids[i][j].push_back(i_val);
                    all_ind = j;
                    n_allowed++;
                }
            }

            if (n_allowed == 0) {
                allowed_ind.back().push_back(-2);
            } else if (n_allowed == 1) {
                allowed_ind.back().push_back(all_ind);
            } else {
                allowed_ind.back().push_back(-1);
            }
        }

        ids_to_allowed.push_back(cur_con);
        ++i_val;
    }
}



bool WordConstraint::propagate_constraint(int modified,
        SolveStack & sstack, const SolveStruc & ss) const {
    bool satisfied = true;
    const vec<trinary> & lookup = ss.value_allowed[supp_var];
    if (modified == supp_var) {
        // Need to map down to nucleotide variables
        for (auto i_var = 0; satisfied && i_var < nuc_ids.size(); i_var++) {
            auto c_var = nuc_ids[i_var];
            for (auto i_val = 0; satisfied && i_val < ss.value_allowed[c_var].size(); i_val++) {
                if (ss.value_allowed[c_var][i_val]) {
                    const vec<int> & tmpvec = varval_to_ids[i_var][i_val];
                    bool cursat = false;
                    for (auto i_supp = 0; !cursat && i_supp < tmpvec.size(); i_supp++) {
                        if (lookup[tmpvec[i_supp]]) cursat = true;
                    }

                    if (!cursat) {
                        if (ss.value_allowed[c_var][i_val] == NUPACK_VV_TRUE) {
                            // NUPACK_DEBUG_PRINT("Failing " << c_var << " " << i_val);
                            satisfied = false;
                        } else {
                            // NUPACK_DEBUG_PRINT("Marking " << c_var << " " << i_val);
                            sstack.push_back(c_var, i_val, NUPACK_VV_FALSE);
                        }
                    }
                }
            }
        }
    } else {
        int i_var = 0, c_var;
        for (i_var = 0; i_var < nuc_ids.size(); i_var++) {
            c_var = nuc_ids[i_var];
            if (c_var == modified) break;
        }
        // Need to map down to supplementary variable
        for (auto i_supp = 0; satisfied && i_supp < ss.value_allowed[supp_var].size(); i_supp++) {
            if (ss.value_allowed[supp_var][i_supp]) {
                bool cur_match = false;
                if (allowed_ind[i_supp][i_var] == -2) {
                    cur_match = false;
                } else if (allowed_ind[i_supp][i_var] >= 0) {
                    cur_match = ss.value_allowed[c_var][allowed_ind[i_supp][i_var]] > 0;
                } else {
                    for (auto i_val = 0; i_val < ids_to_allowed[i_supp][i_var].size()
                            && !cur_match; i_val++) {
                        cur_match = (ids_to_allowed[i_supp][i_var][i_val] &&
                                     ss.value_allowed[c_var][i_val]);
                    }
                }

                if (!cur_match) {
                    if (ss.value_allowed[supp_var][i_supp] == NUPACK_VV_TRUE) {
                        satisfied = false;
                    } else {
                        sstack.push_back(supp_var, i_supp, NUPACK_VV_FALSE);
                    }
                }
            }
        }
    }

    return satisfied;
}



vec<int> WordConstraint::get_constrained_vars() const {
    vec<int> ids = nuc_ids;
    ids.push_back(supp_var);
    return ids;
}



MatchConstraint::MatchConstraint(const vec<int> & vars,
                                 const string & word, vec<double> min_match,
                                 vec<double> max_match) {
    NUPACK_CHECK(min_match.size() == max_match.size(), "min and max match size mismatch");
    NUPACK_CHECK(word.size() == vars.size(), "Word / variable size mismatch");

    try {
        vec<int> nuc_word = SequenceUtils::str_to_nuc(word);

        // add match ranges
        zip(min_match, max_match, [&](auto min, auto max) {
            ranges.emplace_back(min, max);
        });

        // add constrained variables and reference sequence
        zip(nuc_word, vars, [&](auto const & w, auto v) {
            match_nucs.push_back(SequenceUtils::nuc_to_bool(w));
            nuc_ids.push_back(v);
        });

    } catch (NupackException & e) {
        e.print_message(std::cerr);
        NUPACK_DESIGN_ERROR("Error converting nucleotides in match constraint");
    }

    // for (auto & r : ranges) {
    //   NUPACK_LOG_INFO("min: " + to_string(r.first) + ", max: " + to_string(r.second))
    // }
}



bool MatchConstraint::propagate_constraint(int modified,
        SolveStack & sstack, const SolveStruc & ss) const {
    int min_matched = 0;
    int max_matched = 0;

    int n_tot = nuc_ids.size();
    for (auto i = 0; i < n_tot; ++i) {
        auto c_var = nuc_ids[i];
        auto can_match = false;
        auto must_match = true;

        auto cur_size = match_nucs[i].size();
        for (auto c_val = 0; c_val < cur_size; c_val++) {
            if (ss.value_allowed[c_var][c_val] && match_nucs[i][c_val]) {
                can_match = true;
            }
        }

        for (auto c_val = 0; c_val < cur_size; c_val++) {
            if (ss.value_allowed[c_var][c_val] && !match_nucs[i][c_val]) {
                must_match = false;
            }
        }

        if (must_match) min_matched++;
        if (can_match)  max_matched++;
    }

    double min_frac_matched = ((double) min_matched) / n_tot;
    double max_frac_matched = ((double) max_matched) / n_tot;

    double increment = 1.0 / n_tot;

    bool can_add_match = false;
    bool can_sub_match = false;
    bool can_match_tot = false;

    for (auto & r : ranges) {
        if (min_frac_matched < r.second && max_frac_matched > r.first) {
            can_match_tot = true;
            if (min_frac_matched + increment < r.second) can_add_match = true;
            if (max_frac_matched - increment > r.first) can_sub_match = true;
        }
    }

    if (can_match_tot && !can_sub_match) {
        for (auto i = 0; i < n_tot; i++) {
            auto c_var = nuc_ids[i];
            auto cur_match_nuc = match_nucs[i];
            auto can_match = false;
            for (auto c_val = 0; c_val < cur_match_nuc.size(); c_val++) {
                if (ss.value_allowed[c_var][c_val] && cur_match_nuc[c_val]) {
                    can_match = true;
                }
            }

            if (can_match) {
                for (auto c_val = 0; c_val < cur_match_nuc.size(); c_val++) {
                    if (ss.value_allowed[c_var][c_val] && !cur_match_nuc[c_val]) {
                        sstack.push_back(c_var, c_val, NUPACK_VV_FALSE);
                    }
                }
            }
        }
    }

    if (can_match_tot && !can_add_match) {
        for (auto i = 0; i < n_tot; i++) {
            auto c_var = nuc_ids[i];
            auto cur_match_nuc = match_nucs[i];
            auto must_match = true;
            for (auto c_val = 0; c_val < cur_match_nuc.size(); c_val++) {
                if (ss.value_allowed[c_var][c_val] && !cur_match_nuc[c_val]) {
                    must_match = false;
                }
            }

            if (!must_match) {
                for (auto c_val = 0; c_val < cur_match_nuc.size(); c_val++) {
                    if (ss.value_allowed[c_var][c_val] && cur_match_nuc[c_val]) {
                        sstack.push_back(c_var, c_val, NUPACK_VV_FALSE);
                    }
                }
            }
        }
    }

    return can_match_tot;
}



int ConstraintHandler::add_nucleotide_variable(int nuc_con) {
    vec<trinary> cons = SequenceUtils::nuc_to_bool(nuc_con);

    return add_variable(cons);
}



int ConstraintHandler::add_variable(vec<trinary> allowed) {
    for (auto &a: allowed) a = a ? NUPACK_VV_UNSET : NUPACK_VV_FALSE;

    value_allowed.push_back(allowed);
    var_constraint_map.emplace_back();
    return value_allowed.size() - 1;
}



int ConstraintHandler::get_first_allowed(const vec<trinary> & allowed, int i_set) {
    ++i_set; // adapt to generic find_nth interface
    auto index = find_nth(allowed, [](const trinary & a) {return a == NUPACK_VV_UNSET; }, i_set) - allowed.begin();
    if (index < allowed.size()) return index;

    NUPACK_DESIGN_ERROR("Attempting to find allowed variable without a possibility");
}



int ConstraintHandler::get_n_true(const vec<trinary> & allowed) {
    return count(allowed, NUPACK_VV_TRUE);
}



int ConstraintHandler::get_n_unset(const vec<trinary> & allowed) {
    return count(allowed, NUPACK_VV_UNSET);
}



// those that are either NUPACK_VV_TRUE or NUPACK_VV_UNSET
int ConstraintHandler::get_n_allowed(const vec<trinary> & allowed) {
    return allowed.size() - count(allowed, NUPACK_VV_FALSE);
}



vec<int> ConstraintHandler::make_mutation(vec<int> mut_vars, vec<int> start) {
    vec<int> newstart = start;
    vec<int> ret = start;

    NUPACK_CHECK(start.size() == value_allowed.size(), "start (" + to_string(start.size()) + ") / allowed size (" + to_string(value_allowed.size()) + ")mismatch");

    for (auto j = 0; j < ret.size(); j++) {
        if (value_allowed[j].size() == 4) {
            NUPACK_CHECK(ret[j] >= 0 && ret[j] <= 3, "Invalid nucleotide code at start: " + to_string(j) + " : " + to_string(ret[j]));
        }
    }

    for (auto & c_var : mut_vars) {
        NUPACK_CHECK(c_var < value_allowed.size(), "attempting to mutate invalid variable");

        // temporary for holding domains excluding initial values and modified during mutation
        AllowTable allowed = value_allowed;
        int n_allowed = get_n_allowed(allowed[c_var]);

        if (n_allowed > 1) {
            int old_val = newstart[c_var];

            if (newstart[c_var] >= 0) allowed[c_var][newstart[c_var]] = false; // disallow initial value of the variable being mutated
            newstart[c_var] = -1;
            ret = find_closest(newstart, allowed);

            if (ret[0] >= 0) {
                newstart = ret;
            } else {
                newstart[c_var] = old_val;
                ret = newstart;
            }


            // this should only be hit if there is something wrong with find_closest
            for (auto j = 0; j < ret.size(); j++) {
                if (allowed[j].size() == 4) {
                    NUPACK_CHECK(ret[j] >= 0 && ret[j] <= 3,
                                 "Invalid nucleotide code[" + to_string(j) + "]: " + to_string(ret[j]));
                }
            }
        }
    }

    return ret;
}



vec<int> ConstraintHandler::get_possible_nucleotides() const {
    SolveStruc solver;

    solver.value_allowed = value_allowed;
    solver.start.assign(value_allowed.size(), -1);
    solver.weight.assign(solver.value_allowed.size(), 1.0);
    // solver.init_weights(solver.value_allowed);

    auto root = std::make_shared<VariableNode>();
    bool success = propagate_all(root, solver);

    NUPACK_CHECK(success, "No nucleotides found satisfy these constraints");
    return SequenceUtils::bool_to_nuc(solver.value_allowed);
}



bool ConstraintHandler::propagate_all(std::shared_ptr<VariableNode> node,
                                      SolveStruc & ss) const {
    SolveStack sstack;

    bool success = true;
    for (auto const & constraint : constraints) {
        fork(constraint, [&](auto const &c) {
            for (auto const & var : c.get_constrained_vars()) {
                success = c.propagate_constraint(var, sstack, ss) && success;
            }
        });
    }

    success = success && node->add_implications(ss.value_allowed, sstack);
    success = success && propagate(node, ss);

    return success;
}



bool ConstraintHandler::propagate(std::shared_ptr<VariableNode> node,
                                  SolveStruc & ss) const {
    SolveStack sstack(node->get_v());

    bool success = true;
    while (sstack.size() > 0 && success) {
        std::set<int> cur_vars_changed_unique(indirect_view(sstack.v, [](auto const &c) {return c.var;}));
        sstack.clear();

        for (auto c_var : cur_vars_changed_unique) {
            auto n_unset = get_n_unset(ss.value_allowed[c_var]);
            auto n_true = get_n_true(ss.value_allowed[c_var]);

            if (n_true > 1) {
                success = false;
                break;
            } else if (n_true > 0) {
                if (n_unset > 0)
                    for (auto i = 0; i < ss.value_allowed[c_var].size(); i++)
                        if (ss.value_allowed[c_var][i] == NUPACK_VV_UNSET)
                            sstack.push_back(c_var, i, NUPACK_VV_FALSE);
            } else if (n_unset == 0) {
                if (n_true != 1) {
                    success = false;
                    break;
                }
            } else if (n_unset == 1) {
                sstack.push_back(c_var, get_first_allowed(ss.value_allowed[c_var]), NUPACK_VV_TRUE);
            }

            for (auto c_con : var_constraint_map[c_var])
                success = success && fork(constraints[c_con], [&](auto const &c) {
                    return c.propagate_constraint(c_var, sstack, ss);
                });
        }

        if (!success) break;
        success = success && node->add_implications(ss.value_allowed, sstack);
    }

    node->update_cost(ss.start);

    return success;
}



int ConstraintHandler::get_n_possibilities() const {
    int res = 0;
    for_each(value_allowed, [&](auto const & v) {
        res += count(v, NUPACK_VV_UNSET);
    });
    return res;
}



int ConstraintHandler::select_random(const vec<trinary> & allowed) {
    int n_allowed = get_n_allowed(allowed);
    if (n_allowed == 0) {
        return -1;
    } else {
        auto rselect = pick_random_int(0, n_allowed);
        return get_first_allowed(allowed, rselect);
    }
}



vec<int> ConstraintHandler::init_random() const {
    auto vars = find_closest(vec<int>(get_n_variables(), -1), value_allowed);
    NUPACK_CHECK(vars[0] >= 0, "Inconsistent constraints. No viable sequences found");
    return vars;
}



void ConstraintHandler::create_new_branches(std::shared_ptr<VariableNode> parent,
        SolveStruc & solver) const {
    auto pred = [&](auto i) {return get_n_unset(solver.value_allowed[i]) > 0 ? solver.weight[i] : 0;};
    double cur_weight = sum(indices(solver.value_allowed), pred);

    double rval = cur_weight * random_float();
    int last_unset = -1;
    if (rval > 0) last_unset = *(find_cumulative(indices(solver.value_allowed), rval, pred).first);

    if (last_unset >= 0) {
        int min_dist = solver.value_allowed.size();
        int dist_start = parent->get_cost();
        for (auto i = 0; i < solver.value_allowed[last_unset].size(); i++) {
            if (solver.value_allowed[last_unset][i] == NUPACK_VV_UNSET) {
                auto curptr = std::make_shared<VariableNode>(parent);
                SolveStack ss;
                ss.push_back(last_unset, i, NUPACK_VV_TRUE);

                bool success = curptr->add_implications(solver.value_allowed, ss);
                success = success && propagate(curptr, solver);

                if (success) {
                    if (curptr->get_cost() < min_dist) min_dist = curptr->get_cost();
                    if (!solver.min_dist_set || curptr->get_cost() < solver.min_dist)
                        solver.sorter.push(curptr);
                }

                curptr->rollback_variables(solver.value_allowed);
            }
        }
        double decay = 0.5;
        min_dist = std::max(min_dist, dist_start);
        solver.weight[last_unset] = solver.weight[last_unset] * decay
                                    + (1 - decay) * (min_dist - dist_start);
    }
}



vec<int> ConstraintHandler::find_closest(const vec<int> & start,
        const AllowTable & value_allowed_in) const {
    SolveStruc solver;

    auto root = std::make_shared<VariableNode>();

    solver.value_allowed = value_allowed_in;
    solver.start = start;
    solver.weight.assign(solver.value_allowed.size(), 1.0);

    SolveStack stack;
    solver.min_dist = 1e100;
    solver.min_dist_set = false;

    int n_poss = get_n_possibilities();
    auto n_vars = value_allowed.size();
    NUPACK_CHECK(n_vars == (int)solver.value_allowed.size(),
                 "Number of variables is inconsistent " + to_string(n_vars) + " != "
                 + to_string(solver.value_allowed.size()));


     /* initial assignment based purely on incoming domains and their
     differences relative to the initial domains of the solver */
    for (auto i = 0; i < n_vars; i++) {
        NUPACK_CHECK(
            solver.value_allowed[i].size() == value_allowed[i].size(),
            "Invalid value allowed size");

        int n_allowed = 0, n_true = 0;
        for (auto j = 0; j < solver.value_allowed[i].size(); j++) {
            auto sav = solver.value_allowed[i][j];
            auto tav = value_allowed[i][j];

            NUPACK_CHECK(
                tav == NUPACK_VV_UNSET || (tav == NUPACK_VV_TRUE && sav == NUPACK_VV_TRUE)
                || (tav == NUPACK_VV_FALSE && sav == NUPACK_VV_FALSE),
                "Value allowed must be a restriction of base allowed");

            if (sav != tav) stack.push_back(i, j, sav);

            if (sav == NUPACK_VV_UNSET) {n_allowed++;}
            else if (sav == NUPACK_VV_TRUE) {n_true++;}
        }
        NUPACK_CHECK(n_true <= 1, "Only one value allowed per variable");

        if (n_true == 1) {
            for (auto j = 0; j < solver.value_allowed[i].size(); j++) {
                if (solver.value_allowed[i][j] == NUPACK_VV_UNSET &&
                        value_allowed[i][j] == NUPACK_VV_UNSET) {
                    stack.push_back(i, j, NUPACK_VV_FALSE);
                }
            }
        } else if (n_allowed == 1) { // find unset (allowed) value and add true to stack
            auto j = 0;
            for (; j < solver.value_allowed[i].size(); j++) {
                if (solver.value_allowed[i][j] == NUPACK_VV_UNSET
                        && value_allowed[i][j] == NUPACK_VV_UNSET) {
                    break;
                }
            }

            stack.push_back(i, j, NUPACK_VV_TRUE);
        }
    }

    solver.value_allowed = value_allowed;
    // solver.init_weights(solver.value_allowed);

    bool success = root->add_implications(solver.value_allowed, stack);
    success = success && propagate_all(root, solver);

    root->rollback_variables(solver.value_allowed);

    vec<int> res(solver.value_allowed.size(), -1);

    solver.sorter.push(root);
    bool found_sol = false;

    std::shared_ptr<VariableNode> cur_branch = nullptr;
    std::shared_ptr<VariableNode> prev_branch = nullptr;

    // While we haven't assigned all variables and we haven't backtracked
    // to the beginning, guess at a new variable
    while (!found_sol && success && solver.sorter.size() > 0) {
        prev_branch = cur_branch;
        cur_branch = solver.sorter.top();
        solver.sorter.pop();

        VariableNode::change_branch(prev_branch, cur_branch, solver.value_allowed);

        auto n_assigned = cur_branch->get_n_assigned();
        if (n_assigned == n_poss) {
            if (solver.is_better(cur_branch->get_cost())) {
                solver.min_dist_set = true;
                solver.min_dist = cur_branch->get_cost();

                for (auto i = 0; i < solver.value_allowed.size(); i++) {
                    for (auto j = 0; j < solver.value_allowed[i].size(); j++) {
                        NUPACK_CHECK(solver.value_allowed[i][j] != NUPACK_VV_UNSET,
                                     "unset variable " + to_string(i) + ":" + to_string(j) +
                                     "when all " + to_string(n_poss) + " variables are assigned");

                        if (solver.value_allowed[i][j] == NUPACK_VV_TRUE) res[i] = j;
                    }
                }
                found_sol = true;
            }
        }

        if (solver.is_better(cur_branch->get_cost()))
            create_new_branches(cur_branch, solver);
    }

    if (res[0] >= 0) {
        for (auto i = 0; i < res.size(); i++) {
            NUPACK_CHECK(res[i] != -1,
                         "No result found [" + to_string(i) + "]: " + to_string(res[i]));
        }
    }

    return res;
}

}
}


