#pragma once

#include "named_spec.h"
#include "design_debug.h"

#include "types.h"
#include "adapter.h"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <queue>
#include <tuple>
#include <iostream>
#include <typeinfo>

#include "../../reflect/Reflection.h"

#define DISTANCE_COST_NOTEQUAL 1

namespace nupack {
namespace custom_csp {

typedef enum VAR_VALUE {
    NUPACK_VV_FALSE = 0,
    NUPACK_VV_TRUE = 1,
    NUPACK_VV_UNSET = 2,
} trinary_t;

inline void print_table(AllowTable & allow_table, std::ostream & out) {
    for (auto & at : allow_table) for (auto & el : at) out << el;
    out << std::endl;
}

struct VariableTuple {
    VariableTuple() : VariableTuple(-1, -1, NUPACK_VV_UNSET) {};
    VariableTuple(int var, int val, trinary trit) : var(var), val(val), trit(trit) {}

    int var;
    int val;
    trinary trit;
};

struct SolveStack {
    SolveStack() = default;
    SolveStack(const vec<VariableTuple> & in) : v(in) {}

    void clear() { v.clear(); }
    void push_back(int variable, int value, trinary allowed_val) {
        v.emplace_back(variable, value, allowed_val);
    }
    size_t size() const { return v.size(); }

    vec<VariableTuple> v;
};

class VariableNode {
    int depth;
    int cost;
    double tiebreaker;
    int n_assigned;

    std::shared_ptr<VariableNode> parent;

    vec<VariableTuple> v;

    void init_variables();

public:
    VariableNode();
    VariableNode(std::shared_ptr<VariableNode> parent);
    ~VariableNode() { random_float(); } // backwards compatibility for regression testing

    /*
     * @allowed: the current allowed matrix
     * @stack
     *      @variables: variables to add
     *      @values: values to set
     *      @allowed: value to set it to
     */
    bool add_implications(AllowTable & allow_table, const SolveStack & stack);

    /*
     * @allowed: the current allowed matrix, will be rolled back from
     * the current branch (this) and then assigned to according to the
     * new branch (other)
     *
     * Rollback continues until the most recent common ancestor is found
     *
     * Error thrown if allowed is inconsistent with either branch.
     *
     */
    static void change_branch(std::shared_ptr<VariableNode> from, std::shared_ptr<VariableNode> to,
                              AllowTable & allow_table);

    /*
     * @allow_table: the current allowed matrix, will be cleared based
     *  on the current node's implications.
     *
     * Throws an exception if the variables aren't assigned according to
     * the current node
     */
    void rollback_variables(AllowTable & allow_table);

    /*
     * @allow_table: the current allowed matrix, will be assigned to based
     *  on the current node's implications
     *
     * Throws an exception if the variables assigned in the current node are
     * already set.
     */
    void assign_variables(AllowTable & allow_table);


    /* get the number of children for the node */
    int get_depth() const { return depth; }
    /* get the total cost of the current node */
    int get_cost() const { return cost; }
    /* compute the cost based on the vector orig */
    void update_cost(const vec<int> & orig);

    /*
     * compute the total number of variables assigned
     */
    int get_n_assigned() const { return n_assigned; }
    /*
     * get the parent of the current node; only really used in member functions.
     */
    std::shared_ptr<VariableNode> get_parent() const { return parent; }

    const vec<VariableTuple> & get_v() const { return v; }

    // no parent
    bool is_root() const { return !(parent); }

    double get_tiebreaker() const { return tiebreaker; }

    struct VariableNodeComp {
        using var = std::shared_ptr<VariableNode>;
        // lexicographic comparison between VariableNodes in depth, negative cost, and tiebreaker
        bool operator()(const var a, const var b) const {
            return std::make_tuple(a->get_depth(), -(a->get_cost()), a->get_tiebreaker()) <
                   std::make_tuple(b->get_depth(), -(b->get_cost()), b->get_tiebreaker());
        }
    };
};


struct SolveStruc {
    using item = std::shared_ptr<VariableNode>;
    using CVarQueue = std::priority_queue<item, std::deque<item>, VariableNode::VariableNodeComp>;

    vec<int> start;
    AllowTable value_allowed;
    vec<double> weight;
    CVarQueue sorter;

    double min_dist;
    bool min_dist_set;

    /* set weights to the number of unset values in the table (the domain size) */
    void init_weights(AllowTable const &table);
    bool is_better(double x) const {return !min_dist_set || x < min_dist;}
};


typedef enum COMPLEMENT_STRENGTH {
    NUPACK_CS_NONE,
    NUPACK_CS_WEAK,
    NUPACK_CS_STRONG
} complement_strength_t;


class CompConstraint : MemberOrdered {
    int i;
    int j;
    complement_strength_t strength;

public:
    NUPACK_REFLECT(CompConstraint, i, j, strength);
    CompConstraint() = default;
    CompConstraint(int i, int j, complement_strength_t strength) :
        i(i), j(j), strength(strength) {}

    bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

    vec<int> get_constrained_vars() const { return {i, j}; }

    friend std::ostream & operator<<(std::ostream &os, CompConstraint const &c) {
        return os << "i: " << c.i << ", j: " << c.j;
    }
};


class IdentConstraint : MemberOrdered {
    int i;
    int j;

public:
    NUPACK_REFLECT(IdentConstraint, i, j);
    IdentConstraint() = default;
    IdentConstraint(int i, int j) : i(i), j(j) {}

    bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

    vec<int> get_constrained_vars() const { return {i, j}; }
};


class PatternConstraint : MemberOrdered {
    string constraint;

    // Nucleotide ids in order they are prevented in
    vec<int> nuc_ids;
    vec<trinary> starts;

    // map from nucleotide ids to position in nuc_ids
    std::map<int, int> nuc_id_map;
    AllowTable pattern;

public:
    NUPACK_REFLECT(PatternConstraint, constraint, nuc_ids, starts, nuc_id_map, pattern);
    PatternConstraint() = default;
    /**
     * @param vars the constrained variables
     * @param constraint the sequence string representing the prevented pattern
     * @param spec the resolved sequence specification
     */
    PatternConstraint(vec<int> const &vars, const string & constraint, const vec<int> & poss_nucs);

    /* * propagate the constraints for the current pattern constraint */
    bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

    /* get the variables modified by this constraint */
    vec<int> get_constrained_vars() const { return this->nuc_ids; }
};


class WordConstraint : MemberOrdered {
    int supp_var;
    vec<int> nuc_ids;

    vec<vec<int> > allowed_ind;
    vec<vec<vec<int> > > varval_to_ids;
    vec<AllowTable > ids_to_allowed;

    void clear();

public:
    NUPACK_REFLECT(WordConstraint, supp_var, nuc_ids, allowed_ind, varval_to_ids, ids_to_allowed);
    WordConstraint() = default;
    // Create the constraint and add auxiliary variable to the
    // constraint_hander
    WordConstraint(const vec<int> & vars,
                   const vec<string> & words,
                   int additional_var);

    bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;
    vec<int> get_constrained_vars() const;

};


class MatchConstraint;

class MatchConstraint : MemberOrdered {
    // first is min, second is max
    vec<std::pair<double, double>> ranges;
    vec<int> nuc_ids;
    AllowTable match_nucs;

public:
    NUPACK_REFLECT(MatchConstraint, ranges, nuc_ids, match_nucs);
    MatchConstraint() = default;
    MatchConstraint(const vec<int> & vars, const string & words,
                    vec<double> min_match, vec<double> max_match);

    bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

    vec<int> get_constrained_vars() const {return nuc_ids;}
};


inline int pick_random_int(int from, int to) { return ((int)(random_float() * (to - from))) + from; }

using constraint_variant = Variant<CompConstraint, IdentConstraint, PatternConstraint, WordConstraint, MatchConstraint>;

class ConstraintHandler : MemberOrdered {
    vec<constraint_variant> constraints;
    AllowTable value_allowed;
    /** @brief Map from variables to constraints that they take part in (constraints
        to check when the range of the variable changes) */
    vec<vec<int> > var_constraint_map;

    int get_n_possibilities() const;
    bool propagate(std::shared_ptr<VariableNode> cur, SolveStruc & ss) const;
    bool propagate_all(std::shared_ptr<VariableNode> node, SolveStruc & ss) const;

public:
    NUPACK_REFLECT(ConstraintHandler, constraints, value_allowed, var_constraint_map);
    ConstraintHandler() = default;

    int add_variable(vec<trinary> allowed_vals);
    int add_nucleotide_variable(int constraint);

    vec<constraint_variant> const & get_constraints() const {return constraints;}

    template <class C, NUPACK_IF(can_convert<C, constraint_variant>)>
    void add_constraint(C const &con) {
        int c_con = constraints.size();

        constraints.push_back(con);
        for (auto c_var : con.get_constrained_vars()) {
            NUPACK_CHECK(c_var < value_allowed.size(),
                         to_string(c_var) + " is not in the current variable set");
            var_constraint_map[c_var].push_back(c_con);
        }
    }

    vec<int> init_random() const;


    vec<int> find_closest(const vec<int> & start, const AllowTable & value_allowed) const;

    int get_n_variables() const { return this->value_allowed.size(); }


    vec<int> make_mutation(vec<int> mut_vars, vec<int> start);

    void create_new_branches(std::shared_ptr<VariableNode> parent,
                             SolveStruc & solver) const;

    static int get_n_allowed(const vec<trinary> & allowed);
    static int get_n_unset(const vec<trinary> & allowed);
    static int get_n_true(const vec<trinary> & allowed);
    static int get_first_allowed(const vec<trinary> & allowed, int i_set = 0);
    static int select_random(const vec<trinary> & allowed);

    /*
     * @ret: ret[i] = -1 <=> variable i has more than 1 possible value
     *    ret[i] = j variable i must be value j
     *    ret = [] <=> no valid assignments can be made
     */
    vec<int> get_possible_nucleotides() const;

    // friend std::ostream & operator<<(std::ostream &os, ConstraintHandler const &) {return os << "ConstraintHandler";}
};



}
}

