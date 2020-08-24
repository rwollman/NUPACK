#pragma once

#include "named_spec.h"
#include "design_debug.h"

#include "types.h"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <queue>
#include <tuple>
#include <iostream>
#include <typeinfo>

#define DISTANCE_COST_NOTEQUAL 1

namespace nupack {
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
    SolveStack(const std::vector<VariableTuple> & in) : v(in) {}
    
    void clear() { v.clear(); }
    void push_back(int variable, int value, trinary allowed_val) {
      v.emplace_back(variable, value, allowed_val);
    }
    size_t size() const { return v.size(); }

    std::vector<VariableTuple> v;
  };

  class VariableNode {
    public:
      VariableNode();
      VariableNode(std::shared_ptr<VariableNode> parent);
      ~VariableNode() { genrand_real1(); } // backwards compatibility for regression testing

      /*
       * @allowed: the current allowed matrix
       * @variables: variables to add
       * @values: values to set
       * @allowed: value to set it to
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
      void update_cost(const std::vector<int> & orig);

      /*
       * compute the total number of variables assigned
       */
      int get_n_assigned() const { return n_assigned; }
      /*
       * get the parent of the current node; only really used in member functions.
       */
      std::shared_ptr<VariableNode> get_parent() const { return parent; }

      const std::vector<VariableTuple> & get_v() const { return v; }

      // no parent
      bool is_root() const { return !(parent); }

      double get_tiebreaker() const { return tiebreaker; }

    struct VariableNodeComp {
        using var = std::shared_ptr<VariableNode>;
        // lexicographic comparison between VariableNodes in depth, negative cost, and tiebreaker
        bool operator() (const var a, const var b) const {
          return std::make_tuple(a->get_depth(), -(a->get_cost()), a->get_tiebreaker()) <
              std::make_tuple(b->get_depth(), -(b->get_cost()), b->get_tiebreaker());
        }
    };
    
    private:
      void init_variables();

      int depth;
      int cost;
      double tiebreaker;
      int n_assigned;

      std::shared_ptr<VariableNode> parent;

      std::vector<VariableTuple> v;
  };


  class SolveStruc {
    using item = std::shared_ptr<VariableNode>;
    using CVarQueue = std::priority_queue<item, std::deque<item>, VariableNode::VariableNodeComp>;

    public:
      std::vector<int> start;
      AllowTable value_allowed;
      std::vector<double> weight;
      CVarQueue sorter;

      double min_dist;
      bool min_dist_set;
  };

  class Constraint {
    public:
      /* 
       * Propagate based on the variable 'modified', append implied changed
       * variables into ss
       *
       * to denote setting variable i to j, append i to vars, j to vals, and set the
       * value in value_allowed
       *
       * Return true if propagation succeeded and this constraint has a consistent
       * solution, false if there doesn't exist a consistent solution with this 
       * variable state.
       */
      virtual bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const = 0;

      /*
       * Get the variables that can change the status of the constraint
       *
       * Any variable that can change the satisfiability of the constraint is
       * handled in here.
       */
      virtual std::vector<int> get_constrained_vars() const = 0;

      virtual std::shared_ptr<Constraint> clone() const = 0;

      virtual ~Constraint() = default;
  };

  template <typename ConstraintType> 
  class Constraint_CRTP : public Constraint {
    public: 
      virtual std::shared_ptr<Constraint> clone() const {
        return std::make_shared<ConstraintType>(static_cast<ConstraintType const &>(*this));
      };
  };

  typedef enum COMPLEMENT_STRENGTH {
    NUPACK_CS_NONE,
    NUPACK_CS_WEAK,
    NUPACK_CS_STRONG
  } complement_strength_t;

  class CompConstraint : public Constraint_CRTP<CompConstraint> {
    public:
      CompConstraint(int i, int j, complement_strength_t strength) : 
          i(i), j(j), strength(strength) {}
      ~CompConstraint() {}

      bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;
      
      std::vector<int> get_constrained_vars() const { return {i, j}; }

    private:
      int i;
      int j;
      complement_strength_t strength;
  };

  class IdentConstraint : public Constraint_CRTP<IdentConstraint> {
    public:
      IdentConstraint(int i, int j) : i(i), j(j) {}
      ~IdentConstraint() {}

      bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

      std::vector<int> get_constrained_vars() const { return {i, j}; }

    private:
      int i;
      int j;
  };

  class PatternConstraint : public Constraint_CRTP<PatternConstraint> {
    public:
      /*
       * name: the name of the domain or strand
       * constraint: the sequence string representing the prevented pattern
       * spec: the resolved sequence specification
       */
      PatternConstraint(const std::string & name, const std::string & constraint,
          const SequenceSpec & spec, const std::vector<int> & poss_nucs);
      ~PatternConstraint() {}

      /* * propagate the constraints for the current pattern constraint */
      bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

      /* get the variables modified by this constraint */
      std::vector<int> get_constrained_vars() const { return this->nuc_ids; }

    private:
      std::string name;
      std::string constraint;
      
      // Nucleotide ids in order they are prevented in
      std::vector<int> nuc_ids;
      std::vector<trinary> starts;

      // map from nucleotide ids to position in nuc_ids
      std::map<int, int> nuc_id_map;
      AllowTable pattern;
  };

  class WordConstraint : public Constraint_CRTP<WordConstraint> {
    public:
      // Create the constraint and add auxiliary variable to the 
      // constraint_hander
      WordConstraint(const std::vector<int> & vars,
          const std::vector<std::string> & words,
          int additional_var);
      
      bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;
      std::vector<int> get_constrained_vars() const;

      ~WordConstraint() {}

    private:
      void clear();

      int supp_var;
      std::vector<int> nuc_ids;

      std::vector<std::vector<int> > allowed_ind;
      std::vector<std::vector<std::vector<int> > > varval_to_ids;
      std::vector<AllowTable > ids_to_allowed;
  };

  class MatchConstraint;

  struct Match {
    Match(std::string name, std::vector<double> mins, std::vector<double> maxes) : 
        name(name), mins(mins), maxes(maxes) {}

    MatchConstraint make_match_constraint(const SequenceSpec & seqs, const std::string & sequence);

    std::string name;
    std::vector<double> mins;
    std::vector<double> maxes;
  };

  class MatchSpec : public NamedSpec {
    public:
      MatchSpec(const std::string & name, const std::string & seq) :
          NamedSpec(name, "match_spec"), seq(seq) {}
      ~MatchSpec() {}

      void add_match(const std::string & dom_name, const std::vector<double> & mins,
          const std::vector<double> & maxes);

      std::vector<MatchConstraint> make_match_constraints(const SequenceSpec & seqs);

    private:
      std::string seq;
      std::vector<Match> matches;
  };

  class MatchConstraint : public Constraint_CRTP<MatchConstraint> {
    public:
      MatchConstraint(const std::vector<int> & vars, const std::string & words,
          std::vector<double> min_match, std::vector<double> max_match);

      bool propagate_constraint(int modified, SolveStack & sstack, const SolveStruc & ss) const;

      std::vector<int> get_constrained_vars() const;
      ~MatchConstraint() {}

    private:
      // first is min, second is max
      std::vector<std::pair<double, double>> ranges;
      std::vector<int> nuc_ids;
      AllowTable match_nucs;
  };

inline int pick_random_int(int from, int to) { return ((int) (genrand_real1() * (to - from))) + from; }
  
  class ConstraintHandler {
    public:
      ConstraintHandler() {};
      ~ConstraintHandler() {};

      int add_variable(std::vector<trinary> allowed_vals);
      int add_nucleotide_variable(int constraint);
      void add_constraint(const Constraint & con);
      std::vector<int> init_random() const;


      std::vector<int> find_closest(const std::vector<int> & start, 
          const AllowTable & value_allowed) const;

      int get_n_variables() const { return this->value_allowed.size(); }


      std::vector<int> make_mutation(std::vector<int> mut_vars, std::vector<int> start);

      void create_new_branches(std::shared_ptr<VariableNode> parent,
          SolveStruc & solver) const;

      static int get_distance(const std::vector<int> & start, const std::vector<int> & cur);
      static int get_n_allowed(const std::vector<trinary> & allowed);
      static int get_n_unset(const std::vector<trinary> & allowed);
      static int get_n_true(const std::vector<trinary> & allowed);
      static int get_first_allowed(const std::vector<trinary> & allowed, int i_set = 0);
      static int select_random(const std::vector<trinary> & allowed);

      /*
       * @ret: ret[i] = -1 <=> variable i has more than 1 possible value
       *    ret[i] = j variable i must be value j
       *    ret = [] <=> no valid assignments can be made
       */
      std::vector<int> get_possible_nucleotides() const; 


    private:
      int get_n_possibilities() const;
      bool propagate(std::shared_ptr<VariableNode> cur, SolveStruc & ss) const;
      bool propagate_all(std::shared_ptr<VariableNode> node, SolveStruc & ss) const;
      
      std::vector<std::shared_ptr<Constraint>> constraints;
      AllowTable value_allowed;

      // Map from variables to constraints that they take part in 
      // (constraints to check when the range of the variable changes)
      std::vector<std::vector<int> > var_constraint_map;
  };

};

