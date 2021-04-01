#pragma once
#include "../types/Sequence.h"
#include "../standard/Optional.h"
#include "../standard/Variant.h"
#include "custom_csp/constraint_handler.h"
#include <gecode/int.hh>
#include <gecode/set.hh>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>

namespace nupack { namespace newdesign {

using custom_csp::ConstraintHandler;

using Gecode::Space;
using Gecode::IntVarArray;
using Gecode::IntVar;
using Gecode::IntSet;
using Gecode::IntAFC;

vec<int> nuc_values(std::array<bool, 4> const &);
int random_nuc(IntVar);
IntSet nuc_values(Base base);

string multiply_substrings(vec<std::pair<string, int>> const & condensed);

/** to hold pointer to container holding words and index in case of reallocation */
using DictWord = vec<vec<int>>;
using DictWords = vec<DictWord>;
using Dictionary = vec<DictWords>;


struct WordRef {
    Dictionary const * dictionary;
    uint index;

    WordRef() = default;

    DictWords const & words() const { return (*dictionary).at(index); }
    DictWord const & word() const { return words()[0]; }
};


struct NucSpace : public Space {
    using simple_type = True;
    IntVarArray nucs;
    Sequence const * ref;
    IntVarArray extras;

    NucSpace(Sequence const &);
    NucSpace(NucSpace &);

    NucSpace * copy();
    NucSpace * cast_clone();

    void force(int, int);
    void disallow(int, int);

    /* branchers */
    void default_brancher();
    void reference_brancher();
    void cheap_reference_brancher();

    void add_reference(Sequence const &);

    void match_constraint(int, int);
    void complementarity_constraint(int, int, bool wobble);
    void pattern_constraint(vec<int> const &, Sequence const &);
    void diversity_constraint(vec<int> const &, int, int);
    // void word_constraint(vec<int> const &, vec<Sequence> const &); /* library and window */
    void word_constraint(vec<int> const &, WordRef);
    void similarity_constraint(vec<int> const &, WordRef, std::pair<int, int>);

    explicit operator Sequence() const;

    friend std::ostream & operator<<(std::ostream &os, NucSpace const &space) {
        return os << "nucs: " << space.nucs << ", ref: " << space.ref << ", extras: " << space.extras;
    }
};

struct RunningAverage {
    uint count {0};
    real average {0};

    RunningAverage() = default;

    real add_value(real);
};

static_assert(std::is_copy_constructible_v<ConstraintHandler>);

struct Constraints {
    std::unique_ptr<NucSpace> initial {nullptr};
    ConstraintHandler handler;
    /* max time to allow new constraint implementation to search */
    int msec_cutoff {1000};
    /* once old mutation has been used, tie cutoff to time this takes to bound
    extra computation */
    RunningAverage old_mut_time;

    int num_extra_vars {0};
    /** reference sequences for word constraints */
    std::unique_ptr<Dictionary> dictionary {std::make_unique<Dictionary>()};

    Constraints(Sequence const &);

    Constraints() = default;
    Constraints(Constraints const &c) = delete; //: handler(c.handler), msec_cutoff(c.msec_cutoff), old_mut_time(c.old_mut_time), num_extra_vars(c.num_extra_vars)
    Constraints & operator=(Constraints const &) = delete;
    Constraints(Constraints &&) = default;
    Constraints & operator=(Constraints &&) = default;

    /* functions for adding constraints to NucSpace */
    void match_constraint(int, int);
    void complementarity_constraint(int, int, bool wobble);
    void pattern_constraint(vec<int> const &, Sequence const &);
    void diversity_constraint(vec<int> const &, int, int);
    void word_constraint(vec<int> const &, vec<Sequence> const &); /* library and window */
    void similarity_constraint(vec<int> const &, Sequence const &, std::pair<real, real>);

    Optional<Sequence> initial_sequence();
    Optional<Sequence> make_mutation(Sequence const &, vec<int>);

    int sequence_length() const;

    Gecode::Search::Options search_options() const;

    NUPACK_REFLECT(Constraints, initial);

private:
    Variant<Sequence, bool> make_new_mutation(Sequence const &, int);
    Optional<Sequence> make_old_mutation(Sequence const &, int);
    Variant<Sequence, bool> new_initial_sequence();
    Optional<Sequence> old_initial_sequence();

    void update_cutoff(real);
};

}
}
