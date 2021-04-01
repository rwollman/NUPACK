#pragma once
#include <nupack/types/Sequence.h>

namespace nupack {

/******************************************************************************************/

enum class MoveType {Add, Del, Join};

/******************************************************************************************/

struct BasePairAddition : TotallyOrdered {
    BasePairAddition() = default;
    BasePairAddition(BaseIter b1_, BaseIter b2_, iseq s1_, iseq s2_, real dE_, real rate_=0):
        s1(s1_), s2(s2_), b1(b1_), b2(b2_), dE(dE_), rate(rate_) {
        NUPACK_ASSERT(s1 <= s2);
        if (s1 == s2) NUPACK_ASSERT(b1 <= b2);
    };

    iseq s1, s2;
    BaseIter b1, b2;
    real dE, rate;

    NUPACK_REFLECT(BasePairAddition, s1, s2, b1, b2, dE, rate);

    friend bool operator< (BasePairAddition const &m1, BasePairAddition const &m2) {
        return std::make_tuple(m1.s1, m1.s2, m1.b1, m1.b2)
             < std::make_tuple(m2.s1, m2.s2, m2.b1, m2.b2);
    }

    friend bool operator== (BasePairAddition const &m1, BasePairAddition const &m2) {
        return std::make_tuple(m1.s1, m1.s2, m1.b1, m1.b2)
            == std::make_tuple(m2.s1, m2.s2, m2.b1, m2.b2);
    }
};

/******************************************************************************************/

struct JoinMove {
    usize o1, o2; // loop 1, loop 2
    Base b1, b2;  // base 1, base 2
    real margin, scale; // cumulative information for choosing a rate

    NUPACK_REFLECT(JoinMove, o1, o2, b1, b2, margin, scale);
};

/******************************************************************************************/

/// Represents a base pair addition between two bases in two different complexes
/// Does not include information on which loops the bases are in
struct JoinLoc {
    iseq s;         // index of the strand this base is on (I think)
    BaseIter b;     // iterator to the base (dereference to get the base)
    real dE, hrate; // free energy change contributed by this base, partial rate constant contributed by this base
    JoinLoc() {};
    JoinLoc(iseq s_, BaseIter b_, real dE_, real h) : s(s_), b(b_), dE(dE_), hrate(h) {};

    NUPACK_REFLECT(JoinLoc, s, b, dE, hrate);
};

/******************************************************************************************/

/// Represents a base pair addition between two bases in two different complexes
struct ComplexJoinMove {
    usize loop1, loop2; // index of first loop in its complex, second loop in its complex
    JoinLoc loc1, loc2; // base location information for the two joining bases

    NUPACK_REFLECT(ComplexJoinMove, loop1, loop2, loc1, loc2);

    /// Free energy of the complex join
    /// It assumes the energy is decomposable into left and right pieces
    /// This assumption is hard to change
    auto dE() const {return loc1.dE + loc2.dE;}

    /// Base identities of the join (i.e. A, T)
    auto bases() const {return std::make_pair(*loc1.b, *loc2.b);}

    /// Rate constant of this join, it assumes that the rates are decomposable into a product of left and right
    /// It wouldn't be hard to change this so this isn't assumed below though
    template <class RF>
    auto rate(RF const &rf) const {
        // I have no recollection why there's a 2 but it seems to be right
        // Otherwise enumerating (G*10)+(C*10) gives total rate 50 at beta=0
        // Removed the 2, I think it was coincidence that first two factors
        // below multiplied to 0.5
        return 0.5 * loc1.hrate * loc2.hrate * rf.molarity * rf.bimolecular_scaling;
    }
};

/******************************************************************************************/

struct BasePairDeletion {
    real dE = 0;
    real rate = 0;

    NUPACK_REFLECT(BasePairDeletion, dE, rate);
};

/******************************************************************************************/

/// Small container to hold info for last move a state took
struct StateMove {
    Base_Pair bp; /// Last base pair formed or deleted
    real rate; /// Rate of the last step
    NUPACK_REFLECT(StateMove, bp, rate);
};

/******************************************************************************************/

}
