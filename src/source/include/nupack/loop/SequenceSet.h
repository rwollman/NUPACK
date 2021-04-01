#pragma once

#include "../algorithms/Utility.h"
#include "../iteration/Transform.h"
#include "../types/Sequence.h"
#include "../model/Move.h"

namespace nupack {

/******************************************************************************************/

class SequenceSet : TotallyOrdered, public Indexable<SequenceSet> {

protected:
    /// Sequences between base pairs, inclusive
    SubsequenceList seqs;
    /// Index of sequence that nick is before, or NoNick if no nick
    Nick n = NoNick;

public:

    NUPACK_REFLECT(SequenceSet, seqs, n);
    auto & iter() const {return seqs;}
    auto & vec() const {return seqs;}

    void append(Subsequence s) {seqs.emplace_back(s);}
    void set_last(BaseIter b) {seqs.back().set_end(b);}

    auto finalize() {
        auto shift = rotate_min_begin(seqs);
        n = find_nick(seqs);
        return shift;
    }

    /// Default constructor, generally only to be used for container convenience
    SequenceSet() = default;

    SequenceSet(BaseIter b) : seqs{b} {}

    SequenceSet(SubsequenceList v) : seqs(std::move(v)) {finalize();}

    /************************************************************************************/

    /// Is there a strand break?
    bool exterior() const {return n != NoNick;}
    /// Nick getter
    auto nick() const {return n;}

    /************************************************************************************/

    /// Sequences, separated by ','
    string sequence_string(string const sep=", ") const {
        auto ret = sum(seqs, [&](auto const &s){return make_string(s) + sep;});
        ret.pop_back(); return ret;
    }
    /// Classified for loop (hairpin, etc)
    string name() const {
        if (exterior()) return "Exterior";
        else if (seqs.size() == 1) return "Hairpin";
        else if (seqs.size() == 2) return "Interior";
        else return "Multiple";
    }

    /************************************************************************************/

    auto strand_begin() const {NUPACK_ASSERT(exterior()); return seqs[n].begin();}

    /************************************************************************************/

    bool operator==(SequenceSet const &s) const {return seqs == s.seqs;}
    bool operator<(SequenceSet const &s) const {return seqs < s.seqs;}

    /**************************************************************************************/

    std::pair<iseq, iseq> associate(SequenceSet & k, iseq ps, iseq ks, BaseIter pb, BaseIter kb);
    std::pair<iseq, iseq> dissociate(SequenceSet &k, iseq pk, iseq kp);
    std::pair<iseq, iseq> split(BasePairAddition const &m, SequenceSet & d);

    iseq merge(SequenceSet const &k, iseq pk, iseq kp) {
        seqs = merged_seqs(seqs, k.seqs, pk, kp);
        iseq ret = rotate_min_begin(seqs);
        n = find_nick(seqs);
        if (ret) return ret - 1; else return seqs.size() - 1;
    }
};

/******************************************************************************************/

}
