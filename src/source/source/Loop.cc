#include <nupack/loop/SequenceSet.h>
#include <nupack/loop/StaticLoop.h>
#include <nupack/model/StackEnumeration.h>

namespace nupack {

/******************************************************************************************/

std::pair<iseq, iseq> SequenceSet::associate(SequenceSet & k, iseq ps, iseq ks, BaseIter pb, BaseIter kb) {
    auto &p = *this;
    SubsequenceList pseqs = {p.seqs[ps]};
    front(pseqs).set_begin(pb);
    circular_cat(pseqs, p.vec(), next(p, ps + 1), next(p, p.nick()));
    circular_cat(pseqs, k.vec(), next(k, k.nick()), next(k, ks));
    pseqs.emplace_back(k.vec()[ks]);
    back(pseqs).set_end(kb + 1);

    SubsequenceList kseqs = {k.vec()[ks]};
    front(kseqs).set_begin(kb);
    circular_cat(kseqs, k.vec(), next(k, ks + 1), next(k, k.nick()));
    circular_cat(kseqs, p.vec(), next(p, p.nick()), next(p, ps));
    kseqs.emplace_back(p.seqs[ps]);
    back(kseqs).set_end(pb + 1);

    p.seqs = std::move(pseqs); k.seqs = std::move(kseqs);

    auto shift1 = rotate_min_begin(p.seqs), shift2 = rotate_min_begin(k.seqs);

    p.n = find_nick(p.seqs); k.n = find_nick(k.seqs);

    return std::make_pair(shift1, shift2);
}

/******************************************************************************************/

std::pair<iseq, iseq> SequenceSet::dissociate(SequenceSet &k, iseq pk, iseq kp) {
    std::tie(seqs, k.seqs) = get_split_seqs(seqs, k.seqs, static_cast<iseq>(n),
                                            static_cast<iseq>(k.n), pk, kp);
    std::pair<iseq, iseq> ret = {rotate_min_begin(seqs), rotate_min_begin(k.seqs)};
    n = find_nick(seqs); k.n = find_nick(k.seqs);
    return ret;
}

/******************************************************************************************/

std::pair<iseq, iseq> SequenceSet::split(BasePairAddition const &m, SequenceSet & d) {
    std::tie(seqs, d.seqs) = split_midway(seqs, next(seqs, m.s1), next(seqs, m.s2), m.b1, m.b2);
    auto shift1 = rotate_min_begin(seqs);
    auto shift2 = rotate_min_begin(d.seqs);
    n = find_nick(seqs); d.n = find_nick(d.seqs);
    return std::make_pair(shift1, shift2);
}

/******************************************************************************************/

string loop_stack_string(small_vec<LoopStackingState::Stack> const &v) {
    string out;
    out.reserve(v.size());
    for (auto i : v) out.push_back(loop_stack_letter(i));
    return out;
}

/******************************************************************************************/

string loop_stack_sequence_string(small_vec<LoopStackingState::Stack> const &v) {
    string s;
    s.resize(v.size());
    izip(s, [&](auto i, char &c) {
        bool left, right;
        switch (v[i]) {
            case LoopStackingState::None:         {left = false; break;}
            case LoopStackingState::LeftDangle:   {left = false; break;}
            case LoopStackingState::RightDangle:  {left = true;  break;}
            case LoopStackingState::LeftStack:    {left = false; break;}
            case LoopStackingState::RightStack:   {c = 's'; return;}
            case LoopStackingState::BothDangle:   {left = true;  break;}
            case LoopStackingState::Disabled:     {left = false; break;}
        }
        switch (v[(i+1) % v.size()]) {
            case LoopStackingState::None:         {right = false; break;}
            case LoopStackingState::LeftDangle:   {right = true;  break;}
            case LoopStackingState::RightDangle:  {right = false; break;}
            case LoopStackingState::LeftStack:    {c = 's'; return;}
            case LoopStackingState::RightStack:   {right = false; break;}
            case LoopStackingState::BothDangle:   {right = true;  break;}
            case LoopStackingState::Disabled:     {right = false; break;}
        }
        if (!left && !right) c = 'n';
        if (!left &&  right) c = '3';
        if ( left && !right) c = '5';
        if ( left &&  right) c = 'b';
    });
    return s;
}

/******************************************************************************************/

}
