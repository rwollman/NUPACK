/** \file System.h
 * @brief Contains System, representing all strands present in a State or set of States
 */
#pragma once
#include "../types/Sequence.h"
#include "../types/PairList.h"
#include "../iteration/Transform.h"
#include "../standard/Ptr.h"
#include "../standard/Set.h"
#include "ComplexSet.h"

namespace nupack {

/// The system is stored as a contiguous, concatenated array of bases
/// for all the strands. The strands are kept as views into this sequence
/// This is easier for lookup, but worse for dynamic modification
/// of the system contents. However, the system is so small in memory
/// that it is probably easiest to reserve a large chunk of sequence space
/// so that iterators will not be invalidated
class System : public ConstIndexable<System>, TotallyOrdered {

    void make_strands() {
        strands.clear();
        if (total_sequence.empty()) return;
        for (auto i : range(len(nicks) - 1))
            strands.emplace_back(iterator_at(nicks[i]), iterator_at(nicks[i+1]));
    }

public:
    using StrandIter = SubsequenceList::const_iterator;

    System() {}

    /// Make a System from a container of strands
    System(StrandList const &v);

    /// Make a System from a string or container of strings
    template <class S=vec<string>>
    System(S const &s) : System(to_sequences<StrandList>(s)) {}

    /**************************************************************************************/

    /// Contiguous container for concatenated strands
    Sequence total_sequence;
    /// Views into total_sequence representing each strand
    SubsequenceList strands;
    /// Start/end of each strand, strand of each base
    small_vec<iseq> nicks, strand_map;

    NUPACK_REFLECT(System, total_sequence, strands, nicks, strand_map);

    /**************************************************************************************/

    System(System const &y)  : total_sequence(y.total_sequence), nicks(y.nicks),
        strand_map(y.strand_map) {make_strands();}

    System(System &&y) : total_sequence(std::move(y.total_sequence)), nicks(std::move(y.nicks)),
        strand_map(std::move(y.strand_map)) {make_strands();}

    void swap(System &y) {std::swap(*this, y);}

    System & operator=(System const &y) {members() = members_of(y); make_strands(); return *this;}
    System & operator=(System &&y) {members() = members_of(std::move(y)); make_strands(); return *this;}

    /**************************************************************************************/

    SubsequenceList const & iter() const {return strands;}
    bool operator<(System const &o) const {return this != &o && total_sequence < o.total_sequence;}
    bool operator==(System const &o) const {return this == &o || total_sequence == o.total_sequence;}
    auto hash() const {return hash_of(total_sequence);}

    /// Return index corresponding to iterator in the total sequence
    auto index(BaseIter it) const {return it - total_sequence.data();}
    /// Return starting index of a strand from the strand iterator
    auto begin_of_strand(StrandIter it) const {return at(nicks, it - begin_of(strands));}
    /// Return past-the-end index of a strand from the strand iterator
    auto end_of_strand(StrandIter it) const {return at(nicks, it - begin_of(strands) + 1);}
    /// Return next strand within a loop structure recursion
    template <class V>
    StrandIter next_strand_it(value_type_of<V> j, V const &pairs) const {
        while (j != begin_of_strand(strand_it_of(j))) j = pairs[--j]; return strand_it_of(j);}
    /// Return strand index of a sequence iterator
    auto strand_of(BaseIter it) const {return strand_map[it - total_sequence.data()];}
    /// Return strand iterator of a sequence index
    StrandIter strand_it_of(int loc) const {return next(strands, strand_map[loc]);}
    /// Return whether a position i is the past-the-end index of
    bool is_strand_end(iseq i) const {
        if (i == len(total_sequence)) return true;
        return (total_sequence[i] == Base('_') && total_sequence[i - 1] == Base('_'));
    }
    /// Return iterator in sequence from an index
    BaseIter iterator_at(iseq i) const {return total_sequence.data() + i;}
    BaseIter total_begin() const {return total_sequence.data();}
    BaseIter total_end() const {return total_sequence.data() + total_sequence.size();}
    /// Number of nucleotides in the system
    auto n_bases() const {return len(total_sequence) - 2 * len(strands);}

    auto save_repr() const {return vmap<StrandList>(strands, view);}

    void load_repr(StrandList const &seqs) {if (!seqs.empty()) *this = System(seqs);}

    void reserve(iseq n) {total_sequence.reserve(n);}

};

/******************************************************************************************/

template <class V, class ...Ts>
void build_complex(V &loops, System const &s, PairList const &pairs, vec<System::StrandIter> &strands, Ts const &...ts) {
    std::vector<std::tuple<usize, iseq, iseq>> queue;
    if constexpr(has_capacity<V>) queue.reserve(loops.capacity() - loops.size());
    queue.emplace_back(len(loops), s.begin_of_strand(strands.front()), s.end_of_strand(strands.front()) - 1);
    loops.emplace_back(len(loops), Ether, strands.front()->begin(), ts...);

    while (!queue.empty()) {
        auto &q = queue.back(); // LIFO so that a depth first search is used
        auto const index = std::get<0>(q);

        if (loops[index].next_pair(s, len(loops), std::get<1>(q), std::get<2>(q), pairs, strands)) {
            auto d = std::get<1>(q), e = pairs[d];
            std::get<1>(q) = e;
            NUPACK_DREQUIRE(d, !=, e);
            queue.emplace_back(len(loops), d, e);
            loops.emplace_back(len(loops), index, s.iterator_at(d), ts...);
        } else {
            queue.pop_back();
        }
    }
}

/******************************************************************************************/

/// Recurse through a complex of loops calling a callback on each new loop, takes starting index and returns ending index
// Edge build_complex(vec<StrandIter> &its, Edge, PartialLoop, PairList const &pairs, Callback const &o) const;
/// Recurse through each complex in a PairList and return a ComplexSet
template <class V, class ...Ts>
ComplexSet build_complex_set(V &loops, System const &s, PairList const &pairs, Ts const &...ts) {
    if (len(pairs) != sum(s.strands, len))
        NUPACK_ERROR("number of nucleotides doesn't match length of pair list", len(pairs), sum(s.strands, len));

    ComplexSet out{len(s.strands)};

    auto pool = make_set(iterators(s.strands));
    while(!pool.empty()) {
        auto p = *begin_of(pool); // take the first strand not incorporated already
        vec<System::StrandIter> strands = {p};
        build_complex(loops, s, pairs, strands, ts...);
        // reconstruct the order of strands in the complex
        out.emplace_back(vmap<ComplexSet::Indices>(strands, [&](auto it) {
            pool.erase(it);
            return it - begin_of(s.strands);
        }));
    }
    return out;
}

/******************************************************************************************/

}
