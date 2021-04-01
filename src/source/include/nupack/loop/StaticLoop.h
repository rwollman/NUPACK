/** \file Loop.h
 * Contains loop class declarations, then includes definitions from Loop_Impl.h
 */

#pragma once
#include "../Forward.h"
#include "SequenceSet.h"
#include "EdgeSet.h"
#include "../state/System.h"

namespace nupack {

struct Move;

/******************************************************************************************/

/**
 * @brief Represents a nucleic acid loop, its possible moves, and its connections to other loops
 * @details Responsibilities are divided between SequenceSet and EdgeSet. SequenceSet holds the sequences
 * along the loop, independent of the loop's neighbors. EdgeSet holds this loop's and its neighbors'
 * indices only.
 *
 * @tparam Energy=Energy Energy type
 * @tparam Joiner=ProductJoiner Controls calculation of complex join moves between exterior loops
 * @tparam Move_Gen=FullMoveGen<Energy> Controls calculation of intra-loop base pair addition moves
 */
template <class SS>
struct StaticLoop : Iterable<StaticLoop<SS>>, TotallyOrdered {

    /************************************************************************************/

    SS seqs;
    EdgeSet edges;

    NUPACK_REFLECT(StaticLoop, seqs, edges);

    StaticLoop() = default;

    template <class T, NUPACK_IF(can_construct<SS, T const &>)>
    StaticLoop(T const &t) : seqs(t) {}

    template <class ...Ts>
    StaticLoop(Edge i, Edge p, Ts &&...ts) :  seqs(fw<Ts>(ts)...), edges(i, p) {}

    void finalize() {edges.rotate(seqs.finalize());}

    /************************************************************************************/

    auto & iter() {return seqs;};
    bool operator<(StaticLoop const &o) const {return false;}
    bool operator==(StaticLoop const &o) const {return false;}
    /// Return index of sequence that nick is before
    auto nick() const {return seqs.nick();}
    /// Energy of this loop
    auto energy() const {return seqs.energy;}
    /// Sequences delimited by base pairs in this loop
    decltype(auto) sequences() const {return seqs.vec();}
    /// Index of the loop
    auto index() const {return edges.index;};
    /// Index of the parent loop
    auto parent() const {return edges.parent;};
    // Pair of base iterators for the parent base pair
    auto parent_base_pair() const;
    /// Is this loop an exterior loop?
    bool exterior() const {return seqs.exterior();}
    /// Does this loop lack a parent loop?
    bool is_root() const {return edges.is_root();}

    template <class W> static auto edge_getter(W &);
    template <class W> auto strand_index(W const &) const;

    bool next_pair(System const &, usize n, iseq &i, iseq j, PairList const &, vec<SubsequenceList::const_iterator> &);

    /************************************************************************************/
};

/******************************************************************************************/

NUPACK_DEFINE_TEMPLATE(is_loop, StaticLoop, class);

/******************************************************************************************/

template <class SS> template <class W>
auto StaticLoop<SS>::strand_index(W const &w) const {return w.sys->strand_of(seqs.strand_begin());}

/******************************************************************************************/

template <class SS>
auto StaticLoop<SS>::parent_base_pair() const {
    auto b1 = begin_of(seqs.vec()[edges.parent_loc]); decltype(b1) b2;
    if (edges.parent_loc) b2 = end_of(seqs.vec()[edges.parent_loc - 1]) - 1;
    else b2 = end_of(seqs.vec().back()) - 1;
    return std::make_pair(b1, b2);
}

/******************************************************************************************/

template <class SS> template <class W>
auto StaticLoop<SS>::edge_getter(W &w) {
    return [&](Edge i) -> decay<decltype((w[0].edges))> & {return w[i].edges;};
}

/******************************************************************************************/

template <class Loop>
bool is_first_child(Loop const &o, int i) {
    if (o.edges.parent_loc == i) return true; // is first seq
    return (o.exterior() && (o.nick() == i) // seq is after nick
        && (o.edges.parent_loc + 1) % len(o) == o.nick()); // nick is right after parent
}

/******************************************************************************************/

template <class SS, class W>
auto merge_loops(StaticLoop<SS> p, StaticLoop<SS> const &k, W &w) {
    NUPACK_ASSERT(!k.seqs.exterior() || !p.seqs.exterior());
    auto pk_kp = p.edges.template merge<false>(k.edges, StaticLoop<SS>::edge_getter(w));
    auto shift = p.seqs.merge(k.seqs, pk_kp.first, pk_kp.second);
    p.edges.rotate(shift);
    return p;
}

/******************************************************************************************/

// Try to find the next paired base p which is *after* i
// If j is encountered first, return j
template <class SS>
bool StaticLoop<SS>::next_pair(System const &sys, usize n, iseq &i, iseq j, PairList const &pairs, vec<System::StrandIter> &strands) {
    // BEEP("looking for pairs", i, j);
    for (++i; i != j; ++i) {
        if (sys.is_strand_end(i)) {
            seqs.set_last(sys.iterator_at(i)); // Split the views
            if (front(front(seqs)) == Base('_')) {
                // BEEP("root loop done", i, j);
                return false; // single strand, done now
            }
            // BEEP("go to next strand");
            auto s = sys.next_strand_it(j, pairs); // Switch to the new strand
            i = sys.begin_of_strand(s);
            strands.emplace_back(s);
            seqs.append(s->begin());
            edges.append(Ether);
        } else if (pairs[i] != i) { // found a paired base
            seqs.set_last(sys.iterator_at(i + 1)); // Split the views
            seqs.append(sys.iterator_at(pairs[i]));
            edges.append(n);
            // BEEP("found pair", i, j);
            return true;
        }
    }
    // BEEP("found no pairs", i, j);
    seqs.set_last(sys.iterator_at(j + 1));
    return false;
}

/******************************************************************************************/

}
