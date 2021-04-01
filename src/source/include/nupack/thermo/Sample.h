#pragma once
#include "../iteration/Spreadsort.h"
#include "../types/Complex.h"
#include "Algebras.h"
#include "Backtrack.h"
#include "Action.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

using mark_t = vec<uint>;

/** Takes a list of marks and then associates each mark with a random value in
  the interval [0, value] and returns the list sorted by the random values.
*/
inline auto compute_weights(mark_t const & marks, real value) {
    auto w = vmap(marks, [=](auto m) {return std::make_pair(random_float(), m);});
    spreadsort_float(w, first_of, [](auto x, auto y) {return x.first < y.first;});
    for (auto &i : w) i.first *= value;
    return w;
};

/******************************************************************************************/

/** Determines which recursion matrix in block contains element t and then adds
  the corresponding Segment to each of the samples in used.
*/
template <class Block, class Queue, class T, class U>
void push_segment(Block const &block, Queue &queue, T const &t, U used) {
    auto const mems = members_of(block);
    for_each_index(Block::backtracks(), [&](auto I) {
        if (at_c(mems, I).has(t)) {
            auto lims = minmax(at_c(mems, I).indices_of(t));
            Segment s{lims[0], lims[1], at_c(names_of(block), I), -int(decltype(I)::value)};
            // print("    pushing: ", s);
            queue.push(s, mark_t(item_view(used)));
        }
    });
}

/******************************************************************************************/

/** Delegating function to overloaded sample_element that extracts the strands
  included by seg and branches depending on if the segment bridges multiple
  strands
*/
template <class Block, class Queue, class Model>
void sample_element(Block const &block, Complex const &sequence, Model const &model, Queue &queue, Segment const & seg, mark_t const & marks) {
    auto seqs = sequence.strands_included(seg.i, seg.j);
    if (seqs.multi()) sample_element(block, model, queue, seg, marks, MultiStrand(), seqs);
    else sample_element(block, model, queue, seg, marks, SingleStrand(), seqs);
}

/** Core of the sampling algorithm. Finds the matrix element in the block
  refered to by the segment and replays the N4 computation of this element,
  adding new segments to marked samples as the partial sum crosses their
  corresponding weights.
*/
template <class Block, class Queue, class Model, class N, class S>
void sample_element(Block const &block, Model const &model, Queue &queue, Segment const & seg, mark_t const & marks, N, S const &s) {
    using Algebra = SampleAlgebra<typename Model::rig_type>;
    auto elem = get_element(block, seg.i, seg.j, seg.type);
    auto w = compute_weights(marks, mantissa(elem));
    auto weights = view(w);
    for_each_index(Block::backtracks(), [&](auto I) {
        if (at_c(Block::names(), I) != seg.type) return;
        auto const rule = at_c(Block::recursions(), I);

        auto subblock = block.subsquare(span{s.offset, s.offset + len(s)});
        real accum = Model::rig_type::zero();
        Algebra::recurse([&](auto result, auto const &...ts) {
            NUPACK_REQUIRE(len(weights), >, 0);
            auto res = result(-exponent(elem));
            Model::rig_type::plus_eq()(accum, res);
            if (front(weights).first < accum) {
                auto const it = upper_bound(weights.offset(1), accum, first_of);
                auto const used = view(begin_of(weights), it);
                NUPACK_UNPACK(push_segment(block, queue, ts, used));
                weights.set_begin(it);
                return weights.empty(); // short circuit if no more weights
            } else return false;
        }, rule(seg.i-s.offset, seg.j-s.offset, N(), Algebra(), subblock, s, model, DefaultAction()));
    });
}

/******************************************************************************************/

/** Top-level interface function for sampling. This expects a precomputed set of
  recursions matrices (block) corresponding to sequence and a model derived
  from Backtrack_Algebra.
*/
template <class Block, class Model>
std::pair<vec<PairList>, std::size_t> sample_block(Block const &block, Complex const &sequence, Model const &model, uint num_samples=1, bool print_segments=false) {
    if (!num_samples) return {};
    Priority_Queue<Segment, mark_t, typename Segment::Compare> queue;
    vec<PairList> samples(num_samples, PairList(len(sequence)));

    Segment init{0, len(sequence)-1, "Q", -4};
    mark_t marks{indices(samples)};
    queue.push(init, marks);

    std::size_t n = 0;
    for (; !queue.empty(); ++n) {
        auto cur = queue.pop();
        if (print_segments) print("popping: ", cur.first, cur.second);

        /* samples are updated on pop of "B" matrix element. If coaxial stacking and
        dangle states are eventually captured, which will need to be done at
        the time Segments are pushed, this logic should be moved to the same
        point for consistency.
        */
        if (cur.first.type == "B") for (auto i : cur.second) samples[i].add_pair(cur.first.i, cur.first.j);
        sample_element(block, sequence, model, queue, cur.first, cur.second);
    }
    return {std::move(samples), n};
}

/******************************************************************************************/

}}
