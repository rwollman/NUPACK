#pragma once
#include "Backtrack.h"
#include "Action.h"
#include "../types/Complex.h"

// mfe information manipulated externally when pushing segments onto stack

namespace nupack { namespace thermo {

template <class T>
struct Partial_Structure : PairList {
    Priority_Queue<Segment, typename Segment::Compare> segments;
    T mfe;
    real32 tiebreaker {random_float<real32>()};

    Partial_Structure() = default;
    Partial_Structure(iseq s) : PairList(s) {}

    struct Compare {
        bool operator()(Partial_Structure const & a, Partial_Structure const & b) const {
            if (a.no_segments() && b.no_segments())
                return std::make_tuple(a.mfe, a.tiebreaker) < std::make_tuple(b.mfe, b.tiebreaker);
            else if (a.no_segments()) return false;
            else if (b.no_segments()) return true;
            else if (a.segments.top() == b.segments.top()) return a.tiebreaker < b.tiebreaker;
            return Segment::Compare()(a.segments.top(), b.segments.top());
        }
    };

    void update_tiebreaker() {tiebreaker = random_float<real32>();}

    void pop(Segment const &seg, T energy) {
        auto blah = segments.pop();
        NUPACK_REQUIRE(blah, ==, seg, seg, energy);
        mfe -= energy;
    }

    template <class Block, class El>
    void push_segment(Block const &block, El const &t) {
        auto mems = members_of(block);
        for_each_index(Block::backtracks(), [&](auto I) {
            if (at_c(mems, I).has(t)) {
                auto lims = minmax(at_c(mems, I).indices_of(t));
                NUPACK_REQUIRE(*at_c(mems, I)(lims[0], lims[1]), ==, value_of(t));
                Segment s{lims[0], lims[1], at_c(names_of(block), I), -int{decltype(I)::value}};
                segments.push(s);
            }
        });
    }

    bool no_segments() const {return segments.empty();}

    void print_segments() const {
        for (auto const &s : segments) {
            std::cout << s << ", ";
        }
        std::cout << std::endl;
    }
};

/******************************************************************************************/

template <class Block, class Queue, class Finished, class Model, class P>
void subopt_element(Block const &block, Complex const &sequence, Model const &model, Queue &queue, Finished &finished, Segment const & seg, P const &p, real cutoff) {
    auto seqs = sequence.strands_included(seg.i, seg.j);
    if (seqs.multi()) subopt_element(block, model, queue, finished, seg, p, MultiStrand(), seqs, cutoff);
    else subopt_element(block, model, queue, finished, seg, p, SingleStrand(), seqs, cutoff);
}

template <class Block, class Queue, class Model, class Finished, class P, class N, class S>
void subopt_element(Block const &block, Model const &model, Queue &queue, Finished &finished, Segment const & seg, P const &p, N, S const &s, real cutoff) {
    auto partial = view(p);
    using A = SuboptAlgebra<typename Model::rig_type>;
    bool found_one = false;

    for_each_index(Block::backtracks(), [&](auto I) {
        if (at_c(Block::names(), I) != seg.type) return;
        auto const rule = at_c(Block::recursions(), I);

        auto select = [&](auto &&result_f, auto const &...ts) {
            auto result = result_f(Zero());

            NUPACK_REQUIRE(len(partial), >, 0);
            if (front(partial).mfe + result < cutoff) {
                found_one = true;
                // following line needed in general because partial structures
                // with same top element can have different mfe values.
                auto it = lower_bound(partial.offset(1), cutoff, [result](auto const &a) {return a.mfe + result;});
                for (auto v : view(begin_of(partial), it)) {
                    v.mfe += result;
                    v.update_tiebreaker();
                    NUPACK_UNPACK(v.push_segment(block, ts));
                    if (v.no_segments()) {finished.push(v);}
                    else queue.push(v);
                }
            }
            return False();
        };

        auto subblock = block.subsquare(span{s.offset, s.offset + len(s)});
        A::recurse(select, rule(seg.i-s.offset, seg.j-s.offset, N(), A(), subblock, s, model, DefaultAction()));
    });
    NUPACK_ASSERT(found_one, "No substructure matched the intermediate MFE value", seg);
}

/******************************************************************************************/


/******************************************************************************************/

template <template <class...> class Queue, class Block, class Model>
struct Subopt_Iterator {
    using element_type = std::pair<PairList, real>;

    Subopt_Iterator(Block const &_block, Complex const &_sequence,
                    Model const &_model, real gap, bool _print_segments) :
            block(_block), sequence(_sequence), model(_model), cutoff(0.0),
            print_segments(_print_segments) {

        // for allowing structures with energy mfe + gap to be included.
        real bump = 1.0e-3; // 1e-4 in NUPACK 3, found that on large sequence it could sometimes fail
        auto total_mfe = get_element(block, 0, len(sequence)-1, "Q");
        cutoff = total_mfe + gap + bump;

        // begin queue
        Segment init {0, len(sequence)-1, "Q", -4};
        Partial_Structure<real> first(len(sequence));
        first.segments.push(init);
        first.mfe = total_mfe;
        queue.push(first);
    }

    bool done() {return fully_specified.empty() && !can_advance();}

    Subopt_Iterator & operator++() {
        next();
        return *this;
    }

    element_type & operator*() {
        return current;
    }

    static element_type sentinel() {
        return element_type{PairList(), std::numeric_limits<real>::infinity()};
    }

private:
    void next() {
        while (fully_specified.empty() && can_advance()) advance();
        if (fully_specified.empty()) return;

        auto s = fully_specified.pop();
        current = {s, model.complex_result(s.mfe, sequence.views())};
    }

    void advance() {
        auto cur = queue.pop();
        auto seg = cur.segments.top();
        auto energy = get_element(block, seg.i, seg.j, seg.type);
        if (print_segments) {
            print("popping: ", seg, "energy: ", energy);
            print("unfinished structures: ", len(queue));
        }

        vec<Partial_Structure<real>> cur_structures = {cur};

        while (queue(seg, queue)) {
            throw_if_signal();
            auto struc = queue.pop();
            cur_structures.push_back(struc);
        }

        for (auto & c : cur_structures) c.pop(seg, energy);

        if (seg.type == "B") for (auto & s : cur_structures) {s.add_pair(seg.i, seg.j);}
        subopt_element(block, sequence, model, queue, fully_specified, seg, cur_structures, cutoff);
    } // this should be the main loop of subopt_block

    bool can_advance() {return !queue.empty();}

    Queue<Partial_Structure<real>, typename Partial_Structure<real>::Compare> queue;
    Stack<Partial_Structure<real>> fully_specified;
    // element_type current;
    element_type current {sentinel()};

    Block const &block;
    Complex const &sequence;
    Model const &model;

    real cutoff;
    bool print_segments;
};

template <template <class...> class Queue, class Block, class Model>
auto subopt_iterator(Block const &block, Complex const &sequence, Model const &model, real gap=0.0, bool print_segments=false) {
    return Subopt_Iterator<Queue, Block, Model>(block, sequence, model, gap, print_segments);
}

template <template <class...> class Queue, class Block, class Model>
auto subopt_block(Block const &block, Complex const &sequence, Model const &model, real gap=0.0, bool print_segments=false) {
    vec<std::pair<PairList, real>> out;
    if (gap < 0 || !std::isfinite(*block.Q(0, len(sequence) - 1))) return out;

    auto it = Subopt_Iterator<Queue, Block, Model>(block, sequence, model, gap, print_segments);
    while (!it.done()) {
        ++it;
        if (*it != decltype(it)::sentinel()) out.push_back(*it);
    }

    return out;
}

std::map<Structure, std::pair<real, real>> unique_subopt(vec<std::pair<PairList, real>>, Complex const &, Model<float> const &);

}}
