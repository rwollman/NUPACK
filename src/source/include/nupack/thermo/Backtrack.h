/**
 * @brief Backtracking through dynamic programs
 *
 * @file Backtrack.h
 * @author Nick Porubsky
 * @date 2018-05-31
 */
#pragma once
#include "CachedModel.h"
#include "../types/Sequence.h"
#include "../types/PairList.h"
#include "../common/Random.h"
#include "../iteration/Search.h"

#include <map>
#include <set>
#include <utility>

namespace nupack::thermo {

/******************************************************************************************/

/** A priority queue using an underlying map. It expects the value type, V, to
  be Iterable such that if a key exists (by comparison function C), the push
  function will append incoming value to the values in the key already in the
  priority queue.
*/
template <class K, class V, class C=void>
struct Priority_Queue : Iterable<Priority_Queue<K, V, C>> {
    std::map<K, V, C> data;
    using value_type = value_type_of<decltype(data)>;

    void push(K const &key, V const &val) {
        auto it = data.find(key);
        if (it == end_of(data)) data.emplace(key, val);
        else it->second.insert(end_of(it->second), begin_of(val), end_of(val));
    }

    value_type pop() {
        auto el = front(data);
        data.erase(begin_of(data));
        return el;
    }

    value_type const & top() const {return front(data);}

    auto const & iter() const {return data;}
    bool empty() const {return data.empty();}
    void clear() {data.clear();}
    usize size() const {return data.size();}
};


template <class K, class C>
struct Priority_Queue<K, C, void> : Iterable<Priority_Queue<K, C, void>> {
    std::set<K, C> data;
    using value_type = value_type_of<decltype(data)>;

    void push(K const &key) {
        auto it = data.find(key);
        if (it == end_of(data)) data.insert(key);
        // else {NUPACK_ERROR("can't push to queue."); }
    }

    value_type pop() {
        auto el = front(data);
        data.erase(begin_of(data));
        return el;
    }

    value_type const & top() const {return front(data);}

    auto const & iter() const {return data;}
    bool empty() const {return data.empty();}
    void clear() {data.clear();}
    usize size() const {return data.size();}
};

/******************************************************************************************/

template <class K, class ...Ignore>
struct Stack : Iterable<Stack<K, Ignore...>> {
    std::vector<K> data;
    using value_type = value_type_of<decltype(data)>;

    void push(K const &key) {
        data.push_back(key);
    }

    value_type pop() {
        auto el = top();
        data.pop_back();
        return el;
    }

    value_type const & top() const {return back(data);}

    auto const & iter() const {return data;}
    bool empty() const {return data.empty();}
    void clear() {data.clear();}
    usize size() const {return data.size();}
};

/******************************************************************************************/
struct Priority_Queue_Condition {
    template <class S, class Q>
    bool operator()(S const &s, Q const &q) {
        return !q.empty() && s == q.top().segments.top();
    };
};

struct Stack_Condition {
    template <class S, class Q>
    constexpr bool operator()(S const &, Q const &) {return false;};
};

template <class ...Ts>
struct Outer_Priority_Queue : public Priority_Queue<Ts...>, public Priority_Queue_Condition {};

template <class ...Ts>
struct Outer_Stack : public Stack<Ts...>, public Stack_Condition {};
/******************************************************************************************/

/** Represents an element of a recursion matrix that has been sampled. Includes
  matrix type (e.g. "Q", "MS", etc.) priority is determined automatically
  during lookup from the structure of the Block object, and thus Segment works
  without modification for both basic and coaxial matrix sets.
*/
struct Segment {
    usize i, j;
    string type;
    int priority;
    NUPACK_REFLECT(Segment, i, j, type, priority);

    // left base, i, used as an arbitrary tie breaker when comparing segments of the same type and length
    struct Compare {
        bool operator()(Segment const & a, Segment const & b) const {
            return std::make_tuple(-int(len(a)), a.priority, a.i) < std::make_tuple(-int(len(b)), b.priority, b.i);
        }
    };

    bool operator==(Segment const & o) const {
        return std::tie(i, j, type) == std::tie(o.i, o.j, o.type);
    }

    friend std::ostream & operator<<(std::ostream &os, Segment const &t) {
        return os << t.type << ": " << t.i << ", " << t.j << ", priority: " << t.priority;
    }

    usize size() const {return j - i;}

};

/******************************************************************************************/

/** Looks up the matrix entry (i,j) in matrix of type "type" in block.
*/
template <class Block>
auto get_element(Block const &block, int i, int j, string type) {
    auto const mems = members_of(block);
    typename decay<decltype(at_c(mems, size_constant<1>()))>::value_type el {};
    for_each_index(Block::backtracks(), [&](auto I) {
        if (at_c(names_of(block), I) == type) el = value_of(at_c(mems, I)(i, j));
    });
    return el;
}

}
