/** \file Fenwick.h
 * @brief Contains Fenwick class, implementing a Fenwick/binary index tree
 * @todo Fix test error for fenwick.find
 */

#pragma once
#include "../standard/Vec.h"
#include "../common/Error.h"
#include "../common/Random.h"
#include "../algorithms/Utility.h"

namespace nupack {

/******************************************************************************************/

// make the rightmost 1 into a 0
template <class T> inline T pop_one_bit(T t) {return t & plus_one(t);}
// make the rightmost 0 into a 1
template <class T> inline T push_one_bit(T t) {return t | plus_one(t);}


/**
 * @brief Logarithmic Fenwick prefix sum search
 * t: value to search for
 * prefix: typically 0 (the additive identity)
 * total: the total value, used in case of overrun
 * size: number of elements
 * get: map from index to element value
 */
template <class T, class F>
std::pair<uint, T> fenwick_find(T const &t, T prefix, T const &total, int size, F const &get) {
    uint index = 0, mask = next_power_of_two(size);

    while (mask) {
        int cur = index + (mask >>= 1);
        if (cur >= size) continue;

        if (auto next = prefix + get(cur); next < t) {
            index = cur + 1;
            prefix = std::move(next);
        }
    }

    // Possible that a floating point error might push the index over
    if (index == size) {
        // But the entered prefix should have been <= the total still
        NUPACK_REQUIRE(t, <=, total);
        prefix = t;
        --index;
    }
    return {index, t - prefix};
}

/******************************************************************************************/

template <class T, class V=vec<T>>
struct Fenwick : MemberOrdered {
    using is_fenwick = True;
    using container_type = V;
    using value_type = T;
    using size_type = size_type_of<container_type>;
    using difference_type = difference_type_of<container_type>;
    using iterator = iterator_of<container_type>;
    using const_iterator = const_iterator_of<container_type>;

    NUPACK_REFLECT(Fenwick, tree, values, zero_value, total_value);

    /// The actual Fenwick tree values
    container_type tree;
    /// For convenience in updating, hold the marginal values as well
    container_type values;
    /// Value to use as 0, sometimes might not be same as T()
    T zero_value;

protected:
    /// Manually keep track of sum for easy lookup
    T total_value;

    /// Increments value of i-th element by delta. O(log(n))
    void increment_tree(size_type i, T delta) {
        for (; i < len(tree); i = push_one_bit(i)) tree[i] += delta;
    };

    /// Redo position in a tree, assuming it has been set to values[p]
    /// and that tree[:p] has been correctly set (O(1) amortized)
    /// note - must be signed for now since i goes negative
    void extend_to(difference_type const p) {
        for (auto i = minus_one(p); i >= pop_one_bit(p); i = minus_one(pop_one_bit(i))) tree[p] += tree[i];
    }

public:

    // Fenwick() = default;

    /// Initialize tree from zero value
    explicit Fenwick(T t) : zero_value(t), total_value(t) {}

    /// Initialize tree from a vector of values
    explicit Fenwick(T t, container_type v)
        : values(std::move(v)), zero_value(t), total_value(t) {redo_tree();}

    /************************************************************************************/

    /// Bracket operator to marginal values
    value_type operator[](size_type i) const {return values[i];}
    /// First value
    value_type front() const {return values.front();}
    /// Last value
    value_type back() const {return values.back();}
    /// Length of tree
    size_type size() const {return values.size();}
    /// Begin iterator
    iterator begin() {return values.begin();}
    /// Begin const_iterator
    const_iterator begin() const {return values.cbegin();}
    /// End iterator
    iterator end() {return values.end();}
    /// End const_iterator
    const_iterator end() const {return values.cend();}
    /// Resize values and tree
    void resize(size_type s) {values.resize(s, zero_value); tree.resize(s, zero_value);}
    /// Reserve space for values and tree
    void reserve(size_type s) {values.reserve(s); tree.reserve(s);}

    /************************************************************************************/

    /// Add to the value of a position in the tree
    void increment(size_type i, T const &t) {
        values[i] += t;
        total_value += t;
        increment_tree(i, t);
    }
    /// Update the value of a position in the tree
    void update(size_type i, T const &t) {
        if (!Release) NUPACK_REQUIRE(i, <, values.size());
        increment_tree(i, t - values[i]);
        total_value += t - values[i];
        values[i] = t;
    }
    /// Zero the value of a position in the tree
    void zero(size_type i) {update(i, zero_value);}
    /// Swap two positions in the tree
    void swap_pos(size_type i, size_type j) {
        auto delta = values[j] - values[i];
        increment_tree(i, delta); increment_tree(j, -delta);
        swap(values[i], values[j]);
    }

    /************************************************************************************/

    /// Redo the whole tree from the values. O(n)
    void redo_tree() {
        tree.assign(values.begin(), values.end());
        auto const n = len(tree);
        for (size_type i = 0u; i != n; ++i) extend_to(i);
        total_value = sum(values.size());
    };

    /************************************************************************************/

    /// Access sum O(1)
    T total() const {return total_value;}
    /// Sum of all elements of "values" with index < i
    T sum(size_type i) const {
        T out = zero_value;
        while (i--) {out += tree[i]; i = pop_one_bit(i);}
        return out;
    }
    /// Construct all prefix sums - for testing (n log n)
    auto sums() const -> decltype(values) {
        decltype(values) sums;
        auto const n = len(values);
        for (size_type i = 0u; i != n; ++i) sums.emplace_back(sum(plus_one(i)));
        return sums;
    }

    template <class Value, class U=Identity>
    auto find(Value const &t, U const &u={}) const {
        using type = std::decay_t<decltype(u(zero_value))>;
        return fenwick_find<type>(t, u(zero_value), u(total_value), size(), [&](auto i) -> decltype(auto) {return u(tree[i]);});
    }

    /************************************************************************************/

    /// Assign values
    template <class... Args>
    void assign(Args&&... args) {values.assign(fw<Args>(args)...); redo_tree();}
    /// Add a value
    template <class... Args>
    void emplace_back(Args&&... args) {
        auto const p = values.size();
        values.emplace_back(T(fw<Args>(args)...));
        total_value += values[p];
        tree.emplace_back(values[p]); extend_to(p);
    }
    /// Erase last value
    void pop_back() {total_value -= values.back(); values.pop_back(); tree.pop_back();}
    /// Append a range to the back
    template <class It> void extend(It b, It e) {
        values.insert(values.end(), b, e);
        tree.insert(tree.end(), b, e);
        for (int i = len(tree) - e + b; i != len(tree); ++i) {
            extend_to(i);
            total_value += values[i];
        }
    }
    /// Append a range to the back
    template <class V2> void extend(V2 const &v) {extend(v.begin(), v.end());}
    /// Unordered erase by swapping back and position, then popping the back
    void swap_erase(size_type i) {swap_pos(i, minus_one(values.size())); pop_back();}

    template <class ...Ts> void erase(Ts ...ts) {NUPACK_ERROR("not implemented yet");}

    template <class ...Ts> void insert(Ts ...ts) {NUPACK_ERROR("not implemented yet");}

    /******************************************************************************************/

    template <class RNG=decltype(StaticRNG)>
    auto sample(RNG &rng, True replace) const {return find(rng() * total_value).first;}

    // Sample in log(n) time, sort of interesting actually (should check how novel).
    template <class RNG=decltype(StaticRNG)>
    auto sample(RNG &rng, False replace) {auto i = sample(rng, True()); zero(i); return i;}

    /******************************************************************************************/
};

template <class T, class V> void sum(Fenwick<T, V> const &); // not implemented, do not use

/******************************************************************************************/

NUPACK_DETECT(is_fenwick, void_if<T::is_fenwick::value>);

/******************************************************************************************/

}
