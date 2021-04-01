/**
 * @brief Common functors and operations for execution objects
 *
 * @file Operations.h
 * @author Mark Fornace
 * @date 2018-06-01
 */
#pragma once
#include "../standard/Optional.h"

#include <numeric>

namespace nupack {

/// Type representing the parallelization grain size with an integer that must be > 0
struct GrainSize {
    usize value;
    constexpr GrainSize(usize n=1u) : value(n ? n : 1u) {}
};

/******************************************************************************************/

/// Tag for in-order parallelism
struct OrderedSplit {};
/// Tag parallelism that splits range evenly in a fixed manner
struct EvenSplit {};
/// Tag parallelism that splits work to optimize processor affinity
struct AffinitySplit {};

/******************************************************************************************/

// Reduce for optionals; ignores empty values
template <class T, NUPACK_IF(is_optional<T>)>
T reduce(T const &t1, T const &t2) {
    if (t1 && t2) return T(reduce(*t1, *t2));
    else if (t1) return t1;
    else if (t2) return t2;
    else return T();
}

/******************************************************************************************/

namespace traits {
    template <class T, std::size_t ...Is>
    auto reduce_tuple(T const &t, T const &u, indices_t<Is...>) {
        return std::make_tuple(reduce(std::get<Is>(t), std::get<Is>(u))...);
    }
}

/******************************************************************************************/

/// Arithmetic reduction is defined to be addition by default
template <class T>
auto reduce(T const &t, T const &u) -> decltype(t + u) {return t + u;}

/// Construct t from two tuples by calling "reduce" on each element pair
template <class T, NUPACK_IF(is_tuple<T>)>
T reduce(T const &t, T const &u) {return reduce_tuple(t, u, indices_in<T>());}

/******************************************************************************************/

/// Construct t from two tuples by calling "reduce" on each element pair
template <class T, NUPACK_IF(is_pair<T>)>
T reduce(T const &t, T const &u) {
    return {reduce(t.first, u.first), reduce(t.second, u.second)};
}

/******************************************************************************************/

struct DefaultReducer {
    template <class T> T operator() (T const &t1, T const &t2) const {return reduce(t1, t2);};
};

/******************************************************************************************/

struct DefaultAccumulator {
    template <class ...Ts> auto operator() (Ts &&...ts) const {return std::accumulate(fw<Ts>(ts)...);};
};

/******************************************************************************************/

struct DefaultConcatenator {
    template <class T> T operator() (T t1, T const &t2) const {
        t1.insert(end_of(t1), begin_of(t2), end_of(t2)); return t1;
    };
};

/******************************************************************************************/

}
