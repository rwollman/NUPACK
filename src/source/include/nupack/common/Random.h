/**
 * @brief Random number generation and associated functions
 *
 * @file Random.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Threading.h"
#include "../algorithms/Extents.h"
#include "../iteration/Range.h"
#include <algorithm>
#include <random>

namespace nupack {

/******************************************************************************************/

using FastRNG    = std::minstd_rand;
using DefaultRNG = std::mt19937;
using SlowRNG    = std::ranlux48;

/// Default random number generator
extern thread_local DefaultRNG StaticRNG;
/// Whether random device is used in StaticRNG; same value as NUPACK_RANDOM_DEVICE macro
extern bool const RandomDevice;

/******************************************************************************************/

/// Uniform integer distribution in half-open range [b, e)
template <class B, class E, class T=std::common_type_t<B, E>, NUPACK_IF(is_integral<T>)>
auto uniform_distribution(B b, E e) {return std::uniform_int_distribution<T>(b, e-1);}

/// Uniform float distribution in half-open range [b, e)
template <class B, class E, class T=std::common_type_t<B, E>, NUPACK_IF(is_floating_point<T>)>
auto uniform_distribution(B b, E e) {return std::uniform_real_distribution<T>(b, e);}

/// Wrapper for std::discrete_distribution()
template <class T=uint, class V>
auto discrete_distribution(V &&v) {return std::discrete_distribution<T>(begin_of(v), end_of(v));}

template <class R=void, class D, class RNG=decltype(StaticRNG) &>
auto weighted_samples(D &&dist, std::size_t n, RNG &&gen=StaticRNG) {
    nonvoid<R, std::vector<std::size_t>> picks(dist.max() + 1, 0);
    for (std::size_t i = 0; i != n; ++i) picks[dist(gen)] += 1;
    return picks;
}

/// Shuffle a container
template <class V, class RNG=decltype(StaticRNG) &>
void random_shuffle(V &v, RNG &&gen=StaticRNG) {std::shuffle(begin_of(v), end_of(v), gen);}

/// Choose a random iterator from a container
template <class V, class RNG=decltype(StaticRNG) &>
auto random_choice(V &&v, RNG &&gen=StaticRNG) {
    return std::next(begin_of(v), uniform_distribution(0, len(v))(gen));
}

/// Note: I think random_int_distribution
template <class B, class E, class RNG=decltype(StaticRNG) &>
auto random_range(B b, E e, RNG &&gen=StaticRNG) {return uniform_distribution(b, e)(gen);}

/// Random lower case alphabetical string
template <class RNG=decltype(StaticRNG) &>
string random_string(usize n, RNG &&gen=StaticRNG) {
    string ret(n, 0);
    auto dist = uniform_distribution('a', 'z' + char(1));
    for (auto &i : ret) i = dist(gen);
    return ret;
}

/// Note: I think random_int_distribution
template <class T=real, int N=-1, class RNG=decltype(StaticRNG) &>
auto random_float(RNG &&gen=StaticRNG) {
    static constexpr std::size_t M = (N == -1) ? 2 * sizeof(T) : N;
    return std::generate_canonical<T, M>(gen);
}

template <class RNG=decltype(StaticRNG) &>
bool random_bool(RNG &&gen=StaticRNG) {
    return std::uniform_int_distribution<unsigned short>(0u, 1u)(gen);
}

// void sample_no_replace(int N, int n, vector<int> &samples) {
//     int t = 0; // total input records dealt with
//     int m = 0; // number of items selected so far
//
//     while (m < n) {
//         auto u = random_float(); // call a uniform(0,1) random number generator
//         if ((N - t) * u >= n - m) ++t;
//         else {
//             *(it++)
//             samples[m] = t;
//             ++t; ++m;
//         }
//     }
// }


/******************************************************************************************/

/// Instead of making a shuffled vector, make a vector of shuffled iterators to the input one
template <class V, class RNG=decltype(StaticRNG) &>
auto shuffled_view(V const &v, std::size_t n=*inf, RNG &&rng=StaticRNG) {
    std::vector<const_iterator_of<V>> ret{iterators(v)};
    random_shuffle(ret, rng);
    if (n < len(ret)) ret.erase(begin_of(ret) + n, end_of(ret));
    ret.shrink_to_fit();
    return ret;
}

template <class V, class RNG=decltype(StaticRNG) &>
auto shuffled(V const &v, std::size_t n=*inf, RNG &&rng=StaticRNG) {
    std::vector<value_type_of<V>> ret(begin_of(v), end_of(v));
    random_shuffle(ret, rng);
    if (n < len(ret)) ret.erase(begin_of(ret) + n, end_of(ret));
    ret.shrink_to_fit();
    return ret;
}

/******************************************************************************************/

}
