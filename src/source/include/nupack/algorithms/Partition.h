#pragma once
#include "../standard/Vec.h"
#include "../iteration/Range.h"

namespace nupack::parts {

// Return every possible combination of booleans e.g. [000, 001, 010, 011, 100, 101, 110,]
vec<small_vec<bool>> combinations(uint n);

/******************************************************************************************/

void recurse_subsets(small_vec<std::uint32_t> &p, std::size_t n, std::function<void()> const &f);

// For n=3, yield [], [0], [1], [2], [01], [02], [12], [012],
// The first is always [] and each subset is sorted. Total # shoud be 2^n.
template <class F=NoOp>
std::size_t subsets(std::size_t n, F &&f={}) {
    std::size_t count = 0;
    small_vec<std::uint32_t> p;
    p.reserve(n);
    recurse_subsets(p, n, [&f, &count, &p] {
        f(p);
        ++count;
    });
    return count;
}

/******************************************************************************************/

void recurse_partitions(small_vec<small_vec<std::uint32_t>> &p, bool subset, std::size_t n, std::function<void()> const &f);

// Call the callback f with all partitioning of the indices [0, n). Total # should be Bell number.
template <class F=NoOp>
std::size_t partitions(bool subset, std::size_t n, F &&f={}) {
    std::size_t count = 0;
    small_vec<small_vec<std::uint32_t>> p;
    p.reserve(n);
    recurse_partitions(p, subset, n, [&f, &count, &p] {
        f(p);
        ++count;
    });
    return count;
}

/******************************************************************************************/

struct PairPartition {
    small_vec<std::array<std::uint32_t, 2>> pairs;
    small_vec<std::uint32_t> unpaired;

    PairPartition(std::size_t n) {
        pairs.reserve(n / 2);
        unpaired.reserve(n);
    }
};

void recurse_pairings(PairPartition &p, std::uint32_t n, std::function<void()> const &f);

// Call the callback f with all partitioning of the indices [0, n) into partitions of 1 or 2 vertices.
template <class F=NoOp>
std::size_t pairings(std::uint32_t n, F &&f={}) {
    std::size_t count = 0;
    PairPartition p(n);

    recurse_pairings(p, n, [&f, &count, &p] {f(p); ++count;});
    return count;
}

/******************************************************************************************/

struct BipartitePartition {
    small_vec<std::array<std::uint32_t, 2>> pairs;
    small_vec<std::uint32_t> first_unpaired;
    small_vec<std::uint32_t> second_unpaired;

    BipartitePartition(std::uint32_t m, std::uint32_t n) : first_unpaired(range(m)) {
        pairs.reserve(std::min(m, n));
        second_unpaired.reserve(n);
    }
};

void recurse_bipartite(BipartitePartition &p, std::uint32_t n, std::function<void()> const &f);

// Call the callback f with all bipartite matchings of the indices [0, m) and [0, n)
template <class F=NoOp>
std::size_t bipartite(std::uint32_t m, std::uint32_t n, F &&f={}) {
    std::size_t count = 0;
    BipartitePartition p(m, n);
    recurse_bipartite(p, n, [&f, &count, &p] {f(p); ++count;});
    return count;
}

/******************************************************************************************/

// template <std::size_t N>
// using Link = std::array<std::uint32_t, N>;

// template <std::size_t N>
// struct NpartitePartition {
//     std::tuple<std::array<small_vec<Link<#>>, #>> links; // would be e.g. 3 singles, 3 doubles, 1 triple for N=3

//     NpartitePartition(std::array<std::uint32_t, N> n) {
//         for (auto i : range(N)) std::get<0>(links)[i] = range(n[i]);
//         sort(n);
//         for_each_index<N>([](auto i) {
//             std::get<decltype(i)::value>(links).reserve()
//         });
//     }
// };


// void recurse_npartite(NpartitePartition &p, std::uint32_t n, std::function<void()> &&f) {
//     if (n == 0) {f(); return;}

//     recurse_npartite(p, n-1, [t=n-1, f=std::move(f), &p] {
//         for (auto &ref : p.first_unpaired) { // add to one of the existing subsets
//             auto const u = ref;
//             std::swap(ref, p.first_unpaired.back());
//             p.first_unpaired.pop_back();
//             p.pairs.push_back({u, t});
//             f();
//             p.pairs.pop_back();
//             p.first_unpaired.emplace_back(u);
//             std::swap(ref, p.first_unpaired.back());
//         }
//         p.second_unpaired.emplace_back(t); // make a new subset
//         f();
//         p.second_unpaired.pop_back();
//     });
// }

// template <class F>
// std::size_t npartite(std::uint32_t m, std::uint32_t n, F &&f) {
//     std::size_t count = 0;
//     NpartitePartition p(m, n);
//     recurse_npartite(p, n, [&f, &count, &p] {f(p); ++count;});
//     return count;
// }

}
