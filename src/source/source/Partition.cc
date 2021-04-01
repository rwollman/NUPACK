#include <nupack/algorithms/Partition.h>

namespace nupack::parts {

// Return every possible combination of booleans e.g. [000, 001, 010, 011, 100, 101, 110,]
vec<small_vec<bool>> combinations(uint n) {
    vec<small_vec<bool>> out;
    if (n) {
        out.reserve(1 << n);
        out.emplace_back(n, false);
        for (auto i : range(n)) {
            auto it = out.end();
            out.insert(it, out.begin(), it);
            for (; it != out.end(); ++it) (*it)[n-1-i] = true;
        }
    }
    return out;
}

/******************************************************************************************/

void recurse_subsets(small_vec<std::uint32_t> &p, std::size_t n, std::function<void()> const &f) {
    if (n == 0) {
        f();
    } else if (n == 1) {
        f();
        p.emplace_back(0);
        f();
        p.pop_back();
    } else recurse_subsets(p, n-1, [&f, &p, t=n-1] {
        f();
        p.emplace_back(t);
        f();
        p.pop_back();
    });
}

/******************************************************************************************/

void recurse_partitions(small_vec<small_vec<std::uint32_t>> &p, bool subset, std::size_t n, std::function<void()> const &f) {
    if (n == 0) {
        if (subset) f();
    } else if (n == 1) {
        if (subset) f();
        p.emplace_back();
        p.back().emplace_back(0);
        f();
    } else recurse_partitions(p, subset, n-1, [&f, &p, t=n-1, subset] {
        if (subset) f(); // do not add
        for (auto &x : p) { // add to one of the existing subsets
            x.emplace_back(t);
            f();
            x.pop_back();
        }
        p.emplace_back();
        p.back().emplace_back(t); // make a new subset
        f();
        p.pop_back();
    });
}

/******************************************************************************************/

void recurse_pairings(PairPartition &p, std::uint32_t n, std::function<void()> const &f) {
    if (n == 0) {
        // f();
    } else if (n == 1) {
        p.unpaired.emplace_back(0);
        f();
    } else recurse_pairings(p, n-1, [t=n-1, &f, &p] {
        for (auto &u : p.unpaired) { // add to one of the existing subsets
            auto const i = u;
            swap(u, p.unpaired.back());
            p.unpaired.pop_back();
            p.pairs.push_back({i, t});
            f();
            p.pairs.pop_back();
            p.unpaired.emplace_back(i);
            swap(u, p.unpaired.back());
        }
        p.unpaired.emplace_back(t); // make a new subset
        f();
        p.unpaired.pop_back();
    });
}

/******************************************************************************************/

void recurse_bipartite(BipartitePartition &p, std::uint32_t n, std::function<void()> const &f) {
    if (n == 0) {
        f();
    } else recurse_bipartite(p, n-1, [t=n-1, &f, &p] {
        for (auto &ref : p.first_unpaired) { // add to one of the existing subsets
            auto const u = ref;
            swap(ref, p.first_unpaired.back());
            p.first_unpaired.pop_back();
            p.pairs.push_back({u, t});
            f();
            p.pairs.pop_back();
            p.first_unpaired.emplace_back(u);
            swap(ref, p.first_unpaired.back());
        }
        p.second_unpaired.emplace_back(t); // make a new subset
        f();
        p.second_unpaired.pop_back();
    });
}

/******************************************************************************************/

}
