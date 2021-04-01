#include <nupack/common/Costs.h>
#include <nupack/standard/Future.h>
#include <nupack/standard/Map.h>
#include <nupack/iteration/Patterns.h>
#include <nupack/iteration/Range.h>
#include <nupack/execution/Local.h>

namespace nupack {

/// Subblocking and naive evaluation costs of all complexes up to lmax
std::array<std::size_t, 3> unit_evaluation_costs(uint n, std::size_t lmax) {
    std::unordered_map<std::string, uint> blocks;
    std::size_t naive_cost = 0, number = 0;
    for (auto l : range(1, lmax+1)) {
        number += compute_necklaces(std::string(l, '\0'), n, [&](auto const &v) {
            throw_if_signal();
            naive_cost += cube(l);
            for (auto i : range(l)) for (auto c : range(l - i)) {
                ++blocks.try_emplace(v.substr(i, c + 1), 0).first->second;
            }
        });
    }
    std::size_t block_cost = sum(blocks, [](auto const &p) {return unit_subblock_cost(len(p.first));});
    std::size_t check = sum(blocks, [](auto const &p) {return p.second * unit_subblock_cost(len(p.first));});
    NUPACK_REQUIRE(check, ==, naive_cost);

    return {number, block_cost, naive_cost};
}

EvaluationCostTable unit_evaluation_cost_table(uint n, real timeout) {
    EvaluationCostTable out(n);
    Local(0).spread(range(n), 1, [&] (Ignore, Ignore, auto n) {
        auto const t0 = std::chrono::high_resolution_clock::now();
        for (uint lmax = 1;
            std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count() < timeout;
            ++lmax
        ) {
            print(n+1, lmax);
            out[n].emplace_back(unit_evaluation_costs(n+1, lmax));
        }
    });
    return out;
}

}
