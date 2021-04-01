#pragma once
#include "DesignComponents.h"
#include "OutputResult.h"
#include "Result.h"

namespace nupack { namespace newdesign {

struct Archive {
    uint max_size {0};
    vec<Result> results;

    Archive() = default;
    Archive(uint size) : max_size(size) {}

    uint remove_dominated();
    uint remove_dominated_by(Result const &res);
    void reevaluate(Local const &env, Designer &designer, uint depth, EnsemblePartition const &part);
    uint update_estimates(Local const &env, Designer &designer, uint depth, EnsemblePartition const &part);
    std::pair<uint, uint> attempt_add(Result const &res);
    bool full() const { return len(results) >= max_size; }

    std::pair<uint, uint> merge(Archive const &other);

    /* distribution */
    vec<real> densities() const;
    real density(Result const &res) const;
    real distance(Result const &res1, Result const &res2) const;

    auto size() const { return len(results); }

    NUPACK_REFLECT(Archive, max_size, results);
};


using ArchiveState = DesignState<Archive>;

} // newdesign
} // nupack
