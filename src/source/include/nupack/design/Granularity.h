#pragma once
#include "../iteration/Search.h"
#include "../standard/Vec.h"
#include "../common/Error.h"

namespace nupack { namespace newdesign {

struct EnsemblePartition {
    small_vec<bool> mask;
    real deflate;

    EnsemblePartition() = default;

    template <class C>
    EnsemblePartition(vec<C> const &c, real deflate) :
            mask(indirect_view(c, [](auto const &x) { return x.is_on_target(); })),
            deflate(deflate) {};

    auto size() const { return len(mask); }
    auto num_active() const {return count(mask, true);}
    auto num_inactive() const {return count(mask, false);}

    bool all_active() const {return size() == num_active();}

    bool active(uint i) const { return at(mask, i); }
    auto actives() const {
        vec<uint> x;
        izip(mask, [&](auto i, auto m) {if (m) x.emplace_back(i);});
        return x;
    }

    NUPACK_REFLECT(EnsemblePartition, mask, deflate);
};



struct LevelSpecification {
    LevelSpecification() = default;

    std::map<int, uint> exceptions;

    void add_exception(int node, uint depth) {
        if (exceptions.find(node) != end_of(exceptions))
            NUPACK_ERROR("exception already added to LevelSpecification");
        exceptions.emplace(node, depth);
    }

    uint get_depth(int node, uint initial) const {
        auto it = exceptions.find(node);
        if (it == end_of(exceptions)) return initial;   // not an excepted node
        return it->second;                              // node is an exception
    }

    operator bool() const {return len(exceptions) > 0;}

    NUPACK_REFLECT(LevelSpecification, exceptions);
};

struct EnsembleLevelSpecification {
    EnsembleLevelSpecification() = default;

    /** maps complex index to specification of nodes */
    std::map<uint, LevelSpecification> per_complex;
    static LevelSpecification default_spec;

    template <class ...Ts>
    void add_level_spec(uint index, Ts &&...ts) {
        if (per_complex.find(index) != end_of(per_complex))
            NUPACK_ERROR("complex already has a more granular specification.");
        per_complex.emplace(std::piecewise_construct,
                std::forward_as_tuple(index),
                std::forward_as_tuple(fw<Ts>(ts)...));
    }

    LevelSpecification const & get_level_spec(uint index) const {
        auto it = per_complex.find(index);
        if (it == end_of(per_complex)) return default_spec;
        return it->second;
    }

    operator bool() const {return len(per_complex) > 0;}

    NUPACK_REFLECT(EnsembleLevelSpecification, per_complex);
};

}}
