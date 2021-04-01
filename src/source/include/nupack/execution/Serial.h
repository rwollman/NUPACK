#pragma once
#include "Operations.h"
#include "../iteration/Patterns.h"
#include "../iteration/Range.h"

namespace nupack {

/******************************************************************************************/

struct SerialImpl {
    template <class E, class V, class F>
    bool spread(E &&env, V const &v, GrainSize, F const &f, Ignore) const {
        usize i = 0;
        for(auto &&x : v) f(env, static_cast<decltype(x) &&>(x), usize(i++));
        return false;
    }

    template <class R, class V>
    auto reduce(V const &v, R const &r) const {
        return std::accumulate(begin_of(v), end_of(v), value_type_of<V>(), r);
    }

    template <class E, class V, class F>
    void map(E &&env, V &out, GrainSize, F const &fun, Ignore) const {
        usize i = 0;
        for (auto &&x : out) x = fun(env, x, usize(i++));
    }

    constexpr auto n_workers() const {return 1u;}
};

/******************************************************************************************/

}
