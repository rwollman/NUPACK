#pragma once
#include "Operations.h"
#include "../standard/Vec.h"
#include "../iteration/Range.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include "tbb/partitioner.h"
#include <tbb/task_scheduler_init.h>
#include <tbb/cache_aligned_allocator.h>
#include <thread>
#include <atomic>

namespace nupack {

/******************************************************************************************/

inline auto this_thread_id() {return std::this_thread::get_id();}

/******************************************************************************************/

inline auto default_thread_number() {return std::min<int>(TotalCPU, std::thread::hardware_concurrency());}

/******************************************************************************************/

struct SharedImpl {
    struct State {
        using simple_type = True;
        State() = default;
        State(usize n) : contents(std::make_unique<tbb::task_scheduler_init const>(n)), max(n) {}

        std::unique_ptr<tbb::task_scheduler_init const> contents;
        std::mutex mut;
        usize max;

        auto save_repr() const {return make_members(max);}
        void load_repr(usize m) {max = m; contents = std::make_unique<tbb::task_scheduler_init const>(max);}
    };

    SharedImpl(usize n=0) : state(std::make_shared<State>(n ? n : default_thread_number())) {}

    template <class E, class V, class F>
    bool spread(E &&env, V const &v, GrainSize g, F const &f, OrderedSplit) const {
        std::atomic<usize> count{0u};
#       if TBB_INTERFACE_VERSION < 9100 // not sure exactly when static_partitioner was put in
            auto &&p = tbb::simple_partitioner();
#       else
            auto &&p = tbb::static_partitioner();
#       endif
        tbb::parallel_for(usize(0), usize(len(v)), g.value, [&] (usize i) {f(env, v[count++], i);}, p);
        return false;
    }

    /// Parallel for with manually specified granularity
    template <class E, class V, class F, class P, NUPACK_IF(!is_same<P, OrderedSplit>)>
    bool spread(E &&env, V const &v, GrainSize g, F const &f, P) const {
        tbb::parallel_for(tbb::blocked_range<usize>(0u, len(v), g.value),
            [&](auto const &b) {for (auto i : iterators(b)) f(env, v[i], i);}, tbb::auto_partitioner());
        return false;
    }

    template <class R, class V>
    auto reduce(V const &v, R const &r) const {
        auto range = tbb::blocked_range<const_iterator_of<V>>(begin_of(v), end_of(v));
        auto const acc = DefaultAccumulator();
        return tbb::parallel_reduce(range, value_type_of<V>(),
            [&r, &acc](auto const &a, auto const &o) {return acc(begin_of(a), end_of(a), o, r);}, r);
    }

    template <class E, class V, class F, class T>
    void map(E &&env, V &out, GrainSize g, F const &fun, T tag) const {
        spread(env, out, g, [&](auto &&env, auto const &, auto i) {
            out[i] = fun(env, out[i], i);
        }, tag);
    }

    auto n_workers() const {return state->max;}

    NUPACK_REFLECT(SharedImpl, state);

private:
    std::shared_ptr<State> state;
};

/******************************************************************************************/

}
