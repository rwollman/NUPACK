/**
 * @brief Definition of Local executor for shared-memory parallelism
 *
 * @file Local.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Operations.h"

#include "Shared.h"
#include "Serial.h"
#include "../standard/Vec.h"
#include "../standard/Variant.h"

namespace nupack {

/******************************************************************************************/

/// Default constructor Local() makes a serial executor
struct Local {
    constexpr auto even_split() const {return EvenSplit();}

    Variant<SerialImpl, SharedImpl> executor;
    Local(uint n=1) : executor(n == 1 ? SerialImpl() : Variant<SerialImpl, SharedImpl>(SharedImpl(n))) {}

    NUPACK_REFLECT(Local, executor);

    auto n_workers() const {return fork(executor, [](auto const &ex) -> uint {return ex.n_workers();});}

    /**
     * @brief Parallelize a functor across a container
     * @param v: a container, should generally be random access
     * @param g: minimum number of elements to run in a task
     * @param f: functor returning void from (*this, element, index)
     * @param tag: execution tag
     * @return: always false
     */
    template <class V, class F, class T=OrderedSplit, NUPACK_IF(!is_integral<F>)>
    bool spread(V const &v, GrainSize g, F const &f, T tag={}) const {
        return fork(executor, [&](auto const &ex) {return ex.spread(*this, v, g, f, T());});
    }

    /**
     * @brief Parallelize a functor across v which overwrites v
     * @param v: a mutable container, should generally be random access
     * @param g: minimum number of elements to run in a task
     * @param f: functor returning new_element from (*this, element, index)
     * @param tag: execution tag
     */
    template <class V, class F, class T=OrderedSplit, NUPACK_IF(!is_integral<V>)>
    void map(V &out, GrainSize g, F const &fun, T tag={}) const {
        fork(executor, [&](auto const &ex) {ex.map(*this, out, g, fun, tag);});
    }

    /**
     * @brief Parallelize a functor across a range, returning results in a vector
     * @param n: number of elements
     * @param g: minimum number of elements to run in a task
     * @param f: functor returning new_element from (*this, index)
     * @param tag: execution tag
     */
    template <class F, class T=OrderedSplit, NUPACK_IF(!is_integral<F>)>
    auto map(uint n, GrainSize g, F const &fun, T tag={}) const {
        using R = no_ref<decltype(fun(*this, n))>;
        static_assert(!is_same<void, R>, "Functor in nupack::Local::map should not return void");
        vec<R> out(n);
        map(out, g, [&fun](auto &&env, auto &, usize i) {return fun(env, i);}, tag);
        return out;
    }

    template <class R=DefaultReducer, class V>
    auto reduce(V const &v, R const &r=DefaultReducer()) const {
        return fork(executor, [&](auto const &ex) {return ex.reduce(v, r);});
    }

    template <class R=DefaultReducer, class F, class T=OrderedSplit, NUPACK_IF(!is_integral<F>)>
    auto map_reduce(uint n, GrainSize g, F const &f, R const &r={}, T tag={}) const {
        return reduce(map(n, g, f, tag), r);
    }

};


/******************************************************************************************/

}
