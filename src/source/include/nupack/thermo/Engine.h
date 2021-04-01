/**
 * @brief Defines user-level calls to forward and backward dynamic programs
 *
 * @file Engine.h
 * @author Mark Fornace
 * @date 2018-05-31
 * @todo Investigate use of task graph or continuations library to enhance parallelism
 * @todo Look into revamp of LRU--keys could all be added a priori before running calculations
 */
#pragma once

#include "CachedModel.h"
#include "PairProbability.h"
#include "Sample.h"
#include "Subopt.h"
#include "Action.h"

#include "../algorithms/Utility.h"
#include "../reflect/Repr.h"
#include "../types/IO.h"
#include "../thermo/Cache.h"

namespace nupack::thermo {

/**************************************************************************************/

/**
 * @brief Message to send to observer callback on Engine.h functions
 * @tparam T data type of the matrices
 * @tparam Ensemble dangle setting
 * @tparam N complexity
 */
template <class T, class Ensemble, int N>
struct Message {
    using value_type = T;
    decltype(declval<Complex>().views()) sequences; //< reference into the sequences this message was created from
    decltype(declval<BlockMatrix<T, Ensemble, N>>().subsquare(span(0, 0))) block; //< reference into the subblock this message was created from
    real raw_result; //< dynamic program result excluding join penalty and symmetry correction
    real result; //< dynamic program result (including join penalty and symmetry correction)
    char calculation_type; //< enum value of the type of calculation performed
    bool is_subblock; //< is this the complete block or just a subblock of the requested complex?
    NUPACK_REFLECT(Message, sequences, raw_result, result, calculation_type, is_subblock); // block
};

NUPACK_DEFINE_TEMPLATE(is_message, Message, class, class, int);

template <class T, class D, int N> Message<T, D, N> get_message_type(BlockMatrix<T, D, N> const &);
template <class B> using BlockMessage = decltype(get_message_type(declref<B>()));

/**************************************************************************************/

/// Status of calculation holds error descriptor for each sequence on the given diagonal
struct Status {
    small_vec<Stat> errors; //< Stat for each strand: >=0 means that base index diagonal failed
    iseq diagonal = 0; //< the last incomplete diagonal: in terms of sequence indices, not base indices
    Optional<real> result; //< partition function result, if set
    NUPACK_REFLECT(Status, errors, diagonal);
    /// At least one subblock at this diagonal has failed
    bool bad() const {return any_of(errors, [](Stat const &e) {NUPACK_DREQUIRE(e, !=, Stat::ready()); return e != Stat::finished();});}

    void start_diagonal(std::size_t n) {
        errors.resize(n, Stat::ready());
        replace(errors, Stat::finished(), Stat::ready());
    }
    /// Register that a diagonal has finished, return whether an error was detected
    bool finish_diagonal(int o) {
        diagonal = o;
        if (bad()) return true;
        return false;
    }
};

/**************************************************************************************/

/// Calculate a subblock without a cache
/// Return and an int which is -1 if no errors and uplo of what had to be calculated
template <class E, class K, class M, class B, class S, class V, class A>
std::pair<Stat, Region> subblock(E const &env, Stat diag, Region uplo, K const &seq, int i, int j,
                              M const &model, False, B const &, S &sub, V const &pos, A const &action) {
    return {run_block(env, diag, uplo, sub, (j != i), ForwardAlgebra<decltype(model.rig())>(), seq, model, action), uplo};
}

/**
 * @brief Calculate a subblock and return what needed to be calculated: U=upper, L=lower, A=all, C=cached (none)
 * @todo split the mutex out of the cache sometime.
 */
template <class E, class K, class M, class C, class B, class S, class V, class A, NUPACK_IF(!is_same<False, C>)>
std::pair<Stat, Region> subblock(E const &env, Stat diag, Region uplo, K const &seq, int i, int j,
                              M const &model, C &cache, B &block, S &sub, V const &pos, A const &action) {
    if (!cache.limit.satisfiable()) {
        return subblock(env, diag, uplo, seq, i, j, model, False(), block, sub, pos, action);
    } else {
        auto lock = cache.read_lock();
        if (auto const it = cache.find(seq); it != cache.end()) { // computation has finished. read from cache
            auto const done = fork(it->second, [&](auto const &r) -> uint {
                if constexpr (max_log2<value_type_of<decltype(r)>>() > max_log2<value_type_of<B>>()) return 2; // go to ***
                else block.read({pos[i], pos[i+1]}, {pos[j], pos[j+1]}, r);
                return r.complete(); // 1 if record had whole square, 0 if just lower triangle
            });
            if (done == 2) return {Stat(0), uplo}; // *** overflow implied to occur, so say it failed on diagonal 0
            else if (done == 1 || uplo == Region::lower) return {Stat::finished(), Region::cached}; // block is all done
            else uplo = Region::upper; // block is only partially done
        }
    }

    Stat err = run_block(env, diag, uplo, sub, (j != i), ForwardAlgebra<decltype(model.rig())>(), seq, model, action);

    /// Put result in cache if the cache has less information than we computed
    if (err == Stat::finished()) {
        auto lock = cache.write_lock();
        bool already_done = false;
        if (auto const it = cache.find(seq); it != cache.end()) // don't insert if uplo == lower or block is already complete
            already_done = (uplo == Region::lower) || fork(it->second, [](auto const &x) {return x.complete();});
        if (!already_done) {
            cache.insert_or_assign(seq, block.write({pos[i], pos[i+1]}, {pos[j], pos[j+1]}, uplo != Region::lower));
        }
    }
    return {err, uplo};
}

/**************************************************************************************/

namespace detail {
    /// Deduce ensemble type from function, this function is not defined
    template <class T, class Ensemble, int N>
    Ensemble block_ensemble_type(BlockMatrix<T, Ensemble, N> const &);

    template <bool ...Bs, class ...Ms, NUPACK_IF(!sizeof...(Bs))>
    pack<value_type_of<Ms>...> get_data_types(pack<Ms...>);

    template <bool ...Bs, class ...Ms, NUPACK_IF(sizeof...(Ms) == sizeof...(Bs))>
    pack<if_t<Bs, overflow<value_type_of<Ms>>, value_type_of<Ms>>...> get_data_types(pack<Ms...>);

    /// Basically checks that the dangle types of the cache and block agree
    template <class B, class C, NUPACK_IF(!is_same<C, False>)>
    bool_t<(can_assign_from<typename C::mapped_type, BlockRecord<B>>)> check_cache_type(B const &, C const &);

    template <class B, class C, NUPACK_IF(is_same<C, False>)>
    True check_cache_type(B const &, C const);
}

/// Deduce data types to use from a tuple of models and compile-time bools to toggle overflow/nonoverflow
template <class Ms, bool ...Bs>
using DataTypes = decltype(detail::get_data_types<Bs...>(as_pack<decay<Ms>>()));

/**************************************************************************************/

/**
 * @brief Calculate the partition function for a sequence of strands given an initialized block
 * If stat.bad() after this function than overflow occurred
 * @todo eliminate extra computation within a diagonal after a restart
 */
template <class E, class M, class B, class C, class O, class A>
auto run_program(E const &env, Status &stat, Complex const &s, M const &model, B &block, C &cache, O const &observe, A const &action) {
    static_assert(decltype(detail::check_cache_type(block, cache))::value, "Invalid cache type");
    if (!all_of(s, is_canonical))
        NUPACK_ERROR("sequence contains non-canonical nucleotides", s);
    NUPACK_REQUIRE(model.capacity(), >=, len(s));

    using Ens = decltype(detail::block_ensemble_type(block));
    if (!std::holds_alternative<Ens>(ensemble_variant(model.energy_model.ensemble)))
        NUPACK_ERROR("Block and Model have different dangle types");

    auto const list = s.views();
    if (!len(list) || !all_of(list, len)) return model.as_log(model.zero()); // edge cases
    // Prefix sums of the sizes for indexing
    auto const pos = prefixes(true, indirect_view(list, len));
    // Start at diagonal and head for the bottom left
    for (auto o : range(stat.diagonal, s.n_strands())) {
        stat.start_diagonal(s.n_strands() - o); // make room for more blocks
        env.map(stat.errors, 1, [&, observe](auto const &env, Stat &err, auto i) {
            throw_if_signal();
            // Get surrounding block (covers whole square from i to j)
            auto q = block.subsquare({pos[i], pos[i+o+1]});
            auto const k = s.slice(i, i+o+1);

            Region uplo = Region::all;
            std::tie(err, uplo) = subblock(env, err, uplo, k, i, i+o, model, cache, block, q, pos, action);
            if (err == Stat::finished()) {
                auto const r = model.as_log(q.result());
                observe(BlockMessage<B>{k.views(), std::move(q), r,
                    model.complex_result(r, view(list, i, i+o+1)), static_cast<char>(uplo), o+1 != s.n_strands()});
            }
            return err;
        });
        if (stat.finish_diagonal(o)) break;
    }
    auto const q = block.subsquare({0, len(s)}).result();
    auto m = mantissa(q);
    auto const r = model.as_log(q);
    NUPACK_ASSERT(!M::rig_type::prevent_overflow(m), "invalid dynamic program result", s, q);
    return model.complex_result(r, list);
}

/**************************************************************************************/

/**
 * @brief Pair probability (duplicated strands method, see dynamic_program() for common parameters)
 * See dynamic_program() for common arguments.
 * @param stat Calculation status object that holds state for when different data types have to be used
 * @param s Array of the complex sequences and the duplicated complex sequences
 * @return real Partition function or MFE for the complex (with unduplicated sequences)
 */
template <class E, class M, class B, class C, class O, class A>
real run_program(E const &env, Status &stat, std::array<Complex, 2> const &s, M const &model, B &block, C &cache, O const &observe, A const &action) {
    auto const n = s[0].n_strands();

    if (!stat.result) { // PF is not done
        auto result = run_program(env, stat, s[0], model, block, cache, observe, action); // Total matrix block
        if (stat.bad()) return result;
        stat.result.emplace(result);
        stat.diagonal = 0;
        block.copy_square({0, len(s[0])}, {len(s[0]), len(s[1])}); // copy duplicated complex subblock
    }
    fill(stat.errors, Stat::ready());

    auto const list = s[1].views();
    auto const pos = prefixes(true, indirect_view(s[1].views(), len));
    // Start at diagonal, head for the top right corner, stop halfway through
    for (auto o : range(max(1, stat.diagonal), n + 1)) {
        stat.start_diagonal(o); // make room for more blocks
        env.map(stat.errors, 1, [&, observe](auto const &env, Stat &err, auto o2) {
            throw_if_signal();
            auto i = n - o2 - 1;
            auto q = block.subsquare({pos[i], pos[i+o+1]});
            auto const k = s[1].slice(i, i+o+1);
            Region uplo = (o == n) ? Region::lower : Region::all;
            std::tie(err, uplo) = subblock(env, err, uplo, k, i, i+o, model, cache, block, q, pos, action);

            if (err == Stat::finished() && !is_same<O, NoOp>) {
                if (uplo == Region::lower) {
                    observe(BlockMessage<B>{k.views(), std::move(q), 0.0/0.0, 0.0/0.0, static_cast<char>(uplo), true});
                } else {
                    auto const r = model.as_log(q.result());
                    observe(BlockMessage<B>{k.views(), std::move(q), r, model.complex_result(r, view(list, i, i+o+1)), static_cast<char>(uplo), true});
                }
            }

            return err;
        });
        if (stat.finish_diagonal(o)) break;
    }
    return *stat.result;
}

/**************************************************************************************/

/**
 * @brief For each type in Types, calls the given function with a compatible model from models
 * @tparam N Complexity (3 or 4)
 * @tparam Ensemble: recursions to use
 * @tparam Types nupack::pack<> of the possible types
 * @param f Visitor function to apply
 */
template <int N, class Ensemble, class Types, class Ms, class C, class F>
void dispatch_type(Complex const &seq, Types, Ms const &models, C &cache, F &&f) {
    static_assert(Types::size::value >= 1, "Must give at least one data type");
    static_assert(tuple_size<Ms> >= 1, "Must give at least one Model");
    auto Qs = Types::apply([](auto ...ts) { // Make a tuple of Optional<Block>s
        return std::make_tuple(Optional<BlockMatrix<decltype(*ts), Ensemble, N>>()...);
    });
    auto &&c = if_c<is_arithmetic<C>>(cache, [](auto c) {
        return Types::apply([c](auto ...ts) {return Cache<N, Ensemble, decltype(*ts)...>(c);});
    }, Identity());

    Status stat;
    if (while_each_index<Types::size::value>([&](auto I) {
        // Find the model with matching mantissa type to the block type
        auto mod_ok = [](auto x) {return bool_t<(is_same<mantissa_t<decltype(*Types::at(I))>, value_type_of<decltype(*x)>>)>();};
        using M = decltype(find_c(as_pack<decltype(models)>(), mod_ok));
        static_assert(!is_same<M, not_found>, "invalid CachedModel types");
        auto &Q = at_c(Qs, I);
        auto &mod = at_c(models, M());

         // Copy incremental progress from last block if this is not the first type tried
        eval_one(I, [&](auto I) {
            auto &Q0 = at<decltype(I)::value - 1>(Qs);
            Q.emplace(std::move(*Q0));
            Q0.reset();
        }, [&](auto) {
            Q.emplace(seq, mod.zero());
        });

        mod.reserve(len(seq));
        f(stat, *Q, mod, c);
        return stat.bad();
    })) throw std::overflow_error("overflow occurred in dynamic programs for all data types, seq = " + delimited_string(seq.views(), "+"));
}

/**
 * @brief For each type in Types and each possible dangle type, calls the given function with a compatible model from models
 * Essentially dispatches on the dangle type and then calls dispatch_type()
 * @tparam N Complexity (3 or 4)
 * @tparam Types nupack::pack<> of the possible types
 * @param f Visitor function to apply
 */
template <int N, class Types, class Ms, class C, class F>
void dispatch_type_and_dangle(Complex const &seq, Types, Ms const &models, C &cache, F &&f) {
    auto mods = as_tie(models);
    fork(first_of(mods).energy_model.ensemble_type(), [&](auto d) {
        using Ensemble = decltype(d);
        dispatch_type<N, Ensemble>(seq, Types(), mods, cache, [&](auto &stat, auto &Q, auto const &model, auto &cache) {
            using OK = decltype(detail::check_cache_type(Q, cache));
            using requested_type = TypeName<value_type_of<decltype(Q)>>;
            if (!OK::value) NUPACK_ERROR("no provided cache agrees with requested data type and algorithm", requested_type(), type_name(Q), type_name(cache));
            False no_cache;
            f(stat, Q, model, if_c<OK::value>(cache, no_cache));
        });
    });
}

/**************************************************************************************/

/**
 * @brief Return log of partition function or minimum free energy
 *
 * @tparam N=3 Complexity of calculation, i.e. \f$O(N^3)\f$
 * @tparam Bs For each model in models, whether to use overflow type (true) or normal type (false)
 * @param env Environmental execution object (for parallel vs serial execution)
 * @param seq Ordered set of strands
 * @param models A model, std::tie, or std::tuple of models
 * @param cache Either nupack::False() for no cache, or a cache to hold results, e.g. nupack::thermo::LRU
 * @param observe An observer object to call with a nupack::thermo::BlockMessage() from each new subblock calculated
 * @param action An adapter action for the Q_B recursion, e.g. nupack::thermo::PairingAction()
 * @return auto Partition function including join penalties and rotational symmetry correction
 */
template <int N=3, bool ...Bs, class E, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto dynamic_program(E &&env, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    real out = 0;
    dispatch_type_and_dangle<N>(seq, DataTypes<Ms, Bs...>(), models, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
        out = run_program(env, stat, seq, model, Q, cache, observe, action);
    });
    return out;
}

/// Calculate the partition function matrices for a sequence of strands  (see dynamic_program() for common parameters)
template <int N=3, bool ...Bs, class E, class Ensemble, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto block(E &&env, Ensemble, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    auto mods = as_tie(models);
    using Types = DataTypes<Ms, Bs...>;
    auto out = Types::apply([](auto ...ts) {return Optional<MaybeVariant<BlockMatrix<decltype(*ts), Ensemble, N>...>>();});
    real pf;
    dispatch_type<N, Ensemble>(seq, Types(), mods, cache, [&](auto &stat, auto &Q, auto const &model, auto &cache) {
        pf = run_program(env, stat, seq, model, Q, cache, observe, action);
        if (stat.bad()) return;
        out = std::move(Q);
    });
    return std::make_pair(std::move(*out), pf);
}

/**
 * @brief Return structures and their energies (see dynamic_program() for common parameters)
 * @tparam DS=Outer_Stack Algorithm to use
 * @param gap Maximum energy gap
 * @param print_segments Print segments in the stack
 */
template <class Out, int N=3, bool ...Bs, class F, class E, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto block_dependent(E &&env, F &&f, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    std::optional<Out> out;
    dispatch_type_and_dangle<N>(seq, DataTypes<Ms, Bs...>(), models, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
        real pf = run_program(env, stat, seq, model, Q, cache, observe, action);
        if (stat.bad()) return;
        out = f(std::move(Q), seq, model, pf);
    });
    return std::move(out).value();
}

/**
 * @brief Return structures and their energies (see dynamic_program() for common parameters)
 * @tparam DS=Outer_Stack Algorithm to use
 * @param gap Maximum energy gap
 * @param print_segments Print segments in the stack
 */
template <template <class...> class DS=Outer_Stack, int N=3, bool ...Bs, class E, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto subopt(E &&env, real gap, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}, bool print_segments=false) {
    auto out = block_dependent<vec<std::pair<PairList, real>>, N, Bs...>(static_cast<E &&>(env), [&](auto Q, Ignore, auto const &model, Ignore) {
        return subopt_block<DS>(std::move(Q), seq, model, gap, print_segments);
    }, seq, models, cache, observe, action);
    sort(out, [] (auto const &a, auto const &b) {return a.second < b.second; });
    return out;
}

/**
 * @brief Same as subopt() but calls function on each structure as it is found
 * @param f function to call
 */
template <template <class...> class DS=Outer_Stack, int N=3, bool ...Bs, class E, class Ms, class F, class C=False, class O=NoOp, class A=DefaultAction>
real subopt_stream(E &&env, real gap, Complex const &seq, Ms const &models, F &&f, C &&cache={}, O const &observe={}, A const &action={}, bool print_segments=false) {
    real out;
    dispatch_type_and_dangle<N>(seq, DataTypes<Ms, Bs...>(), models, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
        out = run_program(env, stat, seq, model, Q, cache, observe, action);
        if (stat.bad()) return;
        auto it = subopt_iterator<DS>(std::move(Q), seq, model, gap, print_segments);
        while (!it.done()) {
            ++it;
            if (!f(it)) return; // allow early exit
        }
    });
    return out;
}

/**
 * @brief Sample structures from Boltzmann ensemble  (see dynamic_program() for common parameters)
 * @param n_samples number of samples to get
 * @param n_workers number of workers to use in the sampling algorithm
 */
template <int N=3, bool ...Bs, class E, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto sample(E &&env, usize n_samples, usize n_workers, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    if (n_workers == 0) n_workers = env.n_workers();
    std::tuple<vec<PairList>, real, std::size_t> out;
    dispatch_type_and_dangle<N>(seq, DataTypes<Ms, Bs...>(), models, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
        second_of(out) = run_program(env, stat, seq, model, Q, cache, observe, action);
        if (stat.bad()) return;
        if (n_workers == 1) {
            std::tie(first_of(out), third_of(out)) = sample_block(std::move(Q), seq, model, n_samples);
        } else {
            auto v = env.map(n_workers, 1, [&, blk=std::move(Q), n=max(1, ceil(real(n_samples)/n_workers))](auto const &...) {
                return sample_block(blk, seq, model, n);
            });
            first_of(out).reserve(n_samples);
            for (auto &p : v) {
                third_of(out) += p.second;
                first_of(out).insert(first_of(out).end(), move_begin(p.first), move_end(p.first));
            }
        }
    });
    return out;
}

// /**************************************************************************************/

/// Calculate the pair probability for a sequence of strands  (see dynamic_program() for common parameters)
template <int N=3, bool ...Bs, class E, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto pair_probability(E &&env, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    std::pair<Tensor<real, 2>, real> out;
    std::array<Complex, 2> s{seq, seq.duplicated()};

    dispatch_type_and_dangle<N>(s[1], DataTypes<Ms, Bs...>(), models, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
        out.second = run_program(env, stat, s, model, Q, cache, observe, action);
        if (stat.bad()) return;
        out.first = pairs_from_QB<real>(model.rig(), value_of(Q.Q(0, len(s[0])-1)), Q.B.unglued());
    });
    return out;
}


/// Calculate the pair probability for a sequence of strands  (see dynamic_program() for common parameters)
template <int N=3, bool ...Bs, class E, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto bonus_pair_probability(E &&env, Complex const &seq, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}, bool use_B=false) {
    std::pair<Tensor<real, 2>, real> out;
    std::array<Complex, 2> s{seq, seq.duplicated()};

    dispatch_type_and_dangle<N>(s[1], DataTypes<Ms, Bs...>(), models, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
        out.second = run_program(env, stat, s, model, Q, cache, observe, action);
        if (stat.bad()) return;
        NUPACK_ASSERT(std::isfinite(out.second));
        auto divisor = use_B ? value_of(Q.B(0, len(s[0])-1)) : value_of(Q.Q(0, len(s[0])-1));
        out.first = pairs_from_QB<real>(model.rig(), divisor, Q.B.unglued());
    });
    return out;
}

/**************************************************************************************/

/**
 * @brief Calculate the partition functions of all sequences in a given container
 * Even without a cache, there will be savings from subsequence matches
 * See dynamic_program() for common parameters
 * @param sets Ordered sets of strands to measure
 * @param cache Same as dynamic_program(), but if cache is a number, a cache of that max memory size will be used
 */
template <int N=3, bool ...Bs, class E, class V, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto spread(E &&env, V sets, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    using Types = DataTypes<Ms, Bs...>;
    // Start with a map where all sequences have result "undone", sort as biggest first
    auto work = unique_sorted(vmap(sets, [&](auto &i) {
        Complex s{std::move(i)};
        s.rotate_lowest();
        return std::make_pair(std::move(s), std::optional<real>());
    }));
    // Order of iteration so bigger sequences are at front
    auto const order = iter_sort(work, [](auto const &x, auto const &y) {return len(x.first) > len(y.first);});
    // Reserve all models since length is known ahead of time
    auto mods = as_tie(models);
    for_each(mods, [n=order.empty() ? 0 : len(order[0]->first)](auto &m) {m.reserve(n);});

    fork(first_of(mods).energy_model.ensemble_type(), [&](auto d) {
        using Ensemble = decltype(d);
        // Make a new cache if it is given as a max memory size
        // Run this on each set of sequences, putting in subsequence results too
        env.spread(order, 1, [&](auto &&env, auto r, auto const &idx) {
            if (r->second) return;
            dispatch_type<N, Ensemble>(r->first, Types(), mods, cache, [&](auto &stat, auto &Q, auto const &model, auto &&cache) {
                using OK = decltype(detail::check_cache_type(Q, cache));
                if (!OK::value) NUPACK_ERROR("incorrect cache type");
                False no_cache;
                run_program(env, stat, r->first, model, Q, if_c<OK::value>(cache, no_cache), [&](auto const &m) {
                    observe(m);
                    Complex subseqs{m.sequences};
                    subseqs.rotate_lowest();
                    auto it = binary_search(work, subseqs, first_of);
                    if (it != std::cend(work)) it->second.emplace(m.result);
                }, action);
            });
        });
    });
    return vmap(work, [](auto &p) {
        if (!p.second) NUPACK_ERROR("Missing key was not calculated", p.first);
        return std::make_pair(std::move(p.first), *p.second);
    });
}

/**
 * @brief Run all rotationally unique permutations of a set of strands  (see dynamic_program() for common parameters)
 * @param max maximum complex size
 * @param v strands that can be composed (in any order) into the complexes considered
 */
template <int N=3, bool ...Bs, class E, class V, class Ms, class C=False, class O=NoOp, class A=DefaultAction>
auto permutations(E const &env, uint lmax, V const &v, Ms const &models, C &&cache={}, O const &observe={}, A const &action={}) {
    vec<Complex> seqs;
    if (lmax == 0) {
        for_choose_any(false, indices(v), [&](auto const &i) {seqs.emplace_back(indexed_view(i, v));});
    } else {
        while (lmax) compute_necklaces(small_vec<uint>(lmax--), len(v), [&](auto const &i) {seqs.emplace_back(indexed_view(i, v));});
    }
    return spread<N, Bs...>(env, seqs, models, cache, observe, action);
}

/******************************************************************************************/

}

