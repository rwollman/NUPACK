/**
 * @brief Dynamic program matrices for non-coaxial stacking codes
 *
 * @file BasicBlock.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "BasicPF.h"
#include "Block.h"
#include "Overflow.h"
#include "Adapters.h"

#include "../model/ModelVariants.h"

namespace nupack::thermo {

struct Matrices_Base {
    /// o is the current diagonal (j-i), possible i in subblock are [a:b)
    template <class Block, class Seq>
    static auto reserve(Block &Q, Seq const &s, iseq o, span is, bool fresh) {}
};

/******************************************************************************************/

template <class E, class Dangle, int N> struct Storage;

/******************************************************************************************/

template <class E, class Dangle> struct Storage<E, Dangle, 4> {
    Symmetric<E> dangle;
    Upper<E>     MB;
    Lower<E>     B;
    Symmetric<E> T;
    Upper<E>     D;
    Symmetric<E> YA;
    Symmetric<E> YB;
    Lower<E>     MS;
    Upper<E>     M;
    Lower<E>     S;
    Upper<E>     Q;
    NUPACK_REFLECT(Storage, dangle, MB, B, T, D, YA, YB, MS, M, S, Q);
    auto alignable() {return std::tie(MB, B, T, D, YA, YB, MS, M, S, Q);}
};

/******************************************************************************************/

template <class E> struct Storage<E, NoStacking, 4> {
    Upper<E>     MB;
    Lower<E>     B;
    Symmetric<E> T;
    Upper<E>     D;
    Symmetric<E> YA;
    Symmetric<E> YB;
    Lower<E>     MS;
    Upper<E>     M;
    Lower<E>     S;
    Upper<E>     Q;
    NUPACK_REFLECT(Storage, MB, B, T, D, YA, YB, MS, M, S, Q);
    auto alignable() {return std::tie(MB, B, T, D, YA, YB, MS, M, S, Q);}
};

/******************************************************************************************/
template <class E, class Dangle> struct Storage<E, Dangle, 3> {
    XTensor<E>  X;
    Symmetric<E> dangle;
    Upper<E>     MB;
    Lower<E>     B;
    Symmetric<E> T;
    Upper<E>     D;
    Symmetric<E> YA;
    Symmetric<E> YB;
    Lower<E>     MS;
    Upper<E>     M;
    Lower<E>     S;
    Upper<E>     Q;

    NUPACK_REFLECT(Storage, X, dangle, MB, B, T, D, YA, YB, MS, M, S, Q);
    auto alignable() {return std::tie(MB, B, T, D, YA, YB, MS, M, S, Q);}
};

/******************************************************************************************/

template <class E> struct Storage<E, NoStacking, 3> {
    XTensor<E>  X;
    Upper<E>     MB; // U
    Lower<E>     B; // L
    Symmetric<E> T; // S
    Upper<E>     D; // U
    Symmetric<E> YA;// S
    Symmetric<E> YB;// S
    Lower<E>     MS;// L
    Upper<E>     M; // U
    Lower<E>     S; // L
    Upper<E>     Q; // U
    NUPACK_REFLECT(Storage, X, MB, B, T, D, YA, YB, MS, M, S, Q);
    auto alignable() {return std::tie(MB, B, T, D, YA, YB, MS, M, S, Q);}
};

template <class D=MinDangles, int N=3> struct Matrices;

template <class Dangles>
struct Matrices<Dangles, 4> : Matrices_Base {
    template <class E> using storage_type = Storage<E, Dangles, 4>;

    template <class E, class V>
    static storage_type<E> storage(Complex const &s, V value) {
        iseq const n = len(s);
        auto m = [=] {return Tensor<E, 2>(n, n, value);};
        return {m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m()};
    }
                        //                              0    1   2  3  4   5   6   7  8  9  10
    static auto recursions() {return std::make_tuple(dangle, MB, B, T, D, YA, YB, MS, M, S, Q);}
    static auto backtracks() {return indices_between<1, 11>();}
    template <class ...Ts> static constexpr void initialize(Ts const &...) {}
};

/******************************************************************************************/

template <>
struct Matrices<NoStacking, 4> : Matrices_Base {
    template <class E> using storage_type = Storage<E, NoStacking, 4>;

    template <class E, class V>
    static storage_type<E> storage(Complex const &s, V value) {
        iseq const n = len(s);
        auto m = [=] {return Tensor<E, 2>(n, n, value);};
        return {m(), m(), m(), m(), m(), m(), m(), m(), m(), m()};
    }
//                                                   0   1  2  3  4   5   6    7   8   9
    static auto recursions() {return std::make_tuple(MB, B, T, D, YA, YB, MS0, M0, S0, Q);}
    static auto backtracks() {return indices_up_to<10>();}
    template <class ...Ts> static constexpr void initialize(Ts const &...) {}
};

/******************************************************************************************/

template <class Dangles>
struct Matrices<Dangles, 3> : Matrices_Base {
    template <class E> using storage_type = Storage<E, Dangles, 3>;

    template <class E, class V>
    static storage_type<E> storage(Complex const &s, V value) {
        iseq const n = len(s);
        auto m = [=] {return Tensor<E, 2>(n, n, value);};
        return {{s, value}, m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m()};
    }

    template <class Block, class Seq>
    static auto reserve(Block &Q, Seq const &s, int o, span is, bool fresh) {
        Matrices_Base::reserve(Q, s, o, is, fresh);
        if (fresh) Q.X.increment();
    }

    template <class Q, class Seq, class Model>
    static void initialize(Q &q, Seq const &s, Model const &t, bool fresh) {if (fresh) q.X.initialize(s, t.zero());}
//                                                   0    1     2  3   4  5  6   7   8   9  10  11
    static auto recursions() {return std::make_tuple(X, dangle, MB, B, T, D, YA, YB, MS, M, S, Q);}

    static auto backtracks() {return indices_between<2, 12>();}
};

/******************************************************************************************/

template <>
struct Matrices<NoStacking, 3> : Matrices_Base {

    template <class E> using storage_type = Storage<E, NoStacking, 3>;

    template <class E, class V>
    static storage_type<E> storage(Complex const &s, V value) {
        iseq const n = len(s);
        auto m = [=] {return Tensor<E, 2>(n, n, value);};
        return {{s, value}, m(), m(), m(), m(), m(), m(), m(), m(), m(), m()};
    }

    template <class Block, class Seq>
    static auto reserve(Block &Q, Seq const &s, iseq o, span is, bool fresh) {
        Matrices_Base::reserve(Q, s, o, is, fresh);
        if (fresh) Q.X.increment();
    }

    template <class Q, class Seq, class Model>
    static void initialize(Q &q, Seq const &s, Model const &t, bool fresh) {if (fresh) q.X.initialize(s, t.zero());}
                                                //   0  1   2  3  4  5   6   7   8    9   10
    static auto recursions() {return std::make_tuple(X, MB, B, T, D, YA, YB, MS0, M0, S0, Q);}

    static auto backtracks() {return indices_between<1, 11>();}
};

/******************************************************************************************/

template <class T, class D, int N> using BlockMatrix = thermo::Block<T, Matrices<D, N>>;

/******************************************************************************************/

/// Double stranded recursion engine
// diag is the starting diagonal, expected to be -1 if this is a fresh calculation or else the
// diagonal which the calculation should resume on.
template <class E, class Block, class Multi, class Seq, class Model, class P, class A>
Stat run_block_body(E const &env, Stat diag, Region uplo, Block &Q, Multi, A, Seq const &s, Model const &t, P &p) {
    NUPACK_ASSERT(diag == Stat::ready() || diag.value >= 0, diag.value);

    Block::initialize(Q, s, t, diag.value <= 0 && uplo != Region::upper); // reinitialize everything if diag was 0 (no progress before)
    auto reserve = [&] (auto ...ts) {return Block::reserve(Q, s, ts...);};
    auto out = iterate_from_diagonal(env, max(0, diag.value), uplo, Multi(), s, reserve, [&](auto i, auto j) {
        bool err = false;
        auto run = overload([](auto const &M, True) {},
            [&](auto &M, auto rule) {
                if (!err) err = M.set(i, j, A(), rule(i, j, Multi(), A(), Q, s, t, p));
            });
        for_each_zip(members_of(Q), Block::recursions(), run);
        return err;
    });
    return out;
}

template <class E, class Block, class Seq, class Model, class P, class A>
Stat run_block(E const &env, Stat diag, Region uplo, Block &Q, bool multi, A, Seq const &s, Model const &t, P &&p) {
    return multi ? run_block_body(env, diag, uplo, Q, MultiStrand(), A(), s, t, p) :
                   run_block_body(env, diag, uplo, Q, SingleStrand(), A(), s, t, p);
}

}
