#pragma once
#include "BasicBlock.h"
#include "CoaxialPF.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

/// Data has some cached information on the sequence, linear storage complexity
template <class T> struct CoaxialRows : public Rows<T, 144> { // 144 = 16 + 64 + 64
    using base_type = typename CoaxialRows::Rows;

private:
    /// Simple multidimensional array access helpers
    static int index0(int i, int j) {return 4 * i + j;} // 16
    static int index1(int i, int j, int k) {return 16 * (i + 1) + 4 * j + k;} // 64
    static int index2(int i, int j, int k) {return 16 * (i + 5) + 4 * j + k;} // 64

    template <class S>
    decltype(auto) at(iseq i, S s) const {return base_type::operator()(i, s);}

    template <class S>
    decltype(auto) at(iseq i, S s) {return base_type::operator()(i, s);}

public:

    template <class ...Ts>
    CoaxialRows(Ts &&...ts) : base_type(fw<Ts>(ts)...) {}

    template <class Seq>
    auto operator() (int r, int i, int k, Seq const &s) const {
        return at(index2(s[i+1], s[k], s[i]), r);
    }

    template <class Seq>
    auto operator() (span R, int i, int k, Seq const &s) const {
        return at(index2(s[i+1], s[k], s[i]), R);
    }

    template <class Seq>
    auto operator() (int i, span R, int j, Seq const &s) const {
        return at(index0(s[i], s[j]), R);
    }

    template <class Seq>
    auto operator() (int i, int j, span R, Seq const &s) const {
        return at(index1(s[i], s[j], s[j+1]), R);
    }

    template <class Model, class Seq> void initialize(Seq const &s, Model const &t) {
        for (auto b : CanonicalBases) for (auto c : CanonicalBases) {
            for (auto r : range(len(s) - 1)) {
                auto value = t.coaxial(b, s[r], s[r+1], c);
                *at(index0(b, c), r) = value;
            }
            if (t.can_pair(b, c)) for (auto d : CanonicalBases) for (auto r : indices(s)) {
                *at(index1(b, c, d), r) = t.coaxial(b, c, d, s[r]);
                *at(index2(b, c, d), r) = t.coaxial(s[r], d, b, c);
            }
        }
    }
};

/******************************************************************************************/

template <class E> struct Storage<E, Stacking, 4> {
    Symmetric<E>  B;
    Symmetric<E>  T;
    Symmetric<E>  D;
    Symmetric<E>  YA;
    Symmetric<E>  YB;
    Upper<E>      MD;
    Upper<E>      MC;
    Lower<E>      MCS;
    Lower<E>      MS; // MS == CMS + MCS; CMS never referenced alone
    Upper<E>      CD;
    Lower<E>      S;
    Symmetric<E>  M;
    Symmetric<E>  Q;
    Symmetric<E>  N;
    CoaxialRows<copy_qualifier<E, mantissa_t<decay<E>>>> coax;
    NUPACK_REFLECT(Storage, B, T, D, YA, YB, MD, MC, MCS, MS, CD, S, M, Q, N, coax);
    auto alignable() {return std::tie(B, T, D, YA, YB, MD, MC, MCS, MS, CD, S, M, Q, N);}
};
/******************************************************************************************/

template <class E> struct Storage<E, Stacking, 3> {
    XTensor<E>   X;
    Symmetric<E>  B;
    Symmetric<E>  T;
    Symmetric<E>  D;
    Symmetric<E>  YA;
    Symmetric<E>  YB;
    Upper<E>      MD;
    Upper<E>      MC;
    Lower<E>      MCS;
    Lower<E>      MS; // MS == CMS + MCS
    Upper<E>      CD;
    Lower<E>      S;
    Symmetric<E>  M;
    Symmetric<E>  Q;
    Symmetric<E>  N;
    CoaxialRows<copy_qualifier<E, mantissa_t<decay<E>>>> coax;
    NUPACK_REFLECT(Storage, X, B, T, D, YA, YB, MD, MC, MCS, MS, CD, S, M, Q, N, coax);
    auto alignable() {return std::tie(B, T, D, YA, YB, MD, MC, MCS, MS, CD, S, M, Q, N);}
};

template <>
struct Matrices<Stacking, 4> : Matrices_Base {
    template <class Block, class Seq, class Model>
    static void initialize(Block &Q, Seq const &s, Model const &t, bool fresh) {
        if (fresh) {
            Q.coax.initialize(s, t);
        }
    }

    template <class E> using storage_type = Storage<E, Stacking, 4>;

    template <class E, class V>
    static storage_type<E> storage(Complex const &s, V value) {
        iseq const n = len(s);
        auto m = [=] {return Tensor<E, 2>{n, n, value};};
        return {m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), {n, value}};
    }

    static auto recursions() {
        return std::make_tuple(coax::B, T, D, YA, YB, coax::MD, coax::MC, coax::MCS, coax::MS,
                coax::CD, coax::S, coax::M, coax::Q, coax::N, True());
    }

    static auto backtracks() {return indices_up_to<14>();} // all except CoaxialRows
};


template <>
struct Matrices<Stacking, 3> : Matrices_Base {
    template <class Block, class Seq, class Model>
    static void initialize(Block &Q, Seq const &s, Model const &t, bool fresh) {
        if (fresh) Q.coax.initialize(s, t);
        if (fresh) Q.X.initialize(s, t.zero());
    }

    template <class E> using storage_type = Storage<E, Stacking, 3>;

    template <class E, class V>
    static storage_type<E> storage(Complex const &s, V value) {
        iseq const n = len(s);
        auto m = [=] {return Tensor<E, 2>{n, n, value};};
        return {{s, value}, m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), m(), {n, value}};
    }

    template <class Block, class Seq>
    static auto reserve(Block &Q, Seq const &s, iseq o, span is, bool fresh) {
        Matrices_Base::reserve(Q, s, o, is, fresh);
        if (fresh) Q.X.increment();
    }

    static auto recursions() {
        return std::make_tuple(X, coax::B, T, D, YA, YB,
                coax::MD, coax::MC, coax::MCS, coax::MS, coax::CD, coax::S, coax::M,
                coax::Q, coax::N, True());
    }

    // backtrack through each matrix except Q.X and Q.coax
    static auto backtracks() {return indices_t<1,2,3,4,5,6,7,8,9,10,11,12,13,14>();}
};

/******************************************************************************************/

}}
