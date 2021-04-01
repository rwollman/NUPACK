#pragma once
#include "../types/Sequence.h"
#include "../types/Matrix.h"
#include "ParameterStorage.h"

// Loop free energy calculation using a dynamic program-ish method
// Linear time, constant space in the number of base pairs in the loop
// Works by matrix multiplication

namespace nupack {

/******************************************************************************************/

/// Matrix logarithm of 2 x 2 matrix, not tested, simplified in Mathematica, not incorporated yet
/// It could be used to prevent potential overflow on a pathologic multi/exterior loop -- actually not really.
template <class A>
A log2x2(A const &a) {
    using T = value_type_of<A>;
    T const rad = std::sqrt(4 * a(0,1) * a(1,0) + sq(a(0,0) - a(1,1)));
    T const logq = std::log(a(0,0) + a(1,1) - rad);
    T const logp = std::log(a(0,0) + a(1,1) + rad);
    T const dif = (logp - logq) / rad;
    T const add = T(0.5) * ((logp + logq) - T(1.3862943611198906)); // log(4)

    return {{T(0.5) * dif * (a(0,0) - a(1,1)) + add, a(0,1) * dif},
            {a(1,0) * dif, T(0.5) * dif * (a(1,1) - a(0,0)) + add}};
}

/******************************************************************************************/

template <class T, class Plus, class Times>
struct StackMatrix {
    using array = std::array<std::array<T, 2>, 2>;
    array m;

    StackMatrix() = default;
    StackMatrix(std::pair<T, T> const &a, std::pair<T, T> const &b) :
        m{{{{a.first, a.second}}, {{b.first, b.second}}}} {}

    T operator()(bool i, bool j) const {return m[i][j];}

    friend StackMatrix operator*(StackMatrix const &a, StackMatrix const &b) {
        auto dot = [&](auto i, auto j) {
            return Plus()(Times()(a(i, 0), b(0, j)), Times()(a(i, 1), b(1, j)));
        };
        return {{dot(0, 0), dot(0, 1)},
                {dot(1, 0), dot(1, 1)}};
    }

    constexpr T trace() const {return Plus()(m[0][0], m[1][1]);}
};

/******************************************************************************************/

template <class Out, class E>
Out factor(std::size_t l, std::size_t r, bool disable, E const &e) {
    // row (1, 2): this pair (is not / is) stacking towards the left
    // col (1, 2): next pair (is not / is) stacking towards the left
    // I noticed l==3 is same as l > 3 so I combined them.
    if (disable) { // this pair can't stack. used for the fake base pair in exterior loops.
        return              {{e.one(),            e.one()},
                             {e.zero(),           e.zero()}};
    } else if (l == 2) {
        switch (r) {
            case 2:  return {{e.one(),            e.one()},
                             {e.stack(),          e.zero()}};
            case 3:  return {{e.one()+e.right(),  e.one()},
                             {e.stack(),          e.stack()}};
            default: return {{e.one()+e.right(),  e.one()+e.right()},
                             {e.stack(),          e.stack()}};
        }
    } else {
        switch (r) {
            case 2:  return {{e.one(),            e.one()},
                             {e.left(),           e.zero()}};
            case 3:  return {{e.one()+e.right(),  e.one()},
                             {e.left()+e.both(),  e.left()}};
            default: return {{e.one()+e.right(),  e.one()+e.right()},
                             {e.left()+e.both(),  e.left()+e.both()}};
        }
    }
}

/******************************************************************************************/

template <class M>
struct StackModel {
    M const &model;
    Base ip, id, i, j, jd; // base paired to front of first, penultimate of first, last of first, first base of second, second base of second

    real boltz(real e) const  {return model.boltz(e);}
    real stack()       const  {return ip == Base('_') ? zero() : boltz(model.coaxial_stack_energy(ip, id, i, j));}
    real left()        const  {return boltz(model.dG(dangle3, id, i, j));}
    real right()       const  {return boltz(model.dG(dangle5, i, j, jd));}
    real both()        const  {return boltz(model.terminal_mismatch(id, i, j, jd));}

    real one() const {return 1;}
    real zero() const {return 0;}
};

template <class V, class M>
StackModel<M> stack_model(V const &v, std::size_t s0, M const &model) {
    auto const s1 = s0 ? s0 - 1 : v.size() - 1;
    auto const s2 = s1 ? s1 - 1 : v.size() - 1;
    return {model, back(v[s2]), back_index(v[s1], 1), back(v[s1]), front(v[s0]), front(v[s0], 1)};}

/******************************************************************************************/

template <class T, class V, class M>
T stacking_sum(V const &v, int const nick, M const &model) {
    using A = StackMatrix<T, std::plus<T>, std::multiplies<T>>;
    A q = factor<A>(len(back(v)), len(front(v)), nick == 0, stack_model(v, 0, model));
    for (auto s : lrange(1, len(v)))
        q = q * factor<A>(len(v[s-1]), len(v[s]), nick == s, stack_model(v, s, model));
    return q.trace();
}

/******************************************************************************************/

}
