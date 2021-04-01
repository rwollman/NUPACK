#pragma once
#include "../types/Sequence.h"
#include "ParameterStorage.h"

// Classes to iterate through all stacking states of a particular multi/exterior loop
// The way it works is to consider all possibilities for a given base pair, then recursively
// consider the next base pair's possibilities, etc. There is a bit of wrap-around logic for
// the first base pair listed.
// There could be some optimization to cut the sub-ensembles which are separated by more than
// one base, or some mini-dynamic program could also be done presumably too.

namespace nupack {

/******************************************************************************************/

struct LoopStackingState {
    enum Stack : std::uint_fast8_t {None, LeftDangle, RightDangle, BothDangle, LeftStack, RightStack, Disabled};

    small_vec<Stack> state;

    LoopStackingState(std::size_t n, int nick) : state(n, None) {
        if (nick >= 0) state[nick] = Disabled;
    }

    template <class F, class V>
    void enumerate(V const &v, F &&f, unsigned const i=0) {
        // This is the functor to go onto the next base pair
        auto const recurse = [&, i](Stack s) {
            switch (s) {
                case None: {
                    enumerate(v, f, i+1);
                    return;
                }
                case RightStack: {
                    auto tmp = f.add_stack(v, next(v, i));
                    NUPACK_DREQUIRE(state[i], ==, None, i, state);
                    NUPACK_DREQUIRE(state[i+1], ==, None, i, state);
                    state[i] = RightStack; state[i+1] = Disabled;
                    enumerate(v, f, i+1);
                    state[i] = None; state[i+1] = None;
                    f.remove(std::move(tmp));
                    return;
                }
                case LeftStack: {
                    NUPACK_DREQUIRE(i, ==, 0, v, s);
                    auto tmp = f.add_stack(v, v.end()-1);
                    NUPACK_DREQUIRE(state[0], ==, None, i, state);
                    NUPACK_DREQUIRE(state.back(), ==, None, i, state);
                    state[0] = LeftStack; state.back() = Disabled;
                    enumerate(v, f, i+1);
                    state[0] = None; state.back() = None;
                    f.remove(std::move(tmp));
                    return;
                }
                default: {
                    auto tmp = f.add_dangle(v, next(v, i), s);
                    NUPACK_DREQUIRE(state[i], ==, None, i, state);
                    state[i] = s;
                    enumerate(v, f, i+1);
                    state[i] = None;
                    f.remove(std::move(tmp));
                    return;
                }
            }
        };

        if (i == state.size()) {
            f(static_cast<small_vec<Stack> const &>(state));
        } else if (i == 0) {
            recurse(None); // no stacking
            if (state[0] == Disabled) return;

            if (len(v[0]) == 2) {
                if (state[1] != Disabled) recurse(RightStack);
            } else {
                recurse(RightDangle);
            }

            if (len(back(v)) == 2) {
                if (state.back() != Disabled) recurse(LeftStack);
            } else {
                recurse(LeftDangle);
                if (len(v[0]) != 2)
                    recurse(BothDangle);
            }

        } else if (i+1 == state.size()) {
            recurse(None); // Disabled or None
            if (state[i] == Disabled) return;

            if (len(v[i-1]) > 3 || (len(v[i-1]) == 3 && state[i-1] != RightDangle && state[i-1] != BothDangle)) {
                recurse(LeftDangle);
                if (len(v[i]) > 3 || (len(v[i]) == 3 && state[0] != LeftDangle && state[0] != BothDangle))
                    recurse(BothDangle);
            }

            if (len(v[i]) > 3 || (len(v[i]) == 3 && state[0] != LeftDangle && state[0] != BothDangle)) {
                recurse(RightDangle);
            }

        } else {
            recurse(None);

            if (state[i] == Disabled) return;

            if (len(v[i-1]) > 3 || (len(v[i-1]) == 3 && state[i-1] != RightDangle && state[i-1] != BothDangle)) {
                recurse(LeftDangle);
                if (len(v[i]) != 2)
                    recurse(BothDangle);
            }

            if (len(v[i]) != 2)
                recurse(RightDangle);

            if (state[i+1] != Disabled && len(v[i]) == 2)
                recurse(RightStack);
        }
    }

    template <class F>
    struct Simple {
        F &&f;
        bool add_dangle(Ignore, Ignore, Ignore) {return {};}
        bool add_stack(Ignore, Ignore) {return {};}
        void remove(Ignore) {}
        void operator()(small_vec<Stack> const &p) const {f(p);}
    };

    template <class V, class F>
    void enumerate_simple(V const &v, F &&f) {
        return enumerate(v, Simple<F>{std::forward<F>(f)});
    }
};

/******************************************************************************************/

inline std::string_view name(LoopStackingState::Stack s) {
    switch (s) {
        case LoopStackingState::None:        return "0";
        case LoopStackingState::LeftDangle:  return "ld";
        case LoopStackingState::RightDangle: return "rd";
        case LoopStackingState::BothDangle:  return "b";
        case LoopStackingState::LeftStack:   return "ls";
        case LoopStackingState::RightStack:  return "rs";
        case LoopStackingState::Disabled:    return "x";
        default:                             return "";
    }
}

inline std::ostream &operator<<(std::ostream &os, LoopStackingState::Stack s) {return os << name(s);}

/******************************************************************************************/

template <class F, class Model>
struct StackingStateEnergy {
    F emit;
    Model const &model;
    real energy = 0;

    template <class S>
    real add_stack(S const &s, typename S::const_iterator l) {
        NUPACK_DREQUIRE(l->size(), ==, 2, *l, s);
        auto const prev = cyclic_prev(s, l), next = cyclic_next(s, l);
        auto e = model.coaxial_stack_energy(back(*prev), front(*l), back(*l), front(*next));
        energy += e;
        return e;
    }

    template <class S>
    real add_dangle(S const &s, typename S::const_iterator l, LoopStackingState::Stack t) {
        auto const prev = cyclic_prev(s, l);
        real e = 0;
        if (t == LoopStackingState::LeftDangle) {
            NUPACK_DREQUIRE(prev->size(), >, 2, *prev, *l, s, t);
            e = model.dG(dangle3, back_index(*prev, 1), back(*prev), front(*l));
        } else if (t == LoopStackingState::RightDangle) {
            NUPACK_DREQUIRE(l->size(), >, 2, *prev, *l, s, t);
            e = model.dG(dangle5, back(*prev), front(*l), front(*l, 1));
        } else if (t == LoopStackingState::BothDangle) {
            NUPACK_DREQUIRE(prev->size(), >, 2, *prev, *l, s, t);
            NUPACK_DREQUIRE(l->size(), >, 2, *prev, *l, s, t);
            e = model.terminal_mismatch(back_index(*prev, 1), back(*prev), front(*l), front(*l, 1));
        }

        energy += e;
        return e;
    }

    void remove(real r) {energy -= r;}

    template <class V>
    void operator()(V const &v) const {emit(v, energy);}
};

/******************************************************************************************/

template <class V, class F, class M>
void enumerate_stacking_state_energies(V const &sequences, int nick, M const &model, F &&f) {
    LoopStackingState state(sequences.size(), nick);
    StackingStateEnergy<F, M> callback{std::forward<F>(f), model};
    state.enumerate(sequences, callback);
}

/******************************************************************************************/

inline char loop_stack_letter(LoopStackingState::Stack c) {
    switch (c) {
        case LoopStackingState::None:         return 'n';
        case LoopStackingState::LeftDangle:   return 'l';
        case LoopStackingState::RightDangle:  return 'r';
        case LoopStackingState::LeftStack:    return 'l';
        case LoopStackingState::RightStack:   return 'r';
        case LoopStackingState::BothDangle:   return 'b';
        case LoopStackingState::Disabled:     return 'n';
        default: NUPACK_ERROR("invalid stacking state enumerator");
    }
}

/******************************************************************************************/

string loop_stack_string(small_vec<LoopStackingState::Stack> const &v);

string loop_stack_sequence_string(small_vec<LoopStackingState::Stack> const &v);

/******************************************************************************************/

}
