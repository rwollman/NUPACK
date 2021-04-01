/**
 * @brief Dynamic programming recursions for non-coaxial stacking codes
 *
 * @file BasicPF.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include <atomic>
#include "Algebras.h"
#include "../types/IO.h"
#include "../iteration/View.h"
#include "../iteration/Patterns.h"
#include "../types/Sequence.h"

#define NUPACK_WHERE(cond) !(cond) ? A.zero() : A.maybe()

namespace nupack { namespace thermo {

/// Enum for denoting which part of a block has/had to be calculated
enum class Region : char {upper='U', lower='L', all='A', cached='C'};

struct Stat {
    int value;
    constexpr explicit Stat(int i) : value(i) {}
    NUPACK_REFLECT(Stat, value);

    constexpr bool operator==(Stat const &s) const {return value == s.value;}
    constexpr bool operator!=(Stat const &s) const {return value != s.value;}

    static constexpr Stat ready() {return Stat(-1);}
    static constexpr Stat finished() {return Stat(-2);}

    friend std::ostream &operator<<(std::ostream &os, Stat const &s) {
        if (s.value == -1) return os << "ready";
        if (s.value == -2) return os << "finished";
        return os << "failed(" << s.value << ')';
    }
};

/// Single strand top-level partition function iteration
template <class E, class Seq, class F, class G>
Stat iterate_from_diagonal(E const &env, iseq diag, Region uplo, SingleStrand, Seq const &s, G &&g, F &&f) {
    NUPACK_REQUIRE(uplo, ==, Region::all); // no use case for half done single strand right now
    NUPACK_REQUIRE(diag, <, len(s));

    for (auto const o : indices<iseq>(s)) {
        span is{0u, len(s) - o};
        g(o, is, o > diag); // if o > diag will increment X
        if (o >= diag) {
            if (o % 8 == 0) throw_if_signal();
            std::atomic<bool> err{false};
            env.spread(is, min(10, (len(s)-o) / 4), [&](auto &&, auto i, auto) {
                if (unlikely(f(i, i+o))) err.store(true);
            }, env.even_split());
            if (err.load()) return Stat(o);
        }
    }
    return Stat::finished();
}

/// Multiple strand top-level partition function iteration
template <class E, class Seq, class F, class G>
Stat iterate_from_diagonal(E const &env, iseq diag, Region uplo, MultiStrand, Seq const &s, G &&g, F &&f) {
    span os{(uplo == Region::upper ? s.last_nick() : s.last_nick() - s.first_nick() + 1),
            (uplo == Region::lower ? s.last_nick() : len(s))};
    NUPACK_REQUIRE(diag, <, len(s));

    for (auto const o : os) {
        span is{max(o, s.last_nick()) - o, min(s.first_nick(), len(s) - o)};
        g(o, is, o > diag);
        if (o >= diag) {
            if (o % 8 == 0) throw_if_signal();
            std::atomic<bool> err{false};
            env.spread(is, 1, [&](auto &&, auto i, auto) {
                if (unlikely(f(i, i + o))) err.store(true);
            }, env.even_split());
            if (err.load()) return Stat(o);
        }
    }
    return Stat::finished();
}

/******************************************************************************************/

// On the first or last strand of the multiple strand section
template <class V>
bool on_bread(int i, int j, V const &nicks) {return (i < front(nicks) && j >= back(nicks));}

// For i, j, go over all d, d+1 bases which are on the same strand, d ∈ [i, j-1]
template <class V, class F, class Algebra>
auto sandwich(int i, int j, V const &nicks, Algebra A, F &&f) {
    return A.sum(
        NUPACK_WHERE(on_bread(i+1, j, nicks)) & f(span(i, front(nicks)-1)),
        NUPACK_WHERE(on_bread(i, j-1, nicks)) & f(span(back(nicks), j)),
        NUPACK_WHERE(on_bread(i, j, nicks)) & A.total(iterators(nicks, 0, -1), [f, A](auto n) {
            return NUPACK_WHERE(n[1] - n[0] > 1) & f(span(n[0], n[1]-1));
        })
    );
}

/******************************************************************************************/

NUPACK_DETECT(has_fastiloops, decltype(declval<T>().X));

/// Partition function of [i, j] given that i, j close an extensible interior loop
template <class Block, class Alg, class T, NUPACK_IF(traits::has_fastiloops<Block> && Alg::is_forward::value)>
auto x_loops(int i, int j, Alg A, Block const &Q, T const &) {return A.dot(Q.X[0](i, cspan(8, j-i-5)));}

template <class Block, class Alg, class T, NUPACK_IF(!traits::has_fastiloops<Block> && Alg::is_forward::value)>
auto x_loops(int i, int j, Alg A, Block const &Q, T const &t) {
    return A.total(range(5, min(t.int_max-3, j-i-8)), [=, &Q, &t](auto s) {
        cspan R(5, min(t.int_max+2-s, j-i-s-3));
        return A.dot(Q.YA(R+i, j-s), t.int_size(R+s-2), t.int_asym(s, R));
    });
}

template <class Block, class Alg, class T, NUPACK_IF(!Alg::is_forward::value)>
auto x_loops(int i, int j, Alg A, Block const &Q, T const &t) {
    return A.total(range(10, min(t.int_max+2, j-i-3)), [=, &Q, &t](auto z) {
        return A.total(range(5, z-4), [=, &Q, &t](auto s) {
            return A.product(Q.YA(i+z-s, j-s), t.int_size(z-2), t.int_asym(s, z-s));
        });
    });
}

/******************************************************************************************/

template <class Block, class Alg, class T, NUPACK_IF(traits::has_fastiloops<Block> && Alg::is_forward::value)>
auto x_loops(int i, int j, int m, int n, Alg A, Block const &Q, T const &) {return A.dot(Q.X[0](i, cspan(8, j-n+m-i-2)));}

/// r: number of unpaired on left side + 1
/// s: number of unpaired on right side + 1
template <class Block, class Alg, class T, NUPACK_IF(!traits::has_fastiloops<Block> && Alg::is_forward::value)>
auto x_loops(int i, int j, int m, int n, Alg A, Block const &Q, T const &t) {
    return A.total(range(5, min(t.int_max-3, j-n+1)), [=, &Q, &t] (auto s) {
        cspan R(5, min(t.int_max+2-s, m-i));
        return A.dot(Q.YA(R+i, j-s), t.int_size(R+s-2), t.int_asym(s, R));
    });
}

/// z: r + s
template <class Block, class Alg, class T, NUPACK_IF(!Alg::is_forward::value)>
auto x_loops(int i, int j, int m, int n, Alg A, Block const &Q, T const &t) {
    return A.total(lrange(10, std::min<int>(t.int_max+2, j+m+1-n-i)), [=, &Q, &t](int z) {
        return A.total(lrange(std::max<int>(5, z-j+n), std::min<int>(z-4, m-i)), [=, &Q, &t](int r) {
            return A.product(Q.YA(i+r, j+r-z), t.int_size(z-2), t.int_asym(z-r, r));
        });
    });
}

/******************************************************************************************/

// Partition function of [i, j] given that i, j close a non-special case interior loop
static auto B_extensible = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t) {
        auto asymmetric = [&] (auto z, auto const &Y) { // Large asymmetric loops
            cspan R(i+4+z, j-5);
            return NUPACK_WHERE(i+z+9 < j) & A.sum(A.dot(t.int_rsize(R-i-z-z, z), Y(i+1+z, R+1)),
                                             A.dot(t.int_size(R-i-z-z, z), Y(R+1-z, j-1-z)));
        };
        auto const mm = j > i ? t.mismatch(s[j-1], s[j], s[i], s[i+1]) : A.zero();
        cspan R(4, j-i-5);
        return A.sum(
            NUPACK_WHERE(i+9 < j) & A.product(A.sum(A.dot(Q.T(i+1, R+i+1), t.rbulge(R)),
                                                    A.dot(Q.T(R+i+1, j-1), t.bulge(R))), t.terminal(s[j], s[i])),
            A.product(asymmetric(1, Q.YB), t.mismatch(s[j], s[i])),
            A.product(A.sum(asymmetric(2, Q.YA), asymmetric(3, Q.YA),
                            NUPACK_WHERE(13 < j-i) & x_loops(i, j, A, Q, t)), mm)
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t) {
        auto m = s.first_nick(), n = s.last_nick();

        cspan R(4, j-n), S(4, m-i-1);
        return A.sum(
            NUPACK_WHERE(i+1 < m && n+4 < j) & A.product(A.dot(Q.T(i+1, R+n-4), t.rbulge(R)), t.terminal(s[j], s[i])),
            NUPACK_WHERE(j > n && i+5 < m) & A.product(A.dot(Q.T(S+i+1, j-1), t.bulge(S)), t.terminal(s[j], s[i])),
            A.product(A.sum( // Large asymmetric loops 1
                NUPACK_WHERE(i+5 < m && n+1 < j) & A.dot(t.int_size(S-1, 1), Q.YB(S+i+1, j-2)),
                NUPACK_WHERE(n+4 < j && i+2 < m) & A.dot(t.int_rsize(R-1, 1), Q.YB(i+2, R+n-4))),
                t.mismatch(s[j], s[i])
            ),
            A.product(A.sum( // Large asymmetric loops 2
                NUPACK_WHERE(n+4 < j && i+3 < m) & A.dot(t.int_rsize(R-2, 2), Q.YA(i+3, R+n-4)),
                NUPACK_WHERE(n+4 < j && i+4 < m) & A.dot(t.int_rsize(R-3, 3), Q.YA(i+4, R+n-4)),
                NUPACK_WHERE(i+5 < m && n+2 < j) & A.dot(t.int_size(S-2, 2), Q.YA(S+i+1, j-3)),
                NUPACK_WHERE(i+5 < m && n+3 < j) & A.dot(t.int_size(S-3, 3), Q.YA(S+i+1, j-4)),
                NUPACK_WHERE(i+5 < m && n+4 < j) & x_loops(i, j, m, n, A, Q, t)),
                (j > n && i+1 < m) ? t.mismatch(s[j-1], s[j], s[i], s[i+1]) : A.zero()
            )
        );
    }
);

/******************************************************************************************/

NUPACK_LAMBDA(YA) = [](int i, int j, auto, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    auto const mm = (i != 0 && j+1 != len(s)) ? t.mismatch(s[i-1], s[i], s[j], s[j+1]) : A.zero();
    return A.product(Q.B(i, j), mm);
};

NUPACK_LAMBDA(YB) = [](int i, int j, auto, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    return A.product(Q.B(i, j), t.mismatch(s[i], s[j]));
};

NUPACK_LAMBDA(T) = [](int i, int j, auto, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    return A.product(Q.B(i, j), t.terminal(s[i], s[j]));
};

NUPACK_LAMBDA(D) = [](int i, int j, auto, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    return NUPACK_WHERE(t.can_close(s[i], s[j])) & Q.T(i, j);
};

NUPACK_LAMBDA(dangle) = [](int i, int j, auto, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    return t.dangle(i, j, s);
};

// Partition function of [i, j] given that i, j close a small interior loop (both sides have < 4 unpaired nucleotides)
static auto B_inextensible = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t) {
        return A.total(lrange(i + 1, i+5), [=, &Q, &t, &s](auto d) {
            return A.total(lrange(max(j-4, d + 4), j), [=, &Q, &t, &s](auto e) {
                return A.product(Q.B(d, e), t.interior(view(s, i, d+1), view(s, e, j+1)));
            });
        });
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t) {
        // Small x Small interior loops
        return A.total(lrange(i + 1, min(s.first_nick(), i+5)), [=, &t, &s, &Q](auto d) {
            return A.total(lrange(max(s.last_nick(), j-4), j), [=, &t, &s, &Q](auto e) {
                return A.product(Q.B(d, e), t.interior(view(s, i, d+1), view(s, e, j+1)));
            });
        });
    }
);

/******************************************************************************************/

// Fast interior loops X update for single and multiple strands
static auto X = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return [&, i, j, A] (auto &&x, auto const &x2) {
            if (i+15 < j) {
                cspan K(10, j-i-5);

                auto const y1 = Q.YA(i+5, span(0, len(s)));
                auto const y2 = Q.YA(span(0, len(s)), j-5);
                simd::map(x, K, [&](auto k) {
                    auto const &g = t.int_asym(k);
                    auto const e = -simd::max(exponent_at(x2, k-2), exponent_at(y1, j+3-k), exponent_at(y2, k+i-3));
                    auto const m = A.plus(
                        A.ldexp(A.times(mantissa_at(x2, k-2), t.int_scale(k-2)), exponent_at(x2, k-2) + e), // calculation of mantissa in the exponent of x2
                        A.ldexp(A.times(g, mantissa_at(y1, j+3-k)),              exponent_at(y1, j+3-k) + e),
                        A.ldexp(A.times(g, mantissa_at(y2, k+i-3)),              exponent_at(y2, k+i-3) + e)
                    );
                    return std::make_pair(m, -e);
                });
            }

            if (i+13 < j)  *next(x, 8) = A.times(Q.YA(i+5, j-5), t.int_asym(8));
            if (i+14 == j) *next(x, 9) = A.times(Q.YA(i+6, j-5), t.int_asym(9));
            if (i+14 < j)  *next(x, 9) = A.times(A.plus(Q.YA(i+6, j-5), Q.YA(i+5, j-6)), t.int_asym(9));
            bool out = false;

            if (Debug) for (auto &i : mantissa_view(x)) out |= decltype(A)::rig_type::prevent_overflow(i);
            return out;
        };
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return [&, i, j, A] (auto &&x, auto const &x2) {
            auto const m = s.first_nick(), n = s.last_nick();
            auto const y1 = Q.YA(i+5, span(0, len(s)));
            auto const y2 = Q.YA(span(0, len(s)), j-5);

            if (i+6 < m && n+6 <= j) {
                auto const r = j+4-n, w = m-i+3;
                auto const mm = std::minmax(r, w);

                simd::map(x, 10, mm.first, [&](auto s) {
                    auto const &g = t.int_asym(s);
                    auto const e = -simd::max(exponent_at(x2, s-2), exponent_at(y1, j+3-s), exponent_at(y2, s+i-3));
                    auto const m = A.plus(
                        /*C5*/ A.ldexp(A.times(mantissa_at(x2, s-2), t.int_scale(s-2)), exponent_at(x2, s-2) + e),
                        /*C2*/ A.ldexp(A.times(g, mantissa_at(y1, j+3-s)),              exponent_at(y1, j+3-s) + e),
                        /*C3*/ A.ldexp(A.times(g, mantissa_at(y2, s+i-3)),              exponent_at(y2, s+i-3) + e));
                    return std::make_pair(m, -e);
                });

                simd::map(x, r, w, [&](auto s) {
                    auto const e = -simd::max(exponent_at(y2, s+i-3), exponent_at(x2, s-2));
                    auto const m = A.plus(
                        /*C5*/ A.ldexp(A.times(mantissa_at(x2, s-2), t.int_scale(s-2)), exponent_at(x2, s-2) + e),
                        /*C3*/ A.ldexp(A.times(mantissa_at(y2, s+i-3), t.int_asym(s)),  exponent_at(y2, s+i-3) + e));
                    return std::make_pair(m, -e);
                });

                simd::map(x, w, r, [&](auto s) {
                    auto const e = -simd::max(exponent_at(y1, j+3-s), exponent_at(x2, s-2));
                    auto const m = A.plus(
                        /*C5*/ A.ldexp(A.times(mantissa_at(x2, s-2), t.int_scale(s-2)), exponent_at(x2, s-2) + e),
                        /*C2*/ A.ldexp(A.times(mantissa_at(y1, j+3-s), t.int_asym(s)),  exponent_at(y1, j+3-s) + e));
                    return std::make_pair(m, -e);
                });

                simd::map(x, mm.second, j-n+m-i-2, [&](auto s) {
                    /*C5*/ return std::make_pair(A.times(mantissa_at(x2, s-2), t.int_scale(s-2)), exponent_at(x2, s-2));
                });

                *next(x, 9) = A.times(A.plus(Q.YA(i+6, j-5), Q.YA(i+5, j-6)), t.int_asym(9));
            }

            if (i+6 == m && n+6 <= j) {
                cspan R(9, j+4-n);
                simd::map(x, R, [&](auto s) {
                    return std::make_pair(A.times(mantissa_at(y1, j+3-s), t.int_asym(s)), exponent_at(y1, j+3-s));
                });
            }

            if (i+6 < m && n+5 == j) {
                cspan R(9, m-i+3);
                simd::map(x, R, [&](auto s) {
                    return std::make_pair(A.times(mantissa_at(y2, s+i-3), t.int_asym(s)), exponent_at(y2, s+i-3));
                });
            }

            if (i+5 < m && n+5 <= j) *next(x, 8) = A.times(Q.YA(i+5, j-5), t.int_asym(8));

            bool out = false;
            if (Debug) for (auto &i : mantissa_view(x)) out |= decltype(A)::rig_type::prevent_overflow(i);
            return out;
        };
    }
);

/******************************************************************************************/

/// MS recursion for single and multiple strands
static auto MS0 = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i, j-3);
        return NUPACK_WHERE(i + 3 < j) & A.product(A.dot(Q.D(i, R+4), t.multi3r(R-i)), t.multi2); // only for QB non GU
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(0, j-s.last_nick()+1); // no nicks between d and j
        return A.product(A.dot(t.multi3r(R), Q.D(i, R+s.last_nick())), t.multi2);
    }
);

/// MS recursion for single and multiple strands
static auto MS = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i, j-4);
        return A.product(A.sum(
            NUPACK_WHERE(i+4 < j) & A.dot(t.multi3r(R-i+1), Q.D(i, R+4), Q.dangle(R+5, j)), // only for QB non GU
            NUPACK_WHERE(i+3 < j) & Q.D(i, j) // only for QB non GU
        ), t.multi2);
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(0, j-s.last_nick()); // no nicks between d and j
        return A.product(A.sum(
            Q.D(i, j),
            NUPACK_WHERE(s.last_nick() < j) & A.dot(t.multi3r(R+1), Q.D(i, R+s.last_nick()), Q.dangle(R+s.last_nick()+1, j))
        ), t.multi2);
    }
);

/******************************************************************************************/

/// M recursion for single and multiple strands - no dangles
static auto M0 = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i, j-3), S(i, j-4);
        return A.sum(
            NUPACK_WHERE(i+3 < j) & A.dot(Q.MS(R, j), t.multi3(R-i)),
            NUPACK_WHERE(i+4 < j) & A.dot(Q.M(i, S), Q.MS(S+1, j))
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        // no nicks between i and e
        cspan S(i, s.first_nick());
        return A.sum(
            A.dot(t.multi3(S-i), Q.MS(S, j)),
            sandwich(i, j, s.nicks(), A, [=, &Q](auto R) {return A.dot(Q.M(i, R), Q.MS(R+1, j));})
        );
    }
);

/// M recursion for single and multiple strands
static auto M = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan S(i, j-4);
        return A.sum(
            NUPACK_WHERE(i+3 < j) & Q.MS(i, j),
            NUPACK_WHERE(i+4 < j) & A.dot(Q.dangle(i, S), Q.MS(S+1, j), t.multi3(S+1-i)),
            NUPACK_WHERE(i+4 < j) & A.dot(Q.M(i, S), Q.MS(S+1, j))
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan S(i+1, s.first_nick());
        return A.sum(
            Q.MS(i, j),
            NUPACK_WHERE(i+1 < s.first_nick()) & A.dot(t.multi3(S-i), Q.MS(S, j), Q.dangle(i, S-1)),
            sandwich(i, j, s.nicks(), A, [=, &Q](auto R) {return A.dot(Q.M(i, R), Q.MS(R+1, j));})
        );
    }
);

/******************************************************************************************/

/// Partition function of [i, j] given that i, j close a multiloop
static auto MB = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i+5, j-5);
        bool ok = (t.can_close(s[i], s[j])) && i+10 < j;
        return NUPACK_WHERE(ok) & A.product(A.dot(Q.M(i+1, R), Q.MS(R+1, j-1)), t.multi1, t.multi2, t.terminal(s[j], s[i]));
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return NUPACK_WHERE(t.can_close(s[i], s[j])) & A.product(
            sandwich(i+1, j-1, s.nicks(), A, [=, &Q, &t](auto R) {
                return A.product(A.dot(Q.M(i+1, R), Q.MS(R+1, j-1)), t.multi1, t.multi2);
            }), t.terminal(s[j], s[i])
        );
    }
);

/******************************************************************************************/

/// S recursion for single and multiple strands - no dangles
static auto S0 = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        constexpr bool fast = decltype(A)::is_forward::value;
        return NUPACK_WHERE(i+3 < j) & if_c<fast>(A.sum(Q.S(i, j-1), Q.D(i, j)),
                                                  A.dot(Q.D(i, cspan(i+4, j+1)))); // only for non GU
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        constexpr bool fast = decltype(A)::is_forward::value;
        return if_c<fast>(A.sum(NUPACK_WHERE(s.last_nick() < j) & Q.S(i, j-1), Q.D(i, j)),
                          A.dot(Q.D(i, cspan(s.last_nick(), j+1))));
    }
);

/// S recursion for single and multiple strands
static auto S = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i+4, j);
        return A.sum(
            NUPACK_WHERE(i+3 < j) & Q.D(i, j),
            NUPACK_WHERE(i+4 < j) & A.dot(Q.D(i, R), Q.dangle(R+1 , j))
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(s.last_nick(), j);
        return A.sum(
            Q.D(i, j),
            NUPACK_WHERE(s.last_nick() < j) & A.dot(Q.D(i, R), Q.dangle(R+1, j))
        );
    }
);

/******************************************************************************************/

/// Total partition function of [i, j]
static auto Q = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return A.sum(
            overload([&](auto &&Q) -> decltype(Q.dangle(i, j)) {return Q.dangle(i, j);}, [&](auto const &) {return t.one();})(Q),
            NUPACK_WHERE(i+3 < j) & Q.S(i, j),
            NUPACK_WHERE(i+4 < j) & A.dot(Q.Q(i, cspan(i, j-4)), Q.S(cspan(i+1, j-3), j))
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return A.sum(
            Q.S(i, j),
            sandwich(i, max(j-4, s.last_nick()), s.nicks(), A, [=, &Q](auto R) {
                return A.dot(Q.Q(i, R), Q.S(R+1, j));
            })
        );
    }
);

/******************************************************************************************/

/// Partition function of [i, j] given that i, j are paired and have no nested loops
static auto B_single = overload(
    [](int i, int j, SingleStrand, auto A, auto const &, auto const &s, auto const &t) {
        return t.hairpin(view(s, i, j + 1));
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t) {
        auto const m = s.first_nick(), n = s.last_nick();
        return NUPACK_WHERE(t.can_close(s[i], s[j])) & A.product(A.sum(
            NUPACK_WHERE(j != n && i+1 != m) & A.total(s.nicks(), [=, &Q](auto n) {return A.product(Q.Q(i+1, n-1), Q.Q(n, j-1));}),
            NUPACK_WHERE(j != n && i+1 == m) & Q.Q(m, j-1),
            NUPACK_WHERE(j == n && i+1 != m) & Q.Q(i+1, n-1),
            NUPACK_WHERE(i+1 == j) & t.one() // at the end of strands that border each other
        ), t.terminal(s[j], s[i]));
    }
);

/******************************************************************************************/

/// B matrix - B(i, j) is partition function of [i, j] given that i, j are paired
static auto B = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto &&p) {
        auto sum = [&] {return A.sum(
            B_single(i, j, SingleStrand(), A, Q, s, t),
            B_inextensible(i, j, SingleStrand(), A, Q, s, t),
            Q.MB(i, j),
            B_extensible(i, j, SingleStrand(), A, Q, s, t)
        );};
        return p(i, j, t.can_pair(false, next(s, i), next(s, j)), A, Q, s, t, sum);
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto &&p) {
        auto sum = [&] {return A.sum(
            B_single(i, j, MultiStrand(), A, Q, s, t),
            B_inextensible(i, j, MultiStrand(), A, Q, s, t),
            Q.MB(i, j),
            B_extensible(i, j, MultiStrand(), A, Q, s, t)
        );};
        return p(i, j, t.can_pair(true, next(s, i), next(s, j)), A, Q, s, t, sum);
    }
);

/******************************************************************************************/

}}

#undef NUPACK_WHERE

