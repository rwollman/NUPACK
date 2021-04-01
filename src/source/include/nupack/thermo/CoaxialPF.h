#pragma once
#include "BasicPF.h"

#define NUPACK_WHERE(cond) !(cond) ? A.zero() : A.maybe()

namespace nupack { namespace thermo {

namespace coax {

/******************************************************************************************/

template <class Algebra, class F> auto dangle_sum(Algebra A, F &&f) {
    // scans through k, l in range(2), range(2) - compile time constants may or may not help
    return A.sum(f(0, 0), f(0, 1), f(1, 0), f(1, 1));
}

static auto CD = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i, j);
        return A.sum(
            dangle_sum(A, [=, &Q, &s, &t](auto k, auto l) {
                return NUPACK_WHERE(j >= l && i + k < j - l) & A.product(t.dangle(i, i+k, j-l, j, s), Q.D(i+k, j-l));
            }),
            NUPACK_WHERE(i < j) & A.dot(Q.coax(i, R, j, s), Q.D(i, R), Q.D(R+1, j))
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return A.sum(
            sandwich(i, j, s.nicks(), A, [=, &Q, &s](auto R) {
                return A.dot(Q.coax(i, R, j, s), Q.D(i, R), Q.D(R+1, j)); // ICS 6
            }),
            dangle_sum(A, [=, &Q, &s, &t](auto k, auto l) {
                return NUPACK_WHERE(on_bread(i + k, j - l, s.nicks())) & A.product(t.dangle(i, i+k, j-l, j, s), Q.D(i+k, j-l)); // D 5
            })
        );
    }
);


static auto MD = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return dangle_sum(A, [=, &Q, &s, &t](auto k, auto l) {
            return NUPACK_WHERE(j >= l && i + k < j - l) & A.product(Q.D(i+k, j-l), t.multi3(k+l), t.dangle(i, i+k, j-l, j, s), t.multi2);
        });
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return dangle_sum(A, [=, &Q, &s, &t](auto k, auto l) {
            return NUPACK_WHERE(on_bread(i + k, j - l, s.nicks()))
&                 A.product(Q.D(i+k, j-l), t.multi3(k+l), t.dangle(i, i+k, j-l, j, s), t.multi2); // MD 5
        });
    }
);


static auto MC = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i, j);
        return NUPACK_WHERE(i < j) & A.product(A.dot(Q.coax(i, R, j, s), Q.D(i, R), Q.D(R+1, j)), t.multi22);
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return sandwich(i, j, s.nicks(), A, [=, &Q, &s, &t](auto R) {
            return A.product(A.dot(Q.coax(i, R, j, s), Q.D(i, R), Q.D(R+1, j)), t.multi22); // MICS 6
        });
    }
);


NUPACK_LAMBDA(MCS) = [](int i, int j, auto N, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    cspan R(decltype(N)::value ? s.last_nick() : i, j+1);
    return A.dot(t.multi3r(j+1-R), Q.MC(i, R)); // not sure about j+1-R
};


NUPACK_LAMBDA(MS) = [](int i, int j, auto N, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    cspan R(decltype(N)::value ? s.last_nick() : i, j+1);
    return A.sum(Q.MCS(i, j), A.dot(t.multi3r(j+1-R), Q.MD(i, R))); // 29
};

static auto M = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan R(i, j+1), S(i, j);
        return A.sum(
            A.dot(t.multi3(R-i), Q.MS(R, j)),
            NUPACK_WHERE(i < j) & A.dot(Q.M(i, S), Q.MS(S+1, j))
        );
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        cspan S(i, s.first_nick());
        return A.sum(
            A.dot(t.multi3(S-i), Q.MS(S, j)), // 38
            sandwich(i, j, s.nicks(), A, [=, &Q](auto R) {return A.dot(Q.M(i, R), Q.MS(R+1, j));}) // 40
        );
    }
);

/// Q.N (exterior loop) recursions
NUPACK_LAMBDA(N) = [](int i, int j, auto N, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
    return NUPACK_WHERE(decltype(N)::value) & A.total(s.nicks(), [=, &Q](auto n) {
        return A.product(Q.Q(i, n-1), Q.Q(n, j)); // 43
    });
};

/// Q.B contributions from closing pair dangles
NUPACK_LAMBDA(B_cpd) = [](int i, int j, auto A, auto const &Q, auto const &s, auto const &t) {
    return dangle_sum(A, [=, &Q, &s, &t](auto k, auto l) {
        cspan R(i+k+1, s.first_nick());
        auto bread = [ns=s.nicks()](auto i, auto j) {return on_bread(i, j, ns);};
        auto dang = t.dangle(j-l, j, i, i+k, s);
        return A.sum(
            A.product(A.sum(
                NUPACK_WHERE(bread(i+k+1, j-l-1)) & // MCPD 5
                    A.product(A.dot(t.multi3(R+l-i-1), Q.MCS(R, j-l-1)), dang),

                NUPACK_WHERE(bread(i+k+1, j-l-1)) &  // MCPD 7
                    sandwich(i+1, j-1, s.nicks(), A, [=, &Q, &t](auto R) {
                        return A.product(A.dot(Q.M(i+k+1, R), Q.MS(R+1, j-l-1)),
                                       t.multi3(k+l), dang);
                    })
            ), t.multi12),
            // ECPD 6 Both sides dangle
            NUPACK_WHERE(bread(i+k+1, j-l-1))
                & A.product(dang, Q.N(i+k+1, j-l-1)),

            // ECPD 8 Only j side dangles
            NUPACK_WHERE(bread(i+k, j-l-1) && i+k+1 == s.first_nick())
                & A.product(dang, Q.Q(i+k+1, j-l-1)),

            // ECPD 8 Only i side dangles
            NUPACK_WHERE(bread(i+k+1, j-l) && j-l == s.last_nick())
                & A.product(dang, Q.Q(i+k+1, j-l-1)),

            // ECPD 8 Neither side dangles but they're consecutive
            NUPACK_WHERE(i+k+1 == s.first_nick() && i+k+1 == j-l && j-l == s.last_nick())
                & A.product(dang, t.one()) // Q.Q(i+k+1, j-l-1)=1 (empty recursion)
        );
    });
};

// Q.B contributions from closing pair stacking
NUPACK_LAMBDA(B_cps) = [](int i, int j, auto A, auto const &Q, auto const &s, auto const &t) {
    static constexpr int const k = 0, m = 0, l = 0; // have to give sandwich template argument if k != 0 or t != 0
    auto sand = [&, ns=s.nicks()](auto i, auto j, auto f) {return sandwich(i, j, ns, A, f);};
    auto bread = [ns=s.nicks()](auto i, auto j) {return on_bread(i, j, ns);};
    return A.sum(
        sand(i+l+1, j-k-1, [=, &Q, &s, &t](auto R) { // MCPCS 6
            return A.product(A.dot(Q.coax(j, i, R, s), Q.D(i+1, R), Q.M(R+1, j-1)), t.multi122); // multi3(m+k+l)
        }),
        sand(i+m+1, j-l-1, [=, &Q, &s, &t](auto R) { // MCPCS 8
            return A.product(A.dot(Q.coax(R+1, j-1, i, s), Q.M(i+1, R), Q.D(R+1, j-1)), t.multi122); // multi3(m+k+l)
        }),
        sand(i+l+1, j-k-1, [=, &Q, &s](auto R) { // ECPCS 8
            return A.dot(Q.coax(j, i, R, s), Q.D(i+l+1, R), Q.N(R+m+1, j-k-1));
        }),
        sand(i+m+1, j-l-1, [=, &Q, &s](auto R) { // ECPCS 14
            return A.dot(Q.coax(R+1, j-1, i, s), Q.N(i+m+1, R), Q.D(R+k+1, j-l-1));
        }),
        NUPACK_WHERE(bread(i+l+1, j-k-1)) & A.total(s.nicks(), [=, &Q, &s](auto n) { // ECPCS 10B
            return A.product(Q.coax(j, i, n-m-1, s), Q.D(i+l+1, n-m-1), Q.Q(n, j-k-1));
        }),
        NUPACK_WHERE(bread(i+m+1, j-l-1)) & A.total(s.nicks(), [=, &Q, &s](auto n) { // ECPCS 16B
            return A.product(Q.coax(n+k, j-1, i, s), Q.Q(i+m+1, n-1), Q.D(n+k, j-l-1));
        }),
        NUPACK_WHERE(i+m+1 == s.first_nick() && j-l-1 >= s.last_nick()) & sand(i+m, j-l-1, [=, &Q, &s](auto R) {
            return A.dot(Q.coax(R+1, j-1, i, s), Q.Q(i+m+1, R), Q.D(R+k+1, j-l-1)); // ECPCS 16A
        }),
        NUPACK_WHERE(j-k == s.last_nick() && i+l+1 < s.first_nick()) & sand(i+l+1, j-k, [=, &Q, &s](auto R) {
            return A.dot(Q.coax(j, i, R, s), Q.D(i+l+1, R), Q.Q(R+1, j-k-1)); // ECPCS 10A
        }),
        NUPACK_WHERE(j-k == s.last_nick() && i+l+1 < s.first_nick())
            & A.product(Q.coax(j, i, j-k-m-1, s), Q.D(i+l+1, j-k-m-1)), // ECPCS 10C
        NUPACK_WHERE(i+m+1 == s.first_nick() && j-l-1 >= s.last_nick())
            & A.product(Q.coax(i+m+k+1, j-1, i, s), Q.D(i+m+k+1, j-l-1)) // ECPCS 16C
    );
};

static auto S = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return A.dot(Q.CD(i, cspan(i, j+1)));
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return A.dot(Q.CD(i, cspan(s.last_nick(), j+1)) /* * dangle(i, Span, Span + 1)  */);
    }
);

static auto B = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto &&p) {
        cspan R(i+5, j-5);
        auto sum = [&] {return A.sum(
            B_single(i, j, SingleStrand(), A, Q, s, t),
            B_inextensible(i, j, SingleStrand(), A, Q, s, t),
            NUPACK_WHERE(i+11 <= j && t.can_close(s[i], s[j])) & A.product(A.sum(

                A.product(A.dot(Q.coax(j, i, R, s), Q.D(i+1, R), Q.M(R+1, j-1)), t.multi122), // multiloop_init_PAIR_COAXIAL_STACKING
                A.product(A.dot(Q.coax(R+1, j-1, i, s), Q.M(i+1, R), Q.D(R+1, j-1)), t.multi122), // multiloop_init_PAIR_COAXIAL_STACKING
                dangle_sum(A, [=, &Q, &s, &t](auto k, auto l) {
                    // cspan D(i+k+5, j-l-5);
                    // return NUPACK_WHERE(i+k+1 <= j-1) & A.product(A.sum(
                    cspan E(i+k+1, j-l-9);
                    return NUPACK_WHERE(i+k+l+11 <= j) & A.product(A.sum(
                        // A.product(A.dot(t.multi3(R+l-i-1), Q.MCS(R, j-l-1)), t.dangle(j-l, j, i, i+k, s)),
                        A.product(A.dot(t.multi3(E+l-i-1), Q.MCS(E, j-l-1)), t.dangle(j-l, j, i, i+k, s)),
                        A.product(A.dot(Q.M(i+k+1, E+4), Q.MS(E+5, j-l-1)), t.multi3(k+l), t.dangle(j-l, j, i, i+k, s))
                    //     A.product(A.dot(Q.M(i+k+1, R+5), Q.MS(R+6, j-l-1)), t.multi3(k+l), t.dangle(j-l, j, i, i+k, s)),
                    ), t.multi12);
                })
            ), t.terminal(s[j], s[i])),
            B_extensible(i, j, SingleStrand(), A, Q, s, t)
        );};
        return p(i, j, t.can_pair(false, next(s, i), next(s, j)), A, Q, s, t, sum);
    },
    [](int i, int j, MultiStrand, auto A, auto const &Q, auto const &s, auto const &t, auto &&p) {
        auto sum = [&] {return A.sum(
            B_inextensible(i, j, MultiStrand(), A, Q, s, t),
            // B_single is taken care of by B_cpd
            NUPACK_WHERE(t.can_close(s[i], s[j])) & A.product(A.sum(
                B_cps(i, j, A, Q, s, t),
                B_cpd(i, j, A, Q, s, t)
            ), t.terminal(s[j], s[i])),
            B_extensible(i, j, MultiStrand(), A, Q, s, t)
        );};
        return p(i, j, t.can_pair(true, next(s, i), next(s, j)), A, Q, s, t, sum);
    }
);


/// Q recursion for single and multiple strands
static auto Q = overload(
    [](int i, int j, SingleStrand, auto A, auto const &Q, auto const &s, auto const &t, auto const &) {
        return A.sum(
            t.one(),
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

}}}

#undef NUPACK_WHERE

