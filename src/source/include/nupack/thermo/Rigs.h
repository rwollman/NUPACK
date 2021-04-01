/**
 * @brief Algebraic rigs used in NUPACK dynamic programs
 *
 * @file Rigs.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Overflow.h"

namespace nupack { namespace thermo {

/******************************************************************************************/

struct first_arg {
    template <class T, class U>
    constexpr auto operator()(T &&t, U) const {
        static_assert(!is_same<decay<T>, Zero>, "");
        static_assert(is_same<decay<U>, Zero>, "");
        return fw<T>(t);}
};

/******************************************************************************************/

/// In the PF ring, + is plus() and * is times()
struct PF {
    using logarithmic = False;
    static constexpr auto zero() {return *::nupack::zero;}
    static constexpr auto one() {return *::nupack::one;}
    static constexpr auto plus() {return simd::plus;}
    static constexpr auto plus_eq() {return ::nupack::plus_eq;}
    static constexpr auto times() {return simd::times;}
    static constexpr auto invert() {return simd::invert;}
    static constexpr auto sum() {return simd::sum;}
    static constexpr auto ldexp() {return simd::ldexp;}

    // Return if the value overflows,
    // If so, set it to 0 (any finite number OK) so conversion to overflow doesn't break on it
    template <class M>
    static bool prevent_overflow(M &m) {
        // NUPACK_ASSERT(!(m < 0), m); // nan can happen if inf is an intermediate (0 * inf)
        if (unlikely(m < 0 || std::isnan(m) || m == std::numeric_limits<M>::infinity())) {m = 0; return true;}
        else return false;
    }

    // set the mantissa of the output using the initial exponent guess e
    template <class E, class F>
    static auto element_value(bool &err, F &&rule, E e0) {
        if constexpr(std::is_integral_v<E>) {
            for (uint i = 0; i != 512; ++i) {
                auto m = mantissa(rule, -e0);
                E e = exponent(rule, -e0);
                if (!unlikely(prevent_overflow(m))) {
                    auto p = simd::ifrexp(m);
                    if (p.second > 0) {
                        m = p.first;
                        e += p.second;
                    }
                    return std::make_pair(m, e + e0);
                } else {
                    print("nupack: adjusting dynamic program exponent", i, m, e, e0);
                    e0 += 32;//std::numeric_limits<M>::max_exponent / 4; // about half the total exponent range
                }
            }
            err = true;
            return decltype(std::make_pair(mantissa(rule, -e0), E()))();
        } else {
            auto m = simd::ldexp(mantissa(rule, -e0), exponent(rule, -e0));
            err = prevent_overflow(m);
            NUPACK_DASSERT(std::isfinite(m), m);
            return m;
        }
    }
};

// There is an abandoned branch implementing log sum exp algebra
// but even when optimized the code was around 10x slower than PF
struct LSE {
    using logarithmic = True;
    static constexpr auto zero() {return *::nupack::minf;}      // log(0)
    static constexpr auto one() {return *::nupack::zero;}       // log(1)
    static constexpr auto plus() {return simd::lse2;}           // log(exp(a) + exp(b))
    static constexpr auto plus_eq() {return ::nupack::plus_eq;}
    static constexpr auto times() {return simd::plus;}          // log(exp(a) * exp(b)) = a + b
    static constexpr auto invert() {return simd::unary_minus;}  // log(1/exp(a)) = -a
    static constexpr auto sum() {return simd::lse2;}             // log(sum(exp(a[:])))
    static constexpr auto ldexp() {return first_arg();}

    static bool prevent_overflow(Ignore) {return false;}

    template <class M, class E, class F>
    static auto set_element(M &m, E e, F &&rule) {return m = mantissa(fw<F>(rule), -e), False();}
};

/// In the MFE ring, + is min() and * is plus()
struct MFE {

    using logarithmic = True;
    static constexpr auto zero() {return *::nupack::inf;}
    static constexpr auto one() {return *::nupack::zero;}
    static constexpr auto plus() {return simd::min;}
    static constexpr auto plus_eq() {return simd::min_eq;}
    static constexpr auto times() {return simd::plus;}
    static constexpr auto invert() {return simd::unary_minus;}
    static constexpr auto sum() {return simd::minimum;}
    static constexpr auto ldexp() {return first_arg();}

    static bool prevent_overflow(Ignore) {return false;}

    template <class E, class F>
    static auto element_value(bool const &err, F &&rule, E e) {return mantissa(fw<F>(rule), -e);}
};

/******************************************************************************************/

}}
