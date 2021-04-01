/** \file Sequence.h
 * @brief Contains Base class for nucleotides
 */
#pragma once
#include <iostream>

#include "../common/Random.h"
#include "../standard/Array.h"
#include "Matrix.h"
#include "../iteration/Search.h"

namespace nupack {

/******************************************************************************************/

using BaseIndex = std::uint_fast8_t;

struct Base {
    static constexpr const std::array<char, 16> names{
        {'A', 'C', 'G', 'T', 'R', 'M', 'S', 'W', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '_'}
    };

    static constexpr std::array<std::array<bool, 4>, 16> masks = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1},
        {1, 0, 1, 0}, // AG
        {1, 1, 0, 0}, // AC
        {0, 1, 1, 0}, // CG
        {1, 0, 0, 1}, // AU
        {0, 0, 1, 1}, // GU
        {0, 1, 0, 1}, // CU
        {1, 1, 1, 0}, // ACG
        {1, 1, 0, 1}, // ACU
        {1, 0, 1, 1}, // AGU
        {0, 1, 1, 1}, // CGU
        {1, 1, 1, 1}, // ACGU
        {0, 0, 0, 0}  // none
    }};

    static constexpr std::array<BaseIndex, 16> complements {
        3, 2, 1, 0, 9, 8, 6, 7, 5, 4, 13, 12, 11, 10, 14, 15
    };

    static constexpr std::array<BaseIndex, 16> wobble_complements {
        3, 2, 9, 4, 9, 8, 13, 12, 14, 4, 13, 12, 14, 14, 14, 15
    };

    static std::array<std::discrete_distribution<BaseIndex>, 15> const distributions;

    /// Look up the encoded base value in a non-constexpr way - better printout
    static constexpr BaseIndex lookup(char letter) {
        switch(letter) {
            case('A'): return 0;   case('a'): return 0;
            case('C'): return 1;   case('c'): return 1;
            case('G'): return 2;   case('g'): return 2;
            case('T'): return 3;   case('t'): return 3;
            case('U'): return 3;   case('u'): return 3;
            case('R'): return 4;   case('r'): return 4;  // AG
            case('M'): return 5;   case('m'): return 5;  // AC
            case('S'): return 6;   case('s'): return 6;  // CG
            case('W'): return 7;   case('w'): return 7;  // AU
            case('K'): return 8;   case('k'): return 8;  // GU
            case('Y'): return 9;   case('y'): return 9;  // CU
            case('V'): return 10;  case('v'): return 10; // ACG
            case('H'): return 11;  case('h'): return 11; // ACU
            case('D'): return 12;  case('d'): return 12; // AGU
            case('B'): return 13;  case('b'): return 13; // CGU
            case('N'): return 14;  case('n'): return 14; // ACGU
            case('_'): return 15;
            default: NUPACK_ERROR("invalid letter for nucleotide encountered", int(letter), letter);
        }
    }

    BaseIndex value;

    constexpr explicit Base(char letter) : value(lookup(letter)) {};
    constexpr explicit Base(bool, BaseIndex t) : value(t) {}

    Base() = default;

    constexpr operator BaseIndex () const {return value;}
    constexpr operator BaseIndex & () {return value;}

    constexpr bool operator==(Base rhs) const {return value == rhs.value;}
    constexpr bool operator!=(Base rhs) const {return value != rhs.value;}
    constexpr bool operator< (Base rhs) const {return value <  rhs.value;}
    constexpr bool operator> (Base rhs) const {return value >  rhs.value;}
    constexpr bool operator<=(Base rhs) const {return value <= rhs.value;}
    constexpr bool operator>=(Base rhs) const {return value >= rhs.value;}

    struct from_index_t {
        constexpr Base operator()(BaseIndex value) const {
            return (value < 16) ? Base(true, value) : throw Error("Invalid nucleotide index");
        }
    };

    static constexpr auto const from_index = from_index_t{};

    char safe_letter() const {
        try {return Base::names.at(value);}
        catch (...) {NUPACK_BUG("Bad index for nucleotide lookup", int(value), value);}
    }

    template <bool B=true, NUPACK_IF(B && DebugBounds)>
    char letter() const {return safe_letter();}

    template <bool B=true, NUPACK_IF(B && !DebugBounds)>
    constexpr char letter() const noexcept {return Base::names[value];}

    constexpr auto const & mask() const {return Base::masks.at(value);}

    static auto distribution(real A, real C, real G, real U) {
        return discrete_distribution<BaseIndex>(std::array<real, 4>{{A, C, G, U}});
    }

    static auto distribution(real gc=0.5) {return distribution((1-gc)/2, gc/2, gc/2, (1-gc)/2);}

    template <class RNG=decltype(StaticRNG) &>
    Base sample(RNG &&rng=StaticRNG) const {
        if (value < 4 || *this == Base('_')) return *this;
        else return from_index(distributions[value](rng));
    }

    friend constexpr Base complement(Base b) {return from_index(Base::complements[b.value]);}
    friend constexpr Base wobble_complement(Base b) {return from_index(Base::wobble_complements[b.value]);}

    friend std::ostream & operator<<(std::ostream &os, Base c) {return os << c.letter();}

    auto save_repr() const {return value;}
    void load_repr(BaseIndex c) {value = c;}

};

/// Functor for Base or char to say if it is wildcard nucleotide
NUPACK_UNARY_FUNCTOR(is_wildcard, (3 < Base(t).value && t != Base('_')));
/// Functor for Base or char to say if it is not a wildcard nucleotide
NUPACK_UNARY_FUNCTOR(is_determined, (Base(t).value < 4 || t == Base('_')));
/// Functor for Base or char to say if it is one of ACGTU
NUPACK_UNARY_FUNCTOR(is_canonical, (Base(t).value < 4));

static constexpr std::array<Base, 4> CanonicalBases{{Base('A'), Base('C'), Base('G'), Base('T')}};

// numpy.array([[s @ g == s.sum() for s in masks] for g in masks]).astype(int)
static constexpr std::array<std::array<bool, 16>, 16> Specializations = {{
    {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1}},
    {{0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}},
    {{0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1}},
    {{1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1}},
    {{1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1}},
    {{1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1}},
    {{0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1}},
    {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}
}};

/// return if the second argument is a specialization of the first
NUPACK_BINARY_FUNCTOR(is_base_specialization, Specializations[t][u]);

template <> struct traits::is_character_t<Base> : True {};
static_assert(is_character<Base>, "");

static_assert(sizeof(Base) == sizeof(BaseIndex), "Should be same size");
static_assert(std::is_pod<Base>::value, "Should be POD type");

using Base_Pair = std::pair<iseq, iseq>;
using Base_Count = std::array<iseq, 4>;

template <class T> using BaseMat = typename arma::Mat<T>::template fixed<4, 4>;
template <class T> using Base_Array = typename arma::Col<T>::template fixed<4>;

/******************************************************************************************/

inline constexpr bool is_gu(Base b, Base c) {
    return (b == Base('G') && c == Base('U')) || (b == Base('U') && c == Base('G'));
}

inline constexpr bool is_au(Base b, Base c) {
    return (b == Base('A') && c == Base('U')) || (b == Base('U') && c == Base('A'));
}

inline constexpr bool is_gc(Base b, Base c) {
    return (b == Base('G') && c == Base('C')) || (b == Base('C') && c == Base('G'));
}

enum class WobblePairing : bool {off=false , on=true};
enum class WobbleClosing : bool {off=false , on=true};

struct Pairable : MemberOrdered {
    bool wobble_pairing;
    bool wobble_closing;
    NUPACK_REFLECT(Pairable, wobble_pairing, wobble_closing);

    constexpr auto turn() const {return 3;}

    bool can_close(Base b, Base c) const {return (b + c == 3) || (wobble_closing && b + c == 5);}
    bool can_pair(Base b, Base c) const {return (b + c == 3) || (wobble_pairing && b + c == 5);}

    bool operator()(Base b, Base c) const {return can_pair(b, c);}

    template <class Iter>
    bool operator()(bool diff_strand, Iter b, Iter c) const {
        if (!diff_strand) NUPACK_REQUIRE(b, <=, c, "same-strand bases should be ordered for this function");
        return (diff_strand || b + turn() < c) && can_pair(*b, *c);
    }

    template <class V>
    bool check_loop(V const &v) const {
        bool ok = true;
        for_circularly_adjacent(v, [&](auto const &s1, auto const &s2) {
            Base const b = s1.back(), c = s2.front();
            ok = ok && ((b == Base('_') && c == Base('_')) || (*this)(b, c));
        });
        return ok;
    }
};

/******************************************************************************************/

}

namespace std {
    template <> struct hash<nupack::Base> {
        size_t operator()(nupack::Base b) const {return hash<typename nupack::BaseIndex>()(b);}
    };
}
