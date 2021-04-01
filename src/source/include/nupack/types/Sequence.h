/** \file Sequence.h
 * @brief Contains base definitions, enums, Sequence class, SubsequenceList class to hold Sequences
 */
#pragma once
#include <string>
#include <string_view>

#include "Base.h"
#include "../iteration/Search.h"
#include "../iteration/View.h"
#include "../iteration/Patterns.h"

namespace nupack {

/******************************************************************************************/

struct SingleStrand : False {};
struct MultiStrand : True {};

/******************************************************************************************/

template <class V, NUPACK_IF(is_same<value_type_of<V>, Base>)>
string make_string(V const &v) {
    string s;
    s.reserve(len(s));
    for (auto b : v) s.push_back(b.letter());
    return s;
}

/******************************************************************************************/

/// BaseIter was previously Sequence::const_iterator, but this is now elided since they're essentially equivalent
using BaseIter = Base const *;
static_assert(is_dumb_ptr<BaseIter>, "");

using Subsequence = View<Base const *>;

/******************************************************************************************/

struct Sequence : public small_vec<Base, 32> {
    using base_type = small_vec<Base, 32>;
    Sequence() = default;

    using base_type::base_type;

    template <class V, NUPACK_IF(is_same<value_type_of<V>, Base>)>
    explicit Sequence(V const &v) : base_type(begin_of(v), end_of(v)) {}

    explicit Sequence(string_view letters);

    Sequence(iseq n, char b) : base_type(n, Base(b)) {}

    friend std::ostream & operator<<(std::ostream &os, Sequence const &s) {
        for (auto c : s) os << c.letter();
        return os;
    }

    string str() const {return make_string(static_cast<base_type const &>(*this));}
    explicit operator string() const {return str();}

    auto save_repr() const {return str();}
    void load_repr(string const &s) {*this = Sequence(s);}

    operator Subsequence() const & {return {base_type::data(), base_type::data() + base_type::size()};}
};

/******************************************************************************************/

using SubsequenceList = small_vec<Subsequence, 8>;
using SequenceList = small_vec<Sequence, 4>;

/******************************************************************************************/

template <class Out=small_vec<bool>, class S=Sequence>
auto one_hot_sequence(S const &sequence) {
    Out out(CanonicalBases.size() * len(sequence));
    auto it = begin_of(out);
    for (auto &&c : sequence) for (auto b : CanonicalBases) it++ == b;
    return out;
}

/******************************************************************************************/

NUPACK_UNARY_FUNCTOR(all_determined, all_of(Sequence(t), is_determined));
NUPACK_UNARY_FUNCTOR(has_wildcard, any_of(Sequence(t), is_wildcard));

NUPACK_BINARY_FUNCTOR(is_sequence_specialization, equal_ranges(t, u, is_base_specialization));

/******************************************************************************************/

// Basically the same as Sequence but cannot contain null bases or wild cards
struct Strand : Sequence {
    using base_type = Sequence;

    Strand() = default;

    // Forwarded constructor and then remove null bases
    template <class ...Ts, std::enable_if_t<can_construct<Sequence, Ts &&...>, int> = 0>
    Strand(Ts &&...ts) : Sequence(std::forward<Ts>(ts)...) {
        this->erase(std::remove(this->begin(), this->end(), Base('_')), this->end());
        NUPACK_ASSERT(!has_wildcard(*this), Sequence(*this), "Strand may not contain wildcards");
    }
};

using StrandList = small_vec<Strand, 4>;

/******************************************************************************************/

template <class RNG=decltype(StaticRNG) &>
Sequence sample(Sequence s, RNG &&rng=StaticRNG) {
    for (auto &b : s) b = b.sample(rng);
    return s;
}

Sequence reverse_complement(Sequence seq) noexcept;
Sequence reverse_wobble_complement(Sequence seq) noexcept;

inline bool is_palindromic(Sequence const &seq) {
    return seq == reverse_complement(seq);
}

template <class RNG=decltype(StaticRNG) &>
Sequence random_sequence(iseq n, real gc=0.5, RNG &&gen=StaticRNG) {
    auto dist = Base::distribution(gc);
    Sequence out; out.reserve(n);
    for (iseq i=0; i != n; ++i) out.push_back(Base::from_index(dist(gen)));
    return out;
}

template <class Out=StrandList, class RNG=decltype(StaticRNG) &>
Out random_sequences(iseq m, iseq n, real gc=0.5, RNG &&gen=StaticRNG) {
    auto dist = Base::distribution(gc);
    Out out(m);
    for (auto o = begin_of(out); o != end_of(out); ++o) {
        o->reserve(n);
        for (iseq i=0; i != n; ++i) o->push_back(Base::from_index(dist(gen)));
    }
    return out;
}

/******************************************************************************************/

template <> struct memory::impl<Sequence> {
    auto operator()(Sequence const &s) const {return len(s) * sizeof(Sequence::value_type);}
    void erase(Sequence &s) {Sequence s_; swap(s, s_);}
};

template <> struct memory::impl<Strand> : memory::impl<Sequence> {};

/******************************************************************************************/

using Nick = int;
static constexpr Nick const NoNick = -1;

/// Find sequence index of nick: the index is to the sequence after the nick
template <class V>
Nick find_nick(V const &v) {
    auto b = std::begin(v);
    for (Nick i = 0; b != std::end(v); ++b, ++i)
        if (front(*b) == Base('_')) return i;
    return NoNick;
}

/******************************************************************************************/

/// Convert multiple strings e.g. ["ACTGTA", "ACTGAT"] into a vector of strands
template <class V=SequenceList, class S=vec<string>, NUPACK_IF(!can_construct<string_view, S>)>
V to_sequences(S const &strs) {
    return vmap<V>(strs, [](auto const &s) {return value_type_of<V>(s);});
}

/// Convert a single string e.g. "ACTGTA+ACTGAT" into a vector of strings
vec<string> split_sequence_string(string_view s);

/// Convert a single string e.g. "ACTGTA+ACTGAT" into a vector of strands
template <class V=SequenceList>
V to_sequences(string_view s) {
    return vmap<V>(split_sequence_string(s), [](auto const &s) {return value_type_of<V>(s);});
}

/******************************************************************************************/

/**
 * @brief Given a vector of views, split it into 3 at bb of sequence b and ee of sequence e
 * @return std::pair<V, V> representing {1 + 3, 2}
 */
template <class V>
std::pair<V, V> split_midway(V const &v, const_iterator_of<V> b, const_iterator_of<V> e, BaseIter bb, BaseIter ee) {
    NUPACK_ASSERT(e >= b, e - b);
    // v1 should get len(v) - e + b + 1
    V v1; v1.assign(begin_of(v), b + 1);
    v1.back().set_end(bb + 1);
    extend(v1, e, end_of(v));
    v1[b + 1 - begin_of(v)].set_begin(ee);
    // v2 should get e + 1 - b sequences
    V v2; v2.assign(b, e + 1);
    v2.front().set_begin(bb);
    v2.back().set_end(ee + 1);
    return std::make_pair(v1, v2);
}

/******************************************************************************************/

/// Return new loop sequences after a base pair deletion
template <class V> V merged_seqs(V const &p, V const &k, size_type_of<V> pk, size_type_of<V> kp) {
    auto kpm = (!kp ? len(p) : kp) - 1;
    auto pkm = (!pk ? len(k) : pk) - 1;

    V ret; ret.reserve(len(p) + len(k) - 2);
    circular_cat(ret, p, begin_of(p) + kp, begin_of(p) + kpm);
    ret.front().set_begin(begin_of(k[pkm]));
    if (len(k) > 1) {ret.emplace_back(k[pk]); back(ret).set_begin(begin_of(p[kpm]));}
    else ret.front().set_begin(begin_of(p[kpm]));
    circular_cat(ret, k, begin_of(k) + pk + 1, begin_of(k) + pkm);
    return ret;
}

/******************************************************************************************/

/// Return new loop sequences after a dissociation event
template <class V> std::pair<V, V> get_split_seqs(V const &pseqs, V const &kseqs,
    size_type_of<V> pnick, size_type_of<V> knick, size_type_of<V> pk, size_type_of<V> kp) {
    V new_pseqs, new_kseqs;

    circular_cat(new_pseqs, pseqs, begin_of(pseqs) + pnick, begin_of(pseqs) + kp);
    new_pseqs.back().set_end(end_of(kseqs[pk]));
    circular_cat(new_pseqs, kseqs, begin_of(kseqs) + pk + 1, begin_of(kseqs) + knick);

    circular_cat(new_kseqs, kseqs, begin_of(kseqs) + knick, begin_of(kseqs) + pk);
    new_kseqs.back().set_end(end_of(pseqs[kp]));
    circular_cat(new_kseqs, pseqs, begin_of(pseqs) + kp + 1, begin_of(pseqs) + pnick);

    return std::make_pair(std::move(new_pseqs), std::move(new_kseqs));
}

/******************************************************************************************/

using Edge = int;
static constexpr Edge Ether = -1;
using EdgeList = small_vec<Edge>;

/******************************************************************************************/

}

namespace std {
    template <> struct hash<nupack::Sequence> {
        size_t operator()(nupack::Sequence const &s) const;
    };

    template <> struct hash<nupack::Strand> {
        size_t operator()(nupack::Strand const &s) const;
    };
}
