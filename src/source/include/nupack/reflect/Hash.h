#pragma once
#include <boost/functional/hash.hpp>

namespace nupack {

template <class T, class=void>
struct hash;

NUPACK_DETECT(has_hash, decltype(hash<decay<T>>()(declref<decay<T> const>())));
NUPACK_DETECT(has_std_hash, decltype(std::hash<T>()(declref<T const>())));
NUPACK_DETECT(has_member_hash, decltype(declref<T const>().hash()));

/******************************************************************************************/

template <class T>
struct hash<T, void_if<has_std_hash<T>>> : std::hash<T> {};

template <class T>
struct hash<T, void_if<has_member_hash<T>>> {
    constexpr std::size_t operator()(T const &t) const {return t.hash();}
};

/******************************************************************************************/

/// Functor to combine multiple hashes
struct HashCombiner {
    template <class ...Ts>
    std::size_t operator()(std::size_t t, Ts const &...ts) const {
        NUPACK_UNPACK(boost::hash_combine(t, ts));
        return t;
    }
};

/// Functor to get hash of object or multiple objects
struct Hasher {
    template <class T>
    auto operator()(T const &t) const -> decltype(hash<T>()(t)) {return hash<T>()(t);}

    template <class ...Ts, NUPACK_IF(sizeof...(Ts) >= 2 && all_of_c<has_hash<Ts>...>)>
    auto operator()(Ts const &...ts) const {return HashCombiner()(hash<Ts>()(ts)...);}
};

static constexpr auto combine_hashes = HashCombiner();
static constexpr auto hash_of = Hasher();

/******************************************************************************************/

template <class T> struct RangeHash {
    constexpr std::size_t operator()(T const &t) const {
        std::size_t seed = 0;
        for (auto const &i : t) boost::hash_combine(seed, hash_of(i));
        return seed;
    };
};

template <class T>
std::size_t range_hash(T const &t) {return RangeHash<T>()(t);}

/******************************************************************************************/

template <class ...Ts>
struct hash<std::tuple<Ts...>, void_if<(all_of_c<has_hash<Ts>...>)>> {
    constexpr std::size_t operator()(std::tuple<Ts...> const &t) const {return unpack(t, hash_of);}
};

template <class T, class U>
struct hash<std::pair<T, U>, void_if<(has_hash<T> && has_hash<U>)>> {
    constexpr std::size_t operator()(std::pair<T, U> const &t) const {return hash_of(t.first, t.second);}
};

/******************************************************************************************/

struct MemberHash {
    template <class T>
    std::size_t operator()(T const &t) {return unpack(members_of(t), hash_of);}
};

/******************************************************************************************/

}
