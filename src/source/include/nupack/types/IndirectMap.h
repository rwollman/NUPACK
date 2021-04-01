#pragma once
#include "../algorithms/Utility.h"

/******************************************************************************************/

namespace nupack {

/******************************************************************************************/

template <class Map, class F> struct IndirectMap {
    Map const &map;
    F f;

    using iterator = IndirectIter<decltype(begin_of(map)), F>;

    template <class ...Ts>
    IndirectMap(Map const &m, Ts &&...ts) : map(m), f(fw<Ts>(ts)...) {}

    auto begin() const {return iterator(begin_of(map), f);}
    auto end() const {return iterator(end_of(map), f);}

    auto size() const {return len(map);}

    template <class ...Ts> auto find(Ts &&...ts) const {return iterator(map.find(fw<Ts>(ts)...), f);}
};

/******************************************************************************************/

// Make something that looks like a set of keys from a map
template <class Map> IndirectMap<Map, at_t<0>> map_keys(Map const &m) {return {m};}
template <class Map> IndirectMap<Map, at_t<1>> map_values(Map const &m) {return {m};}

/******************************************************************************************/

// Make something that looks like a map from a set
template <class Set, class F>
auto map_from_set(Set const &s, F f) {return IndirectMap<Set, F>(s, f);}

/******************************************************************************************/

}
