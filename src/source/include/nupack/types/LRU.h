/** \file LRU.h
 * @brief Implementation of least-recently-used cache, with specialization for a memory limit
 */
#pragma once

#include "../reflect/Memory.h"
#include "../common/Threading.h"
#include "../algorithms/Operators.h"

#include "../standard/List.h"
#include "../standard/Map.h"

namespace nupack {

/******************************************************************************************/

struct MemoryLimit {
    using length_type = std::size_t;
    length_type length=0, capacity;

    NUPACK_REFLECT(MemoryLimit, length, capacity);

    constexpr MemoryLimit(length_type m=0) : capacity(m) {}

    template <class T>
    void add(T const &t) {length += memory::measure(t);}

    template <class T>
    void remove(T const &t) {length = unsigned_minus(length, memory::measure(t));}

    bool ok() const {return length <= capacity;}

    bool satisfiable() const {return capacity > 0;}

    void clear() {length = 0;}
};

/******************************************************************************************/

/**
 * @brief Least-recently-used cache with partial thread-safety
 * @todo Improve the thread-safety and possibly performance of LRU()
 * @todo Maybe only hold the key once, e.g. hold reference to map element in list instead
 * @todo Fix the members() which is janky right now in excluding the mutex
 * @tparam K Key type
 * @tparam V Value type
 * @tparam L=memory::measure Length operation functor
 * @tparam List_=List List template
 * @tparam Map_=Map Map template
 */
template <class K, class V, class L=MemoryLimit, template <class...> class List_=List, template <class...> class Map_=HashMap>
struct LRU : Iterable<LRU<K, V, L, List_, Map_>> {
    using is_lru = True;
    using value_type = std::pair<K, V>;
    using limit_type = L;
    using mapped_type = V;
    using list_type = List_<value_type>;
    using iterator = iterator_of<list_type>;
    using const_iterator = const_iterator_of<list_type>;

    /// An atomic integer keeps track of any referrents to open Spot()s that haven't been filled
    using map_type = Map_<K, iterator>;
    using key_type = K;
    using size_type = size_type_of<map_type>;
    using difference_type = difference_type_of<map_type>;

    /**************************************************************************************/

    NUPACK_REFLECT(LRU, contents, map, limit);
    auto & iter() {return contents;}

    list_type contents;
    map_type map;
    L limit;
    SharedMutex mutable mut;

    LRU(L lim={}) : limit(std::move(lim)) {};
    LRU(LRU const &o) : contents(o.contents), limit(o.limit) {build_map();}
    LRU(LRU &&o) : contents(std::move(o.contents)), map(std::move(o.map)), limit(std::move(o.limit)) {}
    LRU & operator=(LRU o) {contents = std::move(o.contents); limit = std::move(o.limit); build_map(); return *this;}

    /**************************************************************************************/

    auto read_lock() const {return shared_lock(mut);}
    auto write_lock() const {return unique_lock(mut);}

    auto key_comp() const {return map.key_comp();}

    /**
     * @brief Try to erase the last non-referenced element in contents, then its associated element in map
     * @return true Element was erased
     * @return false Element was not erased (no non-referenced elements in the container)
     */
    bool pop_back() {
        if (contents.empty()) return false;
        auto const m = map.find(contents.back().first);
        NUPACK_DASSERT(m != std::end(map));
        map.erase(m);
        erase_list_iterator(std::prev(contents.end()));
        return true;
    }

    /// Try to erase least recently used element until length is under quota
    void shrink_to_fit() {while (!limit.ok() && pop_back()) {}}

    /// Empty all contents
    void clear() {
        limit.clear();;
        map.clear();
        contents.clear();
    }

    /// Not thread-safe
    void erase(key_type const &k) {
        auto const it = map.find(k);
        if (it->second != contents.end())
            erase_list_iterator(it->second);
        map.erase(it);
    }

    template <class ...Ts>
    bool try_emplace(key_type k, Ts &&...ts) {
        auto map_it = map.emplace(k, std::end(contents)).first;
        if (map_it->second == std::cend(contents)) {
            prepend(std::move(k), static_cast<Ts &&>(ts)...);
            map_it->second = std::begin(contents);
            shrink_to_fit();
            return true;
        }
        if (std::begin(contents) != map_it->second)
            contents.splice(std::begin(contents), contents, map_it->second);
        return false;
    }

    template <class ...Ts>
    bool insert_or_assign(key_type k, Ts &&...ts) {
        auto map_it = map.emplace(k, std::end(contents)).first;
        bool inserted = false;
        prepend(std::move(k), static_cast<Ts &&>(ts)...);
        if (map_it->second == std::cend(contents)) inserted = true;
        else erase_list_iterator(map_it->second);
        map_it->second = std::begin(contents);
        shrink_to_fit();
        return inserted;
    }

    /// Equivalent to std::map::find but moves the found element at the front of the cache (not thread-safe)
    iterator find_and_refresh(key_type const &k) {
        if (auto it = map.find(k); it != std::end(map)) {
            contents.splice(begin_of(contents), contents, it->second);
            return it->second;
        } else return std::end(contents);
    }

    /// Equivalent to std::map::find (thread-safe if LRU not being modified)
    const_iterator find(key_type const &k) const {
        if (auto it = map.find(k); it != std::end(map)) return it->second;
        else return std::end(contents);
    }

    /// Equivalent to std::map::find (thread-safe if LRU not being modified)
    iterator find(key_type const &k) {
        if (auto it = map.find(k); it != std::end(map)) return it->second;
        else return std::end(contents);
    }

    static constexpr auto repr_names() {return make_names("contents", "limit");}

    auto save_repr() const {return make_members(contents, limit);}

    /// Load the contents, then build the map from it
    void load_repr(decltype(contents) l, decltype(limit) o) {
        limit = std::move(o);
        contents = std::move(l);
        build_map();
    }

protected:
    /// Build map from existing contents (not generally for public consumption)
    void build_map() {
        for (auto i : iterators(contents)) {
            map.emplace(i->first, i);
            limit.add(*i);
        }
    }

    template <class ...Ts>
    void prepend(key_type &&k, Ts &&...ts) {
        contents.emplace_front(std::move(k), mapped_type{static_cast<Ts &&>(ts)...});
        limit.add(contents.front());
    }

    /// Erase element in contents and decrement the length counter
    void erase_list_iterator(iterator it) {
        limit.remove(*it);
        contents.erase(it);
    }
};

NUPACK_DETECT(is_lru, void_if<T::is_lru::value>);

/******************************************************************************************/

}
