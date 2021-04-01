#pragma once

#include "adapter.h"

#include <vector>
// #include <unordered_map>
#include <map>

#include <type_traits>

namespace nupack {
namespace custom_csp {

template <bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

template <bool B>
using int_if = enable_if_t<B, int>;

template<typename... Ts>
struct make_void { using type = void;};

template<typename... Ts>
using void_t = typename make_void<Ts...>::type;

template <typename T, typename = void>
struct is_Iterable_t : std::false_type {};

template <typename T>
struct is_Iterable_t<T, void_t<decltype(std::declval<T>().begin()),
       decltype(std::declval<T>().end())>>
       : std::true_type {};

template <typename T, typename = void>
struct is_pair_t : std::false_type {};

template <typename T>
struct is_pair_t<T, void_t<decltype(std::declval<T>().first),
       decltype(std::declval<T>().second)>>
       : std::true_type {};

template <class T> static constexpr bool is_iterable = is_Iterable_t<T>::value;
template <class T> static constexpr bool is_pair = is_pair_t<T>::value;


using vec_structure = vec<int>;

using trinary = uint8_t;
using AllowTable = vec<vec<trinary>>;

struct VecHash {
    std::size_t operator()(const vec<int> & v) const {
        std::size_t seed = 0;
        const std::hash<int> hasher{};
        for (auto i : v) seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

using WeightMap = std::map<vec<int>, real>;

template <class T>
struct PairHash {
    std::size_t operator()(const std::pair<T, T> & v) const {
        std::size_t seed = 0;
        const std::hash<T> hasher{};
        seed ^= hasher(v.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(v.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

}
}
