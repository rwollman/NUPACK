#pragma once

#include <shared.h>

#include <vector>
#include <unordered_map>

#include <type_traits>

namespace nupack {
  template <bool B, class T = void>
  using enable_if_t = typename std::enable_if<B, T>::type;
  
  template <bool B>
  using int_if = enable_if_t<B, int>;
  
  template<typename... Ts> 
  struct make_void { using type = void;};
  
  template<typename... Ts> 
  using void_t = typename make_void<Ts...>::type;
  
  template <typename T, typename = void>
  struct is_iterable : std::false_type {};
  
  template <typename T>
  struct is_iterable<T, void_t<decltype(std::declval<T>().begin()), 
      decltype(std::declval<T>().end())>>
      : std::true_type {};
  
  template <typename T, typename = void>
  struct is_pair : std::false_type {};
  
  template <typename T>
  struct is_pair<T, void_t<decltype(std::declval<T>().first),
      decltype(std::declval<T>().second)>>
      : std::true_type {};
  
    
  using vec_structure = std::vector<int>;
  
  using trinary = uint8_t;
  using AllowTable = std::vector<std::vector<trinary>>;
  
  struct VecHash {
    std::size_t operator()(const std::vector<int> & v) const {
      std::size_t seed = 0;
      const std::hash<int> hasher{};
      for (auto i : v) seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
      return seed;
    }
  };
  
  using Map = std::unordered_map<std::vector<int>, DBL_TYPE, VecHash>;
  
  template <class T>
  struct PairHash {
    std::size_t operator()(const std::pair<T, T> & v) const {
      std::size_t seed = 0;
      const std::hash<T> hasher{};
      seed ^= hasher(v.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
      seed ^= hasher(v.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
      return seed;
    }
  };
}
