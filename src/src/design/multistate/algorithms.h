#pragma once

#include <iterator>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <utility>

namespace nupack {

template<bool b, class T>
using enable_if_t = typename std::enable_if<b, T>::type;

template <bool B>
using int_if = enable_if_t<B, int>;

// http://stackoverflow.com/questions/17865488/find-the-nth-element-satisfying-a-condition
// n >= 1
template <typename Iterator, typename Pred, typename Counter>
Iterator find_nth( Iterator first, Iterator last, Pred closure, Counter n ) {
  typedef typename std::iterator_traits<Iterator>::reference Tref;
  return std::find_if(first, last, [&](Tref x) {
    return closure(x) && !(--n);
  });
}

template <typename Container, typename Pred, typename Counter>
auto find_nth(Container && c, Pred && p, Counter n) -> decltype(c.begin()) {
  return find_nth(c.begin(), c.end(), std::forward<Pred>(p), n);
}

template <class Hier>
bool smaller_depth(const Hier & a, const Hier & b) {
  return a.get_max_depth() < b.get_max_depth();
}

template <class Hier, int_if<std::is_pointer<Hier>::type> = 0>
bool smaller_depth(const Hier & a, const Hier & b) {
  return smaller_depth(*a, *b);
}

template <class Iterator, class Container>
bool contains_it(const Iterator & it, Container & container) {
  return it != container.end();
}

template <class Key, class Map>
bool contains(Key & key, Map & map) {
  return contains_it(map.find(key), map);
}

template <class Key, class Container>
bool lin_contains(Key & key, Container & container) {
  return contains_it(std::find(container.begin(), container.end(), key), container);
}

template <class Container, class Pred>
bool all_of(Container && c, Pred && p) {
  return std::all_of(c.begin(), c.end(), std::forward<Pred>(p));
}

template <class Container, class Value, class Func>
Value accumulate(Container && c, Value init, Func && f) {
  return std::accumulate(c.begin(), c.end(), init, std::forward<Func>(f));
}

template <class A, class B>
void append(A & a, const B & b) {
  a.insert(a.end(), b.begin(), b.end());
}

template <class A>
void sort(A & a) {
  std::sort(a.begin(), a.end());
}

template <class A, class F>
void for_both(std::pair<A, A> & a, F && f) {
  f(a.first);
  f(a.second);
}

}
