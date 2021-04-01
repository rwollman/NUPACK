#pragma once

#include "../common/Config.h"
#include <algorithm>

namespace nupack {

/******************************************************************************************/

/// Heap that cannot have its size changed once constructed
template <class T, template <class...> class V=vec, class Comp=less_t>
class StaticHeap {
protected:
    V<T> c;

public:
    using container_type = V<T>;
    using value_type = T;

    template <class ...Ts> explicit StaticHeap(Ts &&...ts) : c(fw<Ts>(ts)...) {
        std::make_heap(begin_of(c), end_of(c), Comp());
    };

    void pop() {
        std::pop_heap(begin_of(c), end_of(c), Comp());
        c.pop_back();
    }

    container_type const & contents() const {return c;}

    value_type const & top() const {return c.front();}

    container_type sorted() && {
        std::sort_heap(begin_of(c), end_of(c), Comp());
        return std::move(c);
    }

    container_type sorted() const & {return StaticHeap(*this).sorted();}

    auto const & comparator() const {return Comp();}

    auto size() const {return len(c);}
};

NUPACK_DEFINE_TEMPLATE(is_static_heap, StaticHeap, class, template <class ...> class, class);

/******************************************************************************************/

/// Heap that can have its size changed once constructed
template <class T, template <class...> class V=vec, class Comp=less_t>
class Heap : public StaticHeap<T, V, Comp> {
    using StaticHeap<T, V, Comp>::c;

public:
    using value_type = T;
    using container_type = V<T>;

    template <class ...Ts> void emplace(Ts &&...ts) {
        c.emplace_back(fw<Ts>(ts)...);
        std::push_heap(begin_of(c), end_of(c), Comp());
    };

    void clear() {c.clear();}
};

NUPACK_DEFINE_TEMPLATE(isHeap, Heap, class, template <class ...> class, class);

/******************************************************************************************/

/// Heap that can have its size changed up to a maximum size
template <class T, template <class...> class V=vec, class Comp=less_t>
class MaxSizeHeap : public StaticHeap<T, V, Comp> {
    usize max_;
    using StaticHeap<T, V, Comp>::c;

public:
    using value_type = T;
    using container_type = V<T>;
    using StaticHeap<T, V, Comp>::top;

    template <class ...Ts>
    explicit MaxSizeHeap(usize m, Ts &&...ts) : StaticHeap<T, V, Comp>(fw<Ts>(ts)...), max_(m){};

    template <class ...Ts> void emplace_if(Ts &&...ts) {
        if (!max_) return;
        else if (this->size() < max_) {
            c.emplace_back(fw<Ts>(ts)...);
            std::push_heap(begin_of(c), end_of(c), Comp());
        } else {
            value_type t{fw<Ts>(ts)...};
            if (Comp()(t, top())) {
                std::pop_heap(begin_of(c), end_of(c), Comp());
                c.back() = std::move(t);
                std::push_heap(begin_of(c), end_of(c), Comp());
            }
        }
    }

    auto max_elements() const {return max_;}

    void clear() {c.clear();}
};

NUPACK_DEFINE_TEMPLATE(isMaxSizeHeap, MaxSizeHeap, class, template <class...> class, class);

/******************************************************************************************/

}
