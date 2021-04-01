#pragma once
#include "../reflect/Print.h"
#include "../reflect/Memory.h"
#include "../standard/Ptr.h"
#include "../algorithms/Utility.h"
#include "../algorithms/TypeSupport.h"

namespace nupack {

/******************************************************************************************/

template <class T, std::size_t N, bool Stack=true>
class Storage {
    using Array = typename std::aligned_storage<sizeof(T), alignof(T)>::type[N];
    if_t<Stack, StackPtr<Array>, HeapPtr<Array>> data;

public:

    auto begin() {return reinterpret_cast<T *>(*data);}
    auto begin() const {return reinterpret_cast<T const *>(*data);}

    auto end() {return reinterpret_cast<T *>(*data + N);}
    auto end() const {return reinterpret_cast<T const *>(*data + N);}

    T & operator[](std::size_t i) {return begin()[i];}
    T const & operator[](std::size_t i) const {return begin()[i];}

    template <class ...Ts> void emplace(std::size_t i, Ts &&...ts) {new(begin() + i) T(fw<Ts>(ts)...);}
    template <class ...Ts> void destroy(std::size_t i) {begin()[i].~T();}
};

/******************************************************************************************/

/// FIFO queue of a compile-time size
template <class T, std::size_t N, bool Stack=true> class StaticFIFO {

    Storage<T, N, Stack> data;
    std::size_t b = 0, e = 0;

public:

    auto begin() const {return data.begin() + b;}
    auto end() const {return data.begin() + e;}

    template <class ...Ts> void emplace(Ts &&...ts) {
        if (likely(e < N)) {data.emplace(e++, fw<Ts>(ts)...); return;}
        throw std::out_of_range("StaticFIFO::emplace out of bounds");
    }

    bool empty() const {return b == e;}
    bool Full() const {return e == N;}
    void clear() {for (auto i = b; i != e; ++b) data.destroy(i); b = e = 0;}

    T & top() {
        if (unlikely(b >= e)) throw std::out_of_range("StaticFIFO::top out of bounds");
        return data[b];
    }

    void pop() {data.destroy(b++); if (unlikely(b == e)) b = e = 0;}

    StaticFIFO() = default;

    StaticFIFO(StaticFIFO const &o) : b(o.b), e(o.e) {
        std::uninitialized_copy(o.begin(), o.begin(), begin());
    }

    StaticFIFO(StaticFIFO &&o) noexcept : b(o.b), e(o.e) {
        for (auto i = b; i != e; ++b) data.emplace(i, std::move(o.data[i]));
    }

    StaticFIFO & operator=(StaticFIFO const &o) {
        clear(); b = o.b; e = o.e;
        std::uninitialized_copy(o.begin() + b, o.begin() + e, begin());
        return *this;
    }

    StaticFIFO & operator=(StaticFIFO &&o) noexcept {
        clear(); b = o.b; e = o.e;
        for (auto i = b; i != e; ++b) data.emplace(i, std::move(*o[i]));
        return *this;
    }

    ~StaticFIFO() {clear();}
};

/******************************************************************************************/

template <class T, std::size_t N, bool B> struct memory::impl<StaticFIFO<T, N, B>, void> {
    auto operator()(StaticFIFO<T, N> const &v) const {
        return N * sizeof(T) + sum(v, [](auto const &i) {return impl<T>(i) - sizeof(T);});
    }

    void erase(T &t) const {t.clear();}
};

/******************************************************************************************/

}
