#pragma once

#include "FIFO.h"

#ifdef NUPACK_NO_COROUTINE
// namespace nupack {template <class T> struct Coroutine {using pull_type = int; using push_type = int;};}
namespace nupack {template <class T> struct Coroutine {using pull_type = std::function<void(T)>; using push_type = std::function<void(T)>;};}
#else
#include <boost/coroutine2/all.hpp>
namespace nupack {template <class T> using Coroutine = boost::coroutines2::coroutine<T>;}
#endif

namespace nupack {

/******************************************************************************************/

/// Generator from a function and arguments exposing an STL iterator
template <class T> class Generator {
    using Push = typename Coroutine<T>::push_type;

    /// Simple capture of function and arguments
    template <class ...Ts> class Caller {
        std::tuple<Ts...> ts;

        template <std::size_t ...Is>
        constexpr decltype(auto) call(indices_t<Is...>, Push &sink) const {
            return std::get<sizeof...(Ts) - 1>(ts)(std::get<Is>(ts)..., sink);
        }

    public:

        template <class ...Args> Caller(Args &&...args) : ts(fw<Args>(args)...) {}

        constexpr decltype(auto) operator()(Push &sink) const {
            return call(indices_up_to<sizeof...(Ts)-1>(), sink);
        }
    };

public:

    typename Coroutine<T>::pull_type source;

    using value_type = T;
    using simple_type = True;

    template <class F, class ...Ts, NUPACK_IF(!is_same<decay<F>, Generator>)>
    Generator(F &&f, Ts &&...ts) : source(Caller<no_ref<F>, no_ref<Ts>...>(fw<F>(f), fw<Ts>(ts)...)) {}

    T operator()() {T ret = source.get(); source(); return ret;}
    bool empty() const {return !bool(source);}

    decltype(auto) begin() {return begin_of(source);}
    decltype(auto) end() {return end_of(source);}
};

/******************************************************************************************/

NUPACK_DEFINE_TEMPLATE(is_generator, Generator, class);

/******************************************************************************************/

/// Generator that switches context periodically, not after each push-pull. Designed for cases
/// where the elements are not memory-intensive and the context-switching is expensive
template <class T, std::size_t N> class BlockedGenerator {

    static_assert(!is_ref<T>, R"(Use std::ref or similar for BlockedGenerator of reference
types, but remember that the returned objects might alias each other!)");

    using Buffer = StaticFIFO<T, N>;
    using Push = typename Coroutine<void>::push_type;

    static constexpr bool big = (N > 2);

    template <class ...Ts> class Caller {
        std::tuple<Ts...> ts;
        Buffer &buff;

        template <class F, std::size_t ...Is>
        constexpr decltype(auto) call(indices_t<Is...>, F &&f) const {
            return std::get<sizeof...(Ts) - 1>(ts)(std::get<Is>(ts)..., fw<F>(f));
        }

    public:

        template <class ...Args> Caller(Buffer &b, Args &&...args) : buff(b), ts(fw<Args>(args)...) {}

        constexpr decltype(auto) operator()(Push &sink) const {
            call(indices_up_to<sizeof...(Ts)-1>(), [&](auto &&t) {
                buff.emplace(fw<decltype(t)>(t));
                if (unlikely_if<big>(buff.Full())) {sink();}
            });
            if (likely_if<big>(!buff.empty())) sink();
        }
    };

public:

    Buffer buff;
    typename Coroutine<void>::pull_type source;

    using simple_type = True;
    using value_type = T;
    using reference_type = T &&;

    T operator()() {
        T t = std::move(buff.top());
        buff.pop();
        if (unlikely_if<big>(buff.empty()) && likely(source)) {source();}
        return t;
    }

    bool empty() const {return buff.empty() && !bool(source);}

    class iterator {
        struct proxy {T value; T operator*() {return std::move(value);}};

        BlockedGenerator *ptr;

    public:

        constexpr iterator(BlockedGenerator *c) : ptr(c) {}
        iterator & operator++() {
            ptr->buff.pop(); // discard current position
            if (likely_if<big>(!ptr->buff.empty())) return *this;
            if (likely(ptr->source)) ptr->source(); // move on to next if available
            if (ptr->buff.empty()) ptr = nullptr; // kill
            return *this;
        }
        proxy operator++(int) {proxy ret{std::move(ptr->buff.top())}; ++(*this); return ret;}
        T operator*() const {return std::move(ptr->buff.top());}
        bool operator==(iterator o) const {return ptr == o.ptr;}
        bool operator!=(iterator o) const {return ptr != o.ptr;}
    };

    constexpr auto end() {return iterator(nullptr);}
    auto begin() {return empty() ? end() : iterator(this);}

    template <class F, class ...Ts, NUPACK_IF(!is_same<decay<F>, BlockedGenerator>)>
    BlockedGenerator(F &&f, Ts &&...ts) : buff(), source(Caller<no_ref<F>, no_ref<Ts>...>(buff, fw<F>(f), fw<Ts>(ts)...)) {}
};

NUPACK_EXTEND_VARIADIC(is_generator, BlockedGenerator, class, std::size_t);

/******************************************************************************************/

}
