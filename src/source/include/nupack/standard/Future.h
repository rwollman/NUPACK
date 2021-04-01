#pragma once
#include "Optional.h"
#include "../reflect/Memory.h"
#include "../reflect/Print.h"

#include <future>
#include <thread>
#include <csignal>
#include <atomic>

namespace nupack {

/******************************************************************************************/

NUPACK_DEFINE_TEMPLATE(is_future, std::future, class);

/******************************************************************************************/

template <class T>
struct memory::impl<std::future<T>, void> {
    std::size_t operator()(std::future<T> const &t) const {
        return sizeof(t);
    }
    void erase(std::future<T> &t) const {T t_; t.swap(t_);}
};

template <class T> struct io::Printer<std::future<T>, void> {
    void operator()(std::ostream &os, std::future<T> const &) const {
        os << "future(" << TypeName<T>() << ')';
    }
};

template <class D, class F, class ...Ts>
std::optional<std::invoke_result_t<F &&, Ts &&...>> call_with_timeout(D const &duration, F &&f, Ts &&...ts) {
    auto const sig = std::make_shared<std::atomic<int>>(0);
    auto r = std::async(std::launch::async, [&] {
        ThreadLocalSignal = sig;
        return std::invoke(static_cast<F &&>(f), static_cast<Ts &&>(ts)...);
    });

    if (r.wait_for(duration) == std::future_status::ready) {
        return std::move(r).get();
    } else {
        sig->store(SIGINT, std::memory_order_relaxed);
        r.wait();
        return std::nullopt;
    }
}

/******************************************************************************************/

}
