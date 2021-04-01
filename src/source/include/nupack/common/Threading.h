/**
 * @brief Definitions of mutex wrappers and locks
 *
 * @file Threading.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "Config.h"
#include "../algorithms/Constants.h"

#include <random>
#include <mutex>
#include <atomic>
#include <thread>

// shared_mutex is C++14 but some environments might not have it
#ifndef NUPACK_NO_SHARED_MUTEX
#   if __has_include(<shared_mutex>)
#       include <shared_mutex>
#   else
#       define NUPACK_NO_SHARED_MUTEX
#   endif
#endif

namespace nupack {

/******************************************************************************************/

/// Make current thread sleep for a duration
template <class T> void sleep(T const &t) {std::this_thread::sleep_for(t);}

/******************************************************************************************/

#ifdef NUPACK_NO_SHARED_MUTEX
    /// stand-in type for std::shared_timed_mutex
    using SharedMutex = std::mutex;
    inline auto shared_lock(SharedMutex &m) {return std::unique_lock<SharedMutex>(m);}
#else
    using SharedMutex = std::shared_timed_mutex;
    inline auto shared_lock(SharedMutex &m) {return std::shared_lock<SharedMutex>(m);}
#endif
inline auto unique_lock(SharedMutex &m) {return std::unique_lock<SharedMutex>(m);}

/******************************************************************************************/

/// Execute a functor with a lock in place
template <class Mutex, class F, class ...Ts>
decltype(auto) with_lock(Mutex &&t, F &&f, Ts &&...ts) {
    std::lock_guard<decay<Mutex>> lock(fw<Mutex>(t));
    return fw<F>(f)(fw<Ts>(ts)...);
}

/******************************************************************************************/

}
