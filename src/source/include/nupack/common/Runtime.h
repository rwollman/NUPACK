/**
 * @brief Declarations for environmental variable lookup, type name strings, demangling, signal handlers
 *
 * @file Runtime.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once
#include "../algorithms/TypeSupport.h"
#include "Config.h"

#include <mutex>
#include <chrono>
#include <typeinfo>
#include <sstream>
#include <iomanip>
#include <memory>
#include <typeindex>
#include <cctype>

namespace nupack {
using namespace std::chrono_literals;
namespace chrono = std::chrono;

/******************************************************************************************/

void print_backtrace(std::ostream &, std::size_t n=32);
string get_backtrace(std::size_t n=32);

/******************************************************************************************/

/// Prettify a C++ type name
string trim_type_name(string name, int n=30);
/// Demangle a C++ type name
string demangle(string s);

/// Same as type_name but cuts off any template arguments
string short_type_name(string s);

/// A string-convertible class holding the C++ type name for a given type
template <class T> struct TypeName {
    friend std::ostream & operator<<(std::ostream &os, TypeName) {return os << demangle(typeid(T).name());}
    operator string() const {return demangle(typeid(T).name());}
    string operator *() const {return *this;}

    friend auto operator+(TypeName t, char const *s) {return string(t) + string(s);}
    friend auto operator+(char const *s, TypeName t) {return string(s) + string(t);}

    friend auto operator+(TypeName t, string s) {return string(t) + s;}
    friend auto operator+(string s, TypeName t) {return s + string(t);}

    auto operator~() const {return short_type_name(*this);}

    template <class T_> constexpr auto operator()(T_ &&t) const {return TypeName<decltype(t)>();}
};

template <class T>
constexpr TypeName<T> type_name(T const &) {return {};}

/******************************************************************************************/

void time_sink_impl(void const *);

template <class T>
void time_sink(T const &t) {time_sink_impl(std::addressof(t));}

/// Return environment variable
string get_env(string key);

// Run 100% for a specified time
void spin_processor(chrono::duration<real> const t);

/******************************************************************************************/

/// Incomplete type used in SignalRuntime()
struct RuntimeContents;

/// Catch signal errors (SIGINT) and turn them into exceptions while the object is in scope
struct SignalRuntime {
    std::mutex mut;
    std::shared_ptr<RuntimeContents> contents;
    SignalRuntime(SignalRuntime const &) = delete;
    SignalRuntime(SignalRuntime &&) = delete;

    SignalRuntime(); //< Acquire signal handler for "which" -- default SIGINT
    SignalRuntime(int const *begin, int const *end);
    SignalRuntime(std::initializer_list<int> v) : SignalRuntime(v.begin(), v.end()) {}

    ~SignalRuntime(); //< Replace signal handlers with the old ones
};

/******************************************************************************************/

/// Return current time as a string
template <class Clock=chrono::system_clock>
string timestamp(chrono::time_point<Clock> const &t=Clock::now()) {
    auto tt = Clock::to_time_t(t);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&tt), "%Y-%m-%d-%H-%M-%S");
    return ss.str();
}

/******************************************************************************************/

/// Path exists on file system
bool path_exists(string const &);
/// Join two paths together
string path_join(string, string const &);

/******************************************************************************************/

/// Exception for a Unix signal
class SignalError : public std::runtime_error {
    int code_;
    static char const * make_name(int i);
public:
    void raise() const;
    string type() const;
    string name() const {return make_name(code_);}

    explicit SignalError(int i) : std::runtime_error(make_name(i)), code_(i) {}
    static SignalError sigint();

    int code() const {return code_;}
};

extern thread_local std::shared_ptr<std::atomic<int>> ThreadLocalSignal;
/// Throw a SignalError if there are any signals
void throw_if_signal();

void set_static_signal(int code);

/******************************************************************************************/

}
