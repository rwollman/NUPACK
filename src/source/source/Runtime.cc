/**
 * @brief Definitions for environmental variable lookup, type name strings, demangling, signal handlers
 *
 * @file Runtime.cc
 * @author Mark Fornace
 * @date 2018-05-31
 */
#include <nupack/common/Runtime.h>
#include <nupack/common/Error.h>

#include <csignal>
#include <atomic>
#include <backward.hpp>
#include <boost/core/demangle.hpp>

#include <iostream>

namespace nupack {

thread_local std::shared_ptr<std::atomic<int>> ThreadLocalSignal = {};


/******************************************************************************************/

/// get environmental variable or else ""
string get_env(string key) {
    char const * const tmp = std::getenv(key.c_str());
    if (tmp == nullptr) return {};
    else return tmp;
}

bool path_exists(string const &s) {return std::ifstream(s).good();}

string path_join(string s, string const &t) {
    s.push_back(is_windows ? '\\' : '/');
    s += t;
    return s;
}

void spin_processor(chrono::duration<real> const t) {
    auto const t0 = chrono::high_resolution_clock::now();
    real a = 0.5, output = a;
    std::size_t const limit = std::size_t(t.count()) * (std::size_t(2) << 26u);
    for (std::size_t i = 0; i != limit; ++i)
        output = 0.5 * (a / output + output); // Babylonian method for square root
    time_sink(output);
}

/******************************************************************************************/

string short_type_name(string s) {
    // strip off all after colon or non-identifier character
    std::replace(s.begin(), s.end(), ' ', '_');
    std::replace(s.begin(), s.end(), '*', 'p');
    auto it = std::find_if_not(s.begin(), s.end(), [](auto c) {return c == '_' || c == ':' || std::isalnum(c);});
    if (it != s.end()) s.erase(it, s.end());
    // strip off all before last colon
    auto i = s.find_last_of(':');
    if (i != string::npos) s.erase(s.begin(), s.begin()+i+1);
    s.erase(s.begin() + s.find_last_not_of('_') + 1, s.end());
    s.erase(s.begin(), s.begin() + s.find_first_not_of('_'));
    s.shrink_to_fit();
    return s;
}

/******************************************************************************************/

string demangle(string s) {return boost::core::demangle(s.c_str());}

string backtrace_demangle(string s) {
    auto b = s.find('_'); // start of symbol
    auto e = s.find(' ', b); // space after symbol
    if (b == string::npos || b == 0 || e == string::npos) return s;
    s[e] = '\0'; // cut off non-symbol part
    return boost::core::demangle(s.c_str() + b);
}

void print_backtrace(std::ostream &os, std::size_t n) {
    backward::StackTrace st;
    st.load_here(n);
    backward::Printer p;
    p.color_mode = backward::ColorMode::always;
    p.print(st, os);
}

string get_backtrace(std::size_t n) {
    std::ostringstream ss;
    print_backtrace(ss, n);
    return ss.str();
}

Error::Error(std::string s) : std::runtime_error(std::move(s)) {
    if (false && DebugInfo) print_backtrace(std::cerr);
}

/******************************************************************************************/

void SignalError::raise() const {
    if (std::raise(code()))
        throw std::runtime_error("Could not raise signal for code " + std::to_string(code()));
}

SignalError SignalError::sigint() {return SignalError(SIGINT);}

string SignalError::type() const {
    switch(code()) {
        case(SIGINT):  return "SIGINT";
        case(SIGILL):  return "SIGILL";
        case(SIGFPE):  return "SIGFPE";
        case(SIGSEGV): return "SIGSEGV";
        case(SIGTERM): return "SIGTERM";
        case(SIGABRT): return "SIGABRT";
        default: return "Unknown signal";
    }
}

char const * SignalError::make_name(int i) {
    switch(i) {
        case(SIGINT):  return "SIGINT  - Terminal interrupt signal";
        case(SIGILL):  return "SIGILL  - Illegal instruction signal";
        case(SIGFPE):  return "SIGFPE  - Floating point error signal";
        case(SIGSEGV): return "SIGSEGV - Segmentation violation signal";
        case(SIGTERM): return "SIGTERM - Termination request signal";
        case(SIGABRT): return "SIGABRT - Abort (abnormal termination) signal";
        default: return "Unknown signal";
    }
};

/******************************************************************************************/

struct RuntimeContents {
#   ifndef _WIN32
        struct sigaction handler;
#   endif
    backward::SignalHandling sh;
    std::vector<int> signals;
};
std::sig_atomic_t StaticSignal = 0;

void throw_if_signal() {
    int i = 0;
    // Check signal set on this thread
    if (ThreadLocalSignal) i = ThreadLocalSignal->exchange(0, std::memory_order_relaxed);
    // Check system level signal
    if (!i) {i = StaticSignal; StaticSignal = 0;}
    if (i) throw SignalError(i);
}

extern "C" void set_nupack_static_signal(int i) {StaticSignal = i;}

void set_static_signal(int code) {set_nupack_static_signal(code);}

SignalRuntime::SignalRuntime(int const *b, int const *e) : contents(std::make_shared<RuntimeContents>()) {
    contents->signals.assign(b, e);

    if (!contents->sh.loaded()) std::cout << "C++ backtrace handler not loaded" << std::endl;
    std::lock_guard<std::mutex> lock(mut);
    for (auto i : contents->signals) {
#   ifdef _WIN32
        signal(i, set_nupack_static_signal);
#   else
        contents->handler.sa_handler = set_nupack_static_signal;
        sigemptyset(&contents->handler.sa_mask);
        contents->handler.sa_flags = 0;
        if (sigaction(i, &contents->handler, &contents->handler))
            throw std::runtime_error("Error setting sigaction handler");
#   endif
    }
}

constexpr int sigint_value = SIGINT;

SignalRuntime::SignalRuntime() : SignalRuntime(&sigint_value, &sigint_value+1) {}

SignalRuntime::~SignalRuntime() {
    std::lock_guard<std::mutex> lock(mut);
    for (auto i : contents->signals) {
#   ifdef _WIN32
        signal(i, SIG_DFL);
#   else
        if (!sigaction(i, &contents->handler, nullptr)) return;
        std::cerr << "Error resetting sigaction handler" << std::endl;
        std::terminate();
#   endif
    }
}

/******************************************************************************************/

}
