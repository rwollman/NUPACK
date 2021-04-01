#pragma once
// NUPACK debugging macros

#include "pathway_utils.h"

#include <string>
#include <exception>

#ifdef NUPACK_BACKTRACE
#include <execinfo.h>
#endif

namespace nupack {
namespace custom_csp {
class NupackException : public std::exception {
private:
    int type;
    string message;
    string full_message;
    string location;
public:
    NupackException(int type) {
        this->type = type;
        this->message = "";
        this->location = "";
        set_full_message();
    }

    NupackException(string message) {
        this->type = 0;
        this->message = message;
        this->location = "";
        set_full_message();
    }

    NupackException(string location, string message) {
        this->type = 0;
        this->message = message;
        this->location = location;
        set_full_message();
    }

    void set_full_message() {
        std::stringstream ss;
        if (this->type == 0) {
            if (location.size() > 0) {
                ss << location << ": [ERROR] nupack: " << message;
            } else {
                ss << "[ERROR] nupack: " << message;
            }
        } else {
            ss << "NUPACK error code: " << this->type;
        }
        full_message = ss.str();
    }

    void print_message(std::ostream & out = std::cerr) const {
        out << full_message << std::endl << std::flush;
    }

    const char* what() const noexcept {return full_message.c_str();}
};

void print_backtrace(std::ostream &os, std::size_t n = 10);
string get_backtrace(std::size_t n = 10);

#ifdef NUPACK_BACKTRACE
inline void print_backtrace(std::ostream &os = std::cerr, std::size_t n) {
    vec<void *> symbols(n);
    auto size = backtrace(symbols.data(), n);
    auto strings = backtrace_symbols(symbols.data(), size);
    for (std::size_t i = 0; i != n; ++i) os << strings[i] << std::endl;
}

inline string get_backtrace(std::size_t n) {
    std::stringstream ss;
    ss << "\n**** Backtrace ****\n";
    print_backtrace(ss, n); return ss.str();
}
#else
inline void print_backtrace(std::ostream &os, std::size_t n) {};
inline string get_backtrace(std::size_t n) {return "";}
#endif // NUPACK_BACKTRACE
}
}

#define NUPACK_CHECK(condition, message) if (!(condition)) { NUPACK_DESIGN_ERROR(message)}


#ifdef NDEBUG
#define NUPACK_DESIGN_ERROR(message) throw ::nupack::custom_csp::NupackException( \
    string(message) \
    );

#define NUPACK_DEBUG_PRINT(message)

#define NUPACK_DEBUG_CHECK(condition, message)


#define NUPACK_LOG_ERROR(message) std::cerr << "[ERROR] " << message << std::endl << std::flush;
#define NUPACK_LOG_WARN(message) std::cerr << "[WARN] "  << message << std::endl << std::flush;
#define NUPACK_LOG_INFO(message) std::cerr << "[ctx.info] "  << message << std::endl << std::flush;

#else // NDEBUG

#define NUPACK_DESIGN_ERROR(message) throw ::nupack::custom_csp::NupackException( \
      string(__FILE__) + string(":") + ::nupack::custom_csp::to_string(__LINE__), \
      string(message) + ::nupack::custom_csp::get_backtrace() \
    );

#define NUPACK_LOG_ERROR(message) std::cerr << "[ERROR] " \
  << __FILE__ << ":" << __LINE__ << ":" << message << std::endl << std::flush;

#define NUPACK_LOG_WARN(message) std::cerr << "[WARN] " \
  << __FILE__ << ":" << __LINE__ << ":" << message << std::endl << std::flush;

#define NUPACK_LOG_INFO(message) std::cerr << "[ctx.info] " \
  << __FILE__ << ":" << __LINE__ << ":" << message << std::endl << std::flush;

#define NUPACK_DEBUG_PRINT(message) std::cerr << "[DEBUG] " \
  << __FILE__ << ":" << __LINE__ << ":" <<\
       message  << std::endl << std::flush;

#define NUPACK_DEBUG_CHECK(condition, message) NUPACK_CHECK(condition, message)

#endif // NDEBUG

#define NUPACK_EXC_CHECK(fcall, message) \
  try {fcall; } catch (::nupack::custom_csp::NupackException & e) { \
    e.print_message(std::cerr); NUPACK_DESIGN_ERROR(message); }
