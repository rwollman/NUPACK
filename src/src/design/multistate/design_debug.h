#pragma once
// NUPACK debugging macros

#include "pathway_utils.h"

#include <string>
#include <exception>

#ifdef NUPACK_BACKTRACE
#include <execinfo.h>
#endif

namespace nupack {
class NupackException : public std::exception {
  private:
    int type;
    std::string message;
    std::string location;
  public:
    NupackException(int type);
    NupackException(std::string message);
    NupackException(std::string location, std::string message);
    void print_message(std::ostream & = std::cerr);
};

inline NupackException::NupackException(int type) {
  this->type = type;
  this->message = "";
  this->location = "";
}

inline NupackException::NupackException(std::string message) {
  this->type = 0;
  this->message = message;
  this->location = "";
}

inline NupackException::NupackException(std::string location, std::string message) {
  this->type = 0;
  this->message = message;
  this->location = location;
}

inline void NupackException::print_message(std::ostream & out) {
  if (this->type == 0) {
    if (location.size() > 0) {
      out << location << ": [ERROR] nupack: " << message << std::endl << std::flush;
    } else {
      out << "[ERROR] nupack: " << message << std::endl << std::flush; 
    }
  } else {
    out << "NUPACK error code: " << this->type << std::endl << std::flush;
  }
}

void print_backtrace(std::ostream &os, std::size_t n=10);
std::string get_backtrace(std::size_t n=10);

#ifdef NUPACK_BACKTRACE
    inline void print_backtrace(std::ostream &os=std::cerr, std::size_t n) {
        std::vector<void *> symbols(n);
        auto size = backtrace(symbols.data(), n);
        auto strings = backtrace_symbols(symbols.data(), size);
        for (std::size_t i = 0; i != n; ++i) os << strings[i] << std::endl;
    }

    inline std::string get_backtrace(std::size_t n) {
        std::stringstream ss;
        ss << "\n**** Backtrace ****\n";
        print_backtrace(ss, n); return ss.str();
    }
#else
    inline void print_backtrace(std::ostream &os, std::size_t n) {};
    inline std::string get_backtrace(std::size_t n) {return "";}
#endif // NUPACK_BACKTRACE
}

#define NUPACK_CHECK(condition, message) if (!(condition)) { NUPACK_ERROR(message)}


#ifdef NDEBUG
#define NUPACK_ERROR(message) throw ::nupack::NupackException( \
    std::string(message) \
    );

#define NUPACK_DEBUG(message) 

#define NUPACK_DEBUG_CHECK(condition, message) 


#define NUPACK_LOG_ERROR(message) std::cerr << "[ERROR] " << message << std::endl << std::flush;
#define NUPACK_LOG_WARN(message) std::cerr << "[WARN] "  << message << std::endl << std::flush;
#define NUPACK_LOG_INFO(message) std::cerr << "[INFO] "  << message << std::endl << std::flush;

#else // NDEBUG

#define NUPACK_ERROR(message) throw ::nupack::NupackException( \
      std::string(__FILE__) + std::string(":") + ::nupack::to_string(__LINE__), \
      std::string(message) + ::nupack::get_backtrace() \
    );

#define NUPACK_LOG_ERROR(message) std::cerr << "[ERROR] " \
  << __FILE__ << ":" << __LINE__ << ":" << message << std::endl << std::flush;

#define NUPACK_LOG_WARN(message) std::cerr << "[WARN] " \
  << __FILE__ << ":" << __LINE__ << ":" << message << std::endl << std::flush;

#define NUPACK_LOG_INFO(message) std::cerr << "[INFO] " \
  << __FILE__ << ":" << __LINE__ << ":" << message << std::endl << std::flush;

#define NUPACK_DEBUG(message) std::cerr << "[DEBUG] " \
  << __FILE__ << ":" << __LINE__ << ":" <<\
       message  << std::endl << std::flush;

#define NUPACK_DEBUG_CHECK(condition, message) NUPACK_CHECK(condition, message)

#endif // NDEBUG

#define NUPACK_EXC_CHECK(fcall, message) \
  try {fcall; } catch (::nupack::NupackException & e) { \
    e.print_message(std::cerr); NUPACK_ERROR(message); }
