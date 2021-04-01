/**
 * @brief Defines type aliases and a few global constants
 *
 * @file Config.h
 * @author Mark Fornace
 * @date 2018-05-31
 */
#pragma once

#include <cstddef>
#include <cstring>
#include <string>
#include <string_view>
#include <iomanip>
#include <complex>
#include <limits>

/******************************************************************************************/

/// These overloads currently interfere with some printing functionality
#define BOOST_NO_IOSTREAM

/******************************************************************************************/

/// File name macro which excludes directory path
#define NUPACK_FILE (std::strrchr(__FILE__, '/') ? std::strrchr(__FILE__, '/') + 1 : \
    std::strrchr(__FILE__, '\\') ? std::strrchr(__FILE__, '\\') + 1 : __FILE__)

/******************************************************************************************/

namespace nupack {

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    static constexpr bool is_windows = true;
#else
    static constexpr bool is_windows = false;
#endif

/******************************************************************************************/

/// Some simple default types for use across the whole project
using real         = double;
using real32       = float;
using real64       = double;
using string       = std::string;
using string_view  = std::string_view;

using namespace std::string_literals;

using uint         = unsigned int;
using ushort       = unsigned short;
using complex_real = std::complex<real>;

using iseq         = unsigned int;
static_assert(std::numeric_limits<iseq>::max() >= 1e6, "iseq constrains maximum sequence length");

using usize        = std::size_t;
static_assert(sizeof(usize) >= sizeof(std::size_t), "usize should be at least size_t size");

/******************************************************************************************/

/// Some string values that will be set by CMake
extern string const GitRevision, GitBranch, Version;
extern string DefaultParametersPath, DefaultDataPath, MatlabCommand, MathematicaCommand;
/// Number of logical CPU cores - you can change this if desired
extern unsigned int TotalCPU;
/// Total RAM in bytes - you can change this if desired
extern std::size_t TotalRAM;
/// Print backtraces in exceptions
extern bool DebugInfo;

/******************************************************************************************/

/// constexpr Debug, should be optimized out by compiler when put in if()
#if NUPACK_DEBUG
#   if NUPACK_DEBUG == 1
        static constexpr bool const Release = false;
        static constexpr bool const Debug = false;
        static constexpr bool const DebugBounds = false;
#   endif
#   if NUPACK_DEBUG == 2
        static constexpr bool const Release = false;
        static constexpr bool const Debug = true;
        static constexpr bool const DebugBounds = true;
#   endif
#else
    static constexpr bool const Release = true;
    static constexpr bool const Debug = false;
    static constexpr bool const DebugBounds = false;
#endif

/******************************************************************************************/

}
