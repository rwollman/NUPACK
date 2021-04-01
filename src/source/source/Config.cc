/**
 * @brief Some global constant definitions, including version information
 *
 * @file Config.cc
 * @author Mark Fornace
 * @date 2018-05-31
 */
#include <nupack/Version.h>
#include <nupack/common/Config.h>
#include <nupack/common/Error.h>

#include <vector>
#include <regex>
#include <array>
#include <thread>

/******************************************************************************************/

namespace nupack {

void time_sink_impl(void const *) {}

/******************************************************************************************/

#define NUPACK_TMP0(x) #x
#define NUPACK_TMP(x) NUPACK_TMP0(x)

string DefaultParametersPath = NUPACK_TMP(NUPACK_PARAMETERS),
       DefaultDataPath =       NUPACK_TMP(NUPACK_DATA),
       MatlabCommand =         NUPACK_TMP(NUPACK_MATLAB),
       MathematicaCommand =    NUPACK_TMP(NUPACK_MATHEMATICA_PATH);

string const Version =     NUPACK_TMP(NUPACK_VERSION_MAJOR) "." NUPACK_TMP(NUPACK_VERSION_MINOR) "." NUPACK_TMP(NUPACK_VERSION_PATCH),
             GitRevision = NUPACK_TMP(NUPACK_GIT_REVISION),
             GitBranch =   NUPACK_TMP(NUPACK_GIT_BRANCH);

#undef NUPACK_TMP
#undef NUPACK_TMP0

bool DebugInfo = NUPACK_DEBUG_INFO;
unsigned int TotalCPU = std::max<unsigned int>(1u, std::thread::hardware_concurrency());
std::size_t TotalRAM = NUPACK_RAM_IN_MB * 1e6;

/******************************************************************************************/

string message_string(string fn, int line, string msg) {
    std::stringstream buf; buf << std::boolalpha;
    buf << ": \"" << msg << "\" (" << fn << ", line " << line << ")\n";
    return buf.str();
}

/******************************************************************************************/

struct Regex {
    string replace;
    std::regex re;
    Regex(char const *pattern, char const *replace_) : replace(replace_),
        re(pattern, std::regex::ECMAScript | std::regex::optimize) {}
    void operator()(string &s) const {s = std::regex_replace(std::move(s), re, replace);}
};

string trim_type_name(string s, int n) {
    static std::array<Regex, 31> const regexes = {{
        {"(std|__1|boost|container|nupack|detail|python|thermo|arma)::", ""},
        {"( >)",                                            ">"},
        {"(, *char_traits<[^<>]*>)",                        ""},
        {"(, *allocator<[^<>]*>)",                          ""},
        {"(, *allocator<[^<>]*<[^<>]*>[^<>]*>)",            ""},
        {"(, *new_allocator<[^<>]*>)",                      ""},
        {"basic_string<char>",                              "string"},
        {"double",                                          "d"},
        {"float",                                           "f"},
        {"int",                                             "i"},
        {"unsigned int",                                    "u"},
        {"unsigned long",                                   "ul"},
        {"__list_iterator",                                 "list_iter"},
        {"__wrap_iter<([^<>]*)>",                           "$1"},
        {"integral_constant<bool, true>",                   "True"},
        {"integral_constant<bool, false>",                  "False"},
        {"\\b([^,]{3,})((, \\1){15,})",                     "16+ * $1"},
        {"\\b([^,]{3,})((, \\1){14})",                      "15 * $1"},
        {"\\b([^,]{3,})((, \\1){13})",                      "14 * $1"},
        {"\\b([^,]{3,})((, \\1){12})",                      "13 * $1"},
        {"\\b([^,]{3,})((, \\1){11})",                      "12 * $1"},
        {"\\b([^,]{3,})((, \\1){10})",                      "11 * $1"},
        {"\\b([^,]{3,})((, \\1){9})",                       "10 * $1"},
        {"\\b([^,]{3,})((, \\1){8})",                       "9 * $1"},
        {"\\b([^,]{3,})((, \\1){7})",                       "8 * $1"},
        {"\\b([^,]{3,})((, \\1){6})",                       "7 * $1"},
        {"\\b([^,]{3,})((, \\1){5})",                       "6 * $1"},
        {"\\b([^,]{3,})((, \\1){4})",                       "5 * $1"},
        {"\\b([^,]{3,})((, \\1){3})",                       "4 * $1"},
        {"\\b([^,]{3,})((, \\1){2})",                       "3 * $1"},
        {"\\b([^,]{3,})((, \\1){1})",                       "2 * $1"},
    }};
    for (auto const &r : regexes) r(s);
    std::replace(s.begin(), s.end(), '<', '(');
    std::replace(s.begin(), s.end(), '>', ')');
    if (s.size() < n+2) return s;
    auto it = std::find_if(s.begin() + n, s.end(), [](auto c) {return !std::isalnum(c) && c != '_';});
    s.erase(it, s.end());
    return s + "...";
}

/******************************************************************************************/

}
