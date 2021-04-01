option(NUPACK_BUILD_SUPER        "Whether or not a superbuild should be invoked" OFF)
option(NUPACK_BUILD_CXX          "Whether or not to build C++ libraries"         ON)
option(NUPACK_BUILD_TESTS        "Whether or not to build C++ unit tests"        ON)
option(NUPACK_BUILD_PYTHON       "Whether or not to build Python bindings"       ON)
option(NUPACK_BUILD_DOCS         "Whether or not to build website documentation" OFF)
option(NUPACK_SHARED             "Build libnupack as a shared library"           OFF)
option(NUPACK_CCACHE             "Enable CCache"                                 ON)
option(NUPACK_SERIALIZE          "Enable serialization via JSON and msgpack"     ON)
option(NUPACK_PIC                "Use position-independent code linkage"         ${CMAKE_POSITION_INDEPENDENT_CODE})
option(NUPACK_DETERMINISTIC      "Use deterministic random number generator"     OFF)
option(NUPACK_FORTRAN            "Compile Fortran sandbox"                       OFF)
option(NUPACK_IWYU               "Enable include-what-you-use"                   OFF)
option(NUPACK_HDF5               "Enable HDF5"                                   OFF)
option(NUPACK_MATHEMATICA        "Enable Mathematica"                            OFF)
option(NUPACK_MATLAB             "Enable MATLAB"                                 OFF)
option(NUPACK_MLPACK             "Enable MLPACK"                                 OFF)
option(NUPACK_MPI                "Enable MPI compilation"                        OFF)
option(NUPACK_OMP                "Enable OpenMP"                                 OFF)
option(NUPACK_ONLY               "Compile only a given file"                     OFF)
option(NUPACK_DESIGN_ONLY        "Compile only given files with design library"  OFF)
option(NUPACK_PGO                "Enable PGO, can be OFF, READ, or WRITE"        OFF)
option(NUPACK_EXTERNAL_ARMADILLO "Use external version of armadillo"             OFF)

################################################################################

message(STATUS "--------------------------------------------------------------------------------")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" nupack_arch_native)

if(NUPACK_SIMD_FLAGS)
    message(STATUS "-- Using SIMD architecture flags \"${NUPACK_SIMD_FLAGS}\"")
elseif(nupack_arch_native)
    message(STATUS "-- Using \"-march=native\" for SIMD architecture flags")
    set(NUPACK_SIMD_FLAGS "-march=native" CACHE STRING "SIMD architecture flags to use when compiling")
else()
    message(STATUS "-- Not using any SIMD architecture flags")
endif()

################################################################################

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
    "#include <shared_mutex>
    std::shared_timed_mutex mut;
    int main(){
        auto x = std::unique_lock<std::shared_timed_mutex>(mut);
        auto y = std::shared_lock<std::shared_timed_mutex>(mut);
        return 0;
    }"
    can_use_shared_mutex
)

################################################################################

cmake_host_system_information(RESULT ram_in_mb QUERY TOTAL_PHYSICAL_MEMORY)
set(NUPACK_RAM_IN_MB  ${ram_in_mb} CACHE STRING "Maximum amount of RAM")
message(STATUS "-- Assuming maximum RAM to use is ${NUPACK_RAM_IN_MB} MB")

################################################################################

# for each path in a list, add a directory prefix to it
function(prefix_transform var prefix)
    set(listVar "")
    foreach(f ${ARGN})
        list(APPEND listVar "${prefix}/${f}")
    endforeach(f)
    set(${var} "${listVar}" PARENT_SCOPE)
endfunction(prefix_transform)
