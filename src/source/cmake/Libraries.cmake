################################################################################

find_package(TBB CONFIG REQUIRED)

################################################################################

add_library(nupack_lapack INTERFACE)

if(NOT LAPACK_LIBRARIES)
    find_package(LAPACK CONFIG)
endif()

if(NOT LAPACK_LIBRARIES)
    find_library(LAPACK_LIBRARIES lapack)
endif()

if(LAPACK_LIBRARIES)
    message("-- Using LAPACK from ${LAPACK_LIBRARIES}")
    target_link_libraries(nupack_lapack INTERFACE ${LAPACK_LIBRARIES})
else()
    message(FATAL_ERROR "Could not find LAPACK")
endif()

################################################################################

find_package(nlohmann_json CONFIG REQUIRED)
add_library(json INTERFACE)
target_link_libraries(json INTERFACE nlohmann_json nlohmann_json::nlohmann_json)

################################################################################

add_library(armadillo INTERFACE)
find_package(Armadillo CONFIG REQUIRED)
target_link_libraries(armadillo INTERFACE ${ARMADILLO_LIBRARIES})
target_compile_definitions(armadillo INTERFACE ARMA_DONT_USE_WRAPPER=1 ARMA_USE_LAPACK=1)

################################################################################

if(NUPACK_CCACHE)
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
        message("-- Using ccache from ${CCACHE_FOUND}")
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
    else()
        message("-- Not using ccache since it could not be found")
    endif(CCACHE_FOUND)
else()
    message("-- Not using ccache")
endif()

################################################################################

if(${NUPACK_IWYU} STREQUAL "ON")
    find_program(IWYU_PATH NAMES include-what-you-use iwyu PATHS /usr/local/opt/llvm/bin /usr/local/bin)
    if(IWYU_PATH)
        message("-- Using include-what-you-use from ${IWYU_PATH}")
    else()
        message(FATAL_ERROR "Could not find the program include-what-you-use")
    endif()
else()
    message("-- Not using include-what-you-use")
endif()

################################################################################

if(${NUPACK_PGO} STREQUAL "WRITE")
    message("-- Writing profile information")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate")
endif()

if(${NUPACK_PGO} STREQUAL "READ")
    message("-- Using profile guided optimization")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-use=NUPACK_${NUPACK_GIT_REVISION}.profdata")
endif()

################################################################################

add_library(backward INTERFACE)
target_include_directories(backward INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/backward-cpp)

################################################################################

find_package(fmt CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

################################################################################

add_library(nupack-boost INTERFACE)
find_path(BOOST_CORE_INCLUDE_DIRS "boost/checked_delete.hpp" REQUIRED)
find_path(BOOST_PREPROCESSOR_INCLUDE_DIRS "boost/preprocessor.hpp" REQUIRED)
find_path(BOOST_FUNCTIONAL_INCLUDE_DIRS "boost/functional.hpp" REQUIRED)
find_path(BOOST_CONTAINER_INCLUDE_DIRS "boost/container/small_vector.hpp" REQUIRED)
find_path(BOOST_VARIANT_INCLUDE_DIRS "boost/variant.hpp" REQUIRED)
find_path(BOOST_ITERATOR_INCLUDE_DIRS "boost/function_output_iterator.hpp" REQUIRED)
find_path(BOOST_ALIGN_INCLUDE_DIRS "boost/align.hpp" REQUIRED)
find_path(BOOST_SORT_INCLUDE_DIRS "boost/sort/spreadsort/spreadsort.hpp" REQUIRED)

target_include_directories(nupack-boost INTERFACE
    ${BOOST_CORE_INCLUDE_DIRS}
    ${BOOST_PREPROCESSOR_INCLUDE_DIRS}
    ${BOOST_FUNCTIONAL_INCLUDE_DIRS}
    ${BOOST_CONTAINER_INCLUDE_DIRS}
    ${BOOST_FUSION_INCLUDE_DIRS}
    ${BOOST_VARIANT_INCLUDE_DIRS}
    ${BOOST_ITERATOR_INCLUDE_DIRS}
    ${BOOST_ALIGN_INCLUDE_DIRS}
    ${BOOST_SORT_INCLUDE_DIRS}
)

################################################################################

add_library(boostsimd INTERFACE)
target_include_directories(boostsimd INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/boost-simd/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

################################################################################

add_library(nupack-gecode INTERFACE)
find_path(GECODE_INCLUDE_DIR gecode/search.hh REQUIRED)
target_include_directories(nupack-gecode INTERFACE ${GECODE_INCLUDE_DIR})
find_library(GECODE_MINIMODEL gecodeminimodel REQUIRED)
find_library(GECODE_DRIVER gecodedriver REQUIRED)
find_library(GECODE_KERNEL gecodekernel REQUIRED)
find_library(GECODE_INT gecodeint REQUIRED)
find_library(GECODE_SUPPORT gecodesupport REQUIRED)
find_library(GECODE_SET gecodeset REQUIRED)
find_library(GECODE_FLOAT gecodefloat REQUIRED)
find_library(GECODE_SEARCH gecodesearch REQUIRED)
target_link_libraries(nupack-gecode INTERFACE ${GECODE_MINIMODEL} ${GECODE_DRIVER}
    ${GECODE_INT} ${GECODE_SET} ${GECODE_FLOAT} ${GECODE_SEARCH} ${GECODE_KERNEL} ${GECODE_SUPPORT})

################################################################################

message("--------------------------------------------------------------------------------")
