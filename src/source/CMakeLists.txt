cmake_minimum_required(VERSION 3.10)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()

################################################################################

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()
if(POLICY CMP0051)
    cmake_policy(SET CMP0051 NEW)
endif()
if(POLICY CMP0058)
    cmake_policy(SET CMP0058 NEW)
endif()
if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0054 NEW)

################################################################################

set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/external/vcpkg/scripts/buildsystems/vcpkg.cmake
  CACHE STRING "Vcpkg toolchain file")

################################################################################

enable_language(CXX)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(Contents)
include(Options)

################################################################################
project(NUPACK VERSION 4.0.0.23 LANGUAGES CXX C)
################################################################################

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

message("-- CMake build type is \"${CMAKE_BUILD_TYPE}\"")

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

################################################################################

# This function preferentially symlinks a file and records the output in ${link_outputs}
# If on Windows, a copy is done instead of symlinking
function(link_file in out)
    get_filename_component(dir "${CMAKE_BINARY_DIR}/${out}" DIRECTORY)
    if(WIN32)
        set(link_file_flag "copy")
    else()
        set(link_file_flag "create_symlink")
    endif()
    add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/${out}"
        COMMAND ${CMAKE_COMMAND} -E make_directory ${dir}
        COMMAND ${CMAKE_COMMAND} -E ${link_file_flag} "${SOURCE_DIR}/${in}" "${CMAKE_BINARY_DIR}/${out}"
        DEPENDS "${SOURCE_DIR}/${in}"
    )
    set(link_outputs ${link_outputs} "${CMAKE_BINARY_DIR}/${out}" PARENT_SCOPE)
endfunction(link_file)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/rebind EXCLUDE_FROM_ALL)

################################################################################

if(NUPACK_BUILD_CXX)
    include(BuildCXX)
endif()

if(NUPACK_BUILD_PYTHON)
    include(BuildPython)
endif()

if(NUPACK_BUILD_DOCS)
    include(BuildDocs)
endif()

################################################################################
