cmake_minimum_required(VERSION 2.6)
project(NUPACKCORE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(NUPACK_VERSION 3.2.2)
add_definitions(-DCMAKE_NUPACK_VERSION="${NUPACK_VERSION}")


# Packaging Binary and Source
set(CPACK_PACKAGE_VERSION_MAJOR 3)
set(CPACK_PACKAGE_VERSION_MINOR 2)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_DIRECTORY "${PROJECT_BINARY_DIR}/cpack")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "nupack${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES "/build.*/" "/debug/" "/[.].*/" "/[.].*" ${CPACK_SOURCE_IGNORE_FILES})

include(CPack)


include(CheckCXXCompilerFlag)
include(CheckIncludeFiles)

# Set a default build type if none was specified

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build."
            FORCE)
endif()

########################################
# Selection options
option(SAMPLE "SAMPLE" ON)

########################################
# Location options
set(TEMP_LIB_DIR ${PROJECT_BINARY_DIR}/lib)
set(BIN_DIR bin)

########################################
# Include directories
set(INCLUDE_DIRS
    extlib
    src
)

########################################
# Build options
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${TEMP_LIB_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEMP_LIB_DIR})

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O0 -fno-omit-frame-pointer -fsanitize=address")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g")

if(CMAKE_C_COMPILER MATCHES "clang")
message("-- using Clang Compiler")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g")
endif(CMAKE_C_COMPILER MATCHES "clang")

if(CMAKE_CXX_COMPILER MATCHES "clang")
message("-- using Clang++ Compiler")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif(CMAKE_CXX_COMPILER MATCHES "clang")

check_cxx_compiler_flag("-std=c++11" HAS_CXX11)
check_cxx_compiler_flag("-std=c++0x" HAS_CXX0X)

if (HAS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (HAS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else (HAS_CXX11)
    message("Compiler doesn't support -std=c++0x or -std=c++11")
    message("Some features may not build correctly")
endif (HAS_CXX11)

# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_pie -lc")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -O0 -fno-omit-frame-pointer")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fno-omit-frame-pointer")

if(SAMPLE)
    add_definitions(-DNUPACK_SAMPLE)
endif(SAMPLE)

set(CLANG_ASAN "")

# Some strangeness that can be refactored out; minimal changes to allow 
# access to included jsoncpp instead of system.
set(JSONCPP_FOUND 1)
set(JSONCPP_LIBRARY jsoncpp)
set(JSONCPP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/extlib/jsoncpp)
message("-- Found jsoncpp library: " ${JSONCPP_LIBRARY})
message("-- Found jsoncpp include directory: " ${JSONCPP_INCLUDE_DIR})

include_directories(BEFORE ${INCLUDE_DIRS})
link_directories(${TEMP_LIB_DIR})

find_library(MATH_LIB NAMES libm m)
if(MATH_LIB)
    message("-- MATH_LIB is ${MATH_LIB}")
else()
    message("-- MATH_LIB not found. Using -lm linker flag.")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm")
    set(MATH_LIB "")
endif()

# Configure NUPACK_VERSION using Git if it's available

find_package(Git)
if(GIT_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    message("-- Found Git")
    
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} 
            TIMEOUT 3
            OUTPUT_VARIABLE GIT_REVISION
            RESULT_VARIABLE FOUND_REVISION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            )
    
    if(FOUND_REVISION EQUAL 0)
        set(NUPACK_VERSION "${GIT_REVISION}")
    endif(FOUND_REVISION EQUAL 0)
    
endif(GIT_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
message(STATUS "Building version ${NUPACK_VERSION}")

########################################
# Install options
set(PARAMETER_INSTALL_LOCATION share/nupack)
set(INCLUDE_INSTALL_LOCATION include)
set(LIBRARY_INSTALL_LOCATION lib)

install(DIRECTORY "parameters/" DESTINATION ${PARAMETER_INSTALL_LOCATION})

find_package(Doxygen)
if(DOXYGEN_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
    add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR/doc}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)

function(INSTALL_INCLUDE_TREE dir)
    file(GLOB_RECURSE deploy_files "${dir}/*.h" "${dir}/*.hpp")
    if(${ARGN})
        list(REMOVE_ITEM deploy_files ${ARGN})
    endif()
    foreach(file_path ${deploy_files})
        file(RELATIVE_PATH rel_path ${CMAKE_SOURCE_DIR}/src ${file_path})
        get_filename_component(dir_path ${rel_path} PATH)
        install(FILES ${file_path} DESTINATION ${INCLUDE_INSTALL_LOCATION}/nupack/${dir_path})
    endforeach()
endfunction()

# Backtrace stuff
########################################
option(BACKTRACE "Print out backtrace in exceptions" OFF)
if(${BACKTRACE} STREQUAL "ON")
    CHECK_CXX_SOURCE_COMPILES(
        "#include <execinfo.h>
        void *array[10];
        auto size = backtrace(array, 10);
        auto strings = backtrace_symbols(array, size);
        int main() {return 0;}"
        CAN_BACKTRACE
    )
    if(${CAN_BACKTRACE})
        add_definitions(-DNUPACK_BACKTRACE)
        message("-- Printing backtrace in exceptions")
    else()
        message("-- Cannot print backtrace in exceptions")
    endif()
else()
    message("-- Not printing backtrace in exceptions")
endif()
########################################

function(ADD_DIR_IF_EXISTS dir)
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
        add_subdirectory(${dir})
    endif()
endfunction(ADD_DIR_IF_EXISTS)

add_subdirectory(extlib)
add_subdirectory(src)
