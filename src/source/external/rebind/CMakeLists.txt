cmake_minimum_required(VERSION 3.8)
project(rebind VERSION 0.0.0 LANGUAGES CXX)

################################################################################

set(REBIND_PYTHON "python" CACHE STRING "Specified Python executable used to deduce include directory")
set(REBIND_PYTHON_INCLUDE "" CACHE STRING "Specified include directory containing Python.h")

################################################################################

if (${REBIND_PYTHON_INCLUDE})
    message("-- Using specified Python include")
    set_property(GLOBAL PROPERTY rebind_python_include ${REBIND_PYTHON_INCLUDE})
else()
    execute_process(
        COMMAND ${REBIND_PYTHON} -c "import sys, sysconfig; sys.stdout.write(sysconfig.get_paths()['include'])"
        RESULT_VARIABLE python_stat OUTPUT_VARIABLE python_include
    )
    if (python_stat)
        message(FATAL_ERROR "Failed to deduce include directory from '${REBIND_PYTHON}' executable.\nMaybe specify REBIND_PYTHON_INCLUDE directly.")
    endif()
    message("-- Using Python include directory deduced from REBIND_PYTHON=${REBIND_PYTHON}")
    set_property(GLOBAL PROPERTY rebind_python_include ${python_include})
endif()

message("-- Using Python include directory ${python_include}")

################################################################################

add_library(rebind_interface INTERFACE)
target_compile_features(rebind_interface INTERFACE cxx_std_17)
target_include_directories(rebind_interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

################################################################################

# Maybe change in future to user provided interface library?
add_library(librebind STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Source.cc
)
set_target_properties(librebind PROPERTIES OUTPUT_NAME rebind ${REBIND_STATIC_PROPERTIES})
target_link_libraries(librebind PUBLIC rebind_interface)

################################################################################

# Module.cc has to be recompiled based on the exported module name
# Could just build Python.cc as its own library, but here it's built together with Module.cc
set_property(GLOBAL PROPERTY rebind_module_files
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Python.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Module.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Cast.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/source/Globals.cc
)

function(rebind_module module_name output_name)
    get_property(files GLOBAL PROPERTY rebind_module_files)
    add_library(${module_name} MODULE ${files})
    set_target_properties(${module_name} PROPERTIES PREFIX "" OUTPUT_NAME ${output_name})

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set_target_properties(${module_name} PROPERTIES LINK_FLAGS "-Wl,-flat_namespace,-undefined,dynamic_lookup")
    else()
        set_target_properties(${module_name} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()

    target_compile_definitions(${module_name} PRIVATE REBIND_MODULE=${output_name})
    target_link_libraries(${module_name} PRIVATE rebind_interface ${ARGN})

    get_property(python_include GLOBAL PROPERTY rebind_python_include)
    target_include_directories(${module_name} PRIVATE ${python_include})
endfunction(rebind_module)

################################################################################

add_library(librebindtest SHARED EXCLUDE_FROM_ALL source/Test.cc)
set_target_properties(librebindtest PROPERTIES OUTPUT_NAME "rebindtest")
target_link_libraries(librebindtest PRIVATE librebind)
rebind_module(rebindtest rebindtest librebindtest)

################################################################################

set(REBIND_PYTHON_FILES
    rebind/__init__.py
    rebind/blank.py
    rebind/common.py
    rebind/dispatch.py
    rebind/render.py
    rebind/types.py
    CACHE INTERNAL "List of Python files in the rebind module"
)
