################################################################################

set(link_outputs)

foreach(x ${NUPACK_PYTHON_FILES})
    link_file(python/${x} nupack/${x})
endforeach()

foreach(x dna04.json rna95.json rna06.json rna99.json dna04-nupack3.json rna95-nupack3.json rna99-nupack3.json)
    link_file(parameters/${x} nupack/parameters/${x})
endforeach()

foreach(x ${NUPACK_PYTHON_TEST_FILES})
    link_file(test/python/${x} nupack/test/${x})
endforeach()

foreach(x ${REBIND_PYTHON_FILES})
    link_file(external/rebind/${REBIND_PYTHON_ROOT}/${x} nupack/${x})
endforeach()

configure_file(package/setup.py ${CMAKE_BINARY_DIR}/setup.py)
configure_file(LICENSE.txt ${BUILD_DIR}/LICENSE.txt)
# file(WRITE ${BUILD_DIR}/build.sh "echo 'Installing nupack module'\n$PYTHON -m pip install --no-deps --ignore-installed .\n\n")

configure_file(package/info.py ${CMAKE_BINARY_DIR}/nupack/info.py)

add_custom_target(copy_python ALL DEPENDS ${link_outputs})

################################################################################

if(NUPACK_BUILD_CXX)
    rebind_module(nupack-python cpp nupack-bind)
else()
    set(NUPACK_BIND_LIBRARY)
    if(NUPACK_BIND_LIBRARY)
        set(NUPACK_BIND_PATH ${NUPACK_BIND_LIBRARY})
        message("-- Using binding library ${NUPACK_BIND_PATH}")
    else()
        find_library(NUPACK_BIND_PATH NAMES nupack-bind)
        if(NUPACK_BIND_PATH)
            message("-- Found binding library ${NUPACK_BIND_PATH}")
        else()
            message(FATAL_ERROR "libnupack-bind cannot be found")
        endif()
    endif()
    add_library(nupack-bind-import INTERFACE)
    target_link_libraries(nupack-bind-import INTERFACE ${NUPACK_BIND_PATH})
    rebind_module(nupack-python cpp nupack-bind-import)
endif()

set_target_properties(nupack-python PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/nupack) # CXX_VISIBILITY_PRESET hidden)
list(APPEND nupack_cxx_targets nupack-python)
add_dependencies(nupack-python copy_python)
add_library(nupack::python ALIAS nupack-python)
