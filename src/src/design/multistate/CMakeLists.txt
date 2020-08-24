cmake_minimum_required(VERSION 2.6)

include(CheckIncludeFileCXX)
check_include_file_cxx(cstdint HAVE_CSTDINT_H)
if(${HAVE_CSTDINT_H})
    add_definitions(-DHAVE_CSTDINT_H)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
if(${JSONCPP_FOUND})
    add_definitions(-DJSONCPP_FOUND)
    set(JSONCPP_LIB ${JSONCPP_LIBRARY})
    include_directories(${JSONCPP_INCLUDE_DIR})
endif()


if(BISONFLEX)
    BISON_TARGET(PATHWAYPARSER pathway_parser.y 
                 ${CMAKE_CURRENT_SOURCE_DIR}/pathway_parser.c
                 COMPILE_FLAGS 
                 "--warnings=all --verbose --report-file=${CMAKE_CURRENT_BINARY_DIR}/verbose_parser_output.txt"
    )
    
    FLEX_TARGET(PATHWAYSCANNER pathway_parser.l 
                ${CMAKE_CURRENT_SOURCE_DIR}/pathway_lexer.c 
                COMPILE_FLAGS --header-file=${CMAKE_CURRENT_SOURCE_DIR}/pathway_lexer.h
    )
    ADD_FLEX_BISON_DEPENDENCY(PATHWAYSCANNER PATHWAYPARSER)
    
    string(REGEX REPLACE "/" "_" PATH_STRING ${CMAKE_CURRENT_SOURCE_DIR})
    string(TOUPPER ${PATH_STRING} PATH_STRING)
    
    add_custom_command(OUTPUT mock.out1
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${FLEX_PATHWAYSCANNER_OUTPUTS}
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${CMAKE_CURRENT_SOURCE_DIR}/pathway_lexer.h
        DEPENDS ${FLEX_PATHWAYSCANNER_OUTPUTS} ${CMAKE_CURRENT_SOURCE_DIR}/pathway_lexer.h
        )
    add_custom_command(OUTPUT mock.out2
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${BISON_PATHWAYPARSER_OUTPUT_SOURCE}
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${BISON_PATHWAYPARSER_OUTPUT_HEADER}
        COMMAND sed -e \"s:${PATH_STRING}::g\" -i \"\" ${BISON_PATHWAYPARSER_OUTPUT_SOURCE}
        COMMAND sed -e \"s:${PATH_STRING}::g\" -i \"\" ${BISON_PATHWAYPARSER_OUTPUT_HEADER}
        DEPENDS ${BISON_PATHWAYPARSER_OUTPUT_SOURCE} ${BISON_PATHWAYPARSER_OUTPUT_HEADER}
        )
    add_custom_target(msdesign_rename ALL DEPENDS mock.out1 mock.out2)
else(BISONFLEX)
    set(BISON_PATHWAYPARSER_OUTPUTS pathway_parser.c)
    set(FLEX_PATHWAYSCANNER_OUTPUTS pathway_lexer.c)
endif(BISONFLEX)


include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(/usr/local/include)

set(FILELIST nupack_invariants.cc physical_spec.cc 
        pathway_input.cc constraint_handler.cc objective_handler.cc
        design_result.cc designer.cc pathway_utils.cc design_spec.cc
        symmetry_calc.cc weight_spec.cc node_spec.cc pair_probabilities.cc
        sequence_spec.cc sequence_utils.cc split_set.cc structure_spec.cc
        complex_spec.cc structure_utils.cc tube_spec.cc
        ${BISON_PATHWAYPARSER_OUTPUTS} ${FLEX_PATHWAYSCANNER_OUTPUTS}
        complex_result.cc complex_spec.cc node_result.cc sequence_state.cc 
        structure_result.cc tube_result.cc physical_result.cc
        equilibrium_concentrations.c utils.c)

add_library(msdesign OBJECT ${FILELIST})

add_executable(multitubedesign multistatedesign.cc $<TARGET_OBJECTS:msdesign>)
target_link_libraries(multitubedesign ${GPERF_PROFILER_LIBRARY} 
        nupackpfunc nupackutils ${MATH_LIB} ${JSONCPP_LIB} ${MPI_C_LIBRARIES} )

add_executable(multitubedefect multistatedefect.cc $<TARGET_OBJECTS:msdesign>)
target_link_libraries(multitubedefect ${GPERF_PROFILER_LIBRARY} 
        nupackpfunc nupackutils ${MATH_LIB} ${JSONCPP_LIB} ${MPI_C_LIBRARIES} )

if(BISONFLEX)
add_dependencies(multitubedesign msdesign_rename)
add_dependencies(multitubedefect msdesign_rename)
endif(BISONFLEX)

install(TARGETS multitubedesign DESTINATION ${BIN_DIR})
install(TARGETS multitubedefect DESTINATION ${BIN_DIR})

