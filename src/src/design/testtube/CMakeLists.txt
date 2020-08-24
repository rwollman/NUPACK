
cmake_minimum_required(VERSION 2.6)

set(BISON_TUBEPARSER_OUTPUTS npparser.c)
set(FLEX_TUBESCANNER_OUTPUTS nplexer.c)

if(BISONFLEX)
    BISON_TARGET(TUBEPARSER tubeparser.y ${CMAKE_CURRENT_SOURCE_DIR}/npparser.c)
    FLEX_TARGET(TUBESCANNER tubeparser.l ${CMAKE_CURRENT_SOURCE_DIR}/nplexer.c COMPILE_FLAGS --header-file=${CMAKE_CURRENT_SOURCE_DIR}/nplexer.h)
    ADD_FLEX_BISON_DEPENDENCY(TUBESCANNER TUBEPARSER)
    
        
    string(REGEX REPLACE "/" "_" PATH_STRING ${CMAKE_CURRENT_SOURCE_DIR})
    string(TOUPPER ${PATH_STRING} PATH_STRING)
    
    add_custom_command(OUTPUT mock.out1
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${FLEX_TUBESCANNER_OUTPUTS}
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${CMAKE_CURRENT_SOURCE_DIR}/nplexer.h
        DEPENDS ${FLEX_TUBESCANNER_OUTPUTS} ${CMAKE_CURRENT_SOURCE_DIR}/nplexer.h
        )
    add_custom_command(OUTPUT mock.out2
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${BISON_TUBEPARSER_OUTPUT_SOURCE}
        COMMAND sed -e \"s:${CMAKE_CURRENT_SOURCE_DIR}/::g\" -i \"\" ${BISON_TUBEPARSER_OUTPUT_HEADER}
        COMMAND sed -e \"s:${PATH_STRING}::g\" -i \"\" ${BISON_TUBEPARSER_OUTPUT_SOURCE}
        COMMAND sed -e \"s:${PATH_STRING}::g\" -i \"\" ${BISON_TUBEPARSER_OUTPUT_HEADER}
        DEPENDS ${BISON_TUBEPARSER_OUTPUT_SOURCE} ${BISON_TUBEPARSER_OUTPUT_HEADER}
        )
    add_custom_target(tubedesign_rename ALL DEPENDS mock.out1 mock.out2)
endif(BISONFLEX)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(tubedesign main.c pathway_design.c pathway_utils.c 
    pathway_output.c equilibrium_concentrations.c utils.c pathway_input.c
    ${BISON_TUBEPARSER_OUTPUTS} ${FLEX_TUBESCANNER_OUTPUTS})

add_executable(tubedefect defect_main.c pathway_design.c pathway_utils.c 
    pathway_output.c equilibrium_concentrations.c utils.c pathway_input.c
    ${BISON_TUBEPARSER_OUTPUTS} ${FLEX_TUBESCANNER_OUTPUTS})

# add_executable(decomp decomp.c pathway_design.c pathway_utils.c 
#     pathway_output.c equilibrium_concentrations.c utils.c pathway_input.c
#     ${BISON_TUBEPARSER_OUTPUTS} ${FLEX_TUBESCANNER_OUTPUTS})

target_link_libraries(tubedesign nupackpfunc nupackutils ${MATH_LIB} ${CLANG_ASAN} ${MPI_C_LIBRARIES})
target_link_libraries(tubedefect nupackpfunc nupackutils ${MATH_LIB} ${CLANG_ASAN} ${MPI_C_LIBRARIES})

if(BISONFLEX)
add_dependencies(tubedesign tubedesign_rename)
add_dependencies(tubedefect tubedesign_rename)
endif(BISONFLEX)

# target_link_libraries(decomp nupackpfunc nupackutils  ${CLANG_ASAN})

# add_executable(parsetest parsetest.c pathway_input.c 
# ${BISON_TUBEPARSER_OUTPUTS} ${FLEX_TUBESCANNER_OUTPUTS})

install(TARGETS tubedesign DESTINATION ${BIN_DIR})
install(TARGETS tubedefect DESTINATION ${BIN_DIR})
