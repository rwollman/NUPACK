# This function is to add executables that depend on the partition function utilities
function(np_add_executable NAME)
    add_executable(${NAME} ${ARGN})
    target_link_libraries(${NAME} 
                          nupackpfunc 
                          nupackutils 
                          ${CLANG_ASAN} 
                          ${MATH_LIB}
                          ${MPI_C_LIBRARIES}
    )
    install(TARGETS ${NAME} RUNTIME DESTINATION ${BIN_DIR})
endfunction(np_add_executable)

add_subdirectory(concentrations)
add_subdirectory(core)
add_subdirectory(complexes)
add_subdirectory(basics)
add_subdirectory(distributions)
add_dir_if_exists(centroid)

install_include_tree(centroid)
install_include_tree(complexes)
install_include_tree(concentrations)
install_include_tree(distributions)
install_include_tree(core)
