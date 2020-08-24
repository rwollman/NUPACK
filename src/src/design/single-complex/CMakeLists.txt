add_executable(complexdesign  design.c design_engine.c design_utils.c 
    read_command_line.c nsStar_psStar pairPrStruct)

target_link_libraries(complexdesign nupackpfunc nupackutils m ${MPI_C_LIBRARIES})

install(TARGETS complexdesign DESTINATION ${BIN_DIR})
