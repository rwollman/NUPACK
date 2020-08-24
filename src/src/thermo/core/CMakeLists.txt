file(GLOB FILELIST *.c)
if(NOT MPI)
    list(REMOVE_ITEM FILELIST "${CMAKE_CURRENT_SOURCE_DIR}/nupack_mpi.c")
endif()

add_library(nupackpfunc ${FILELIST})
install(TARGETS nupackpfunc DESTINATION ${LIBRARY_INSTALL_LOCATION})
