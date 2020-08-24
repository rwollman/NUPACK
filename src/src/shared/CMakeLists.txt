configure_file(externals.c.in "${CMAKE_CURRENT_BINARY_DIR}/externals.c")

add_library(nupackutils hash.c mt19937ar.c functions.c "${CMAKE_CURRENT_BINARY_DIR}/externals.c")

install(TARGETS nupackutils DESTINATION ${LIBRARY_INSTALL_LOCATION})

file(GLOB SHARED_INCLUDES *.h)
install(FILES ${SHARED_INCLUDES} DESTINATION ${INCLUDE_INSTALL_LOCATION}/nupack/shared)
