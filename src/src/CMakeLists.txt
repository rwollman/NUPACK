add_subdirectory(shared)
add_subdirectory(thermo)
add_subdirectory(design)

install(FILES "nupack.h" DESTINATION ${INCLUDE_INSTALL_LOCATION})
install(FILES "design.h" DESTINATION ${INCLUDE_INSTALL_LOCATION}/nupack)
install(FILES "shared.h" DESTINATION ${INCLUDE_INSTALL_LOCATION}/nupack)
install(FILES "thermo.h" DESTINATION ${INCLUDE_INSTALL_LOCATION}/nupack)
