add_subdirectory(single-complex)
add_subdirectory(testtube)
add_dir_if_exists(multistate)

if (${JSONCPP_FOUND})
    add_dir_if_exists(pathway_design)
    install_include_tree(pathway_design)
endif()

install_include_tree(single-complex)
install_include_tree(testtube)
install_include_tree(multistate)
