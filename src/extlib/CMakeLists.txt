if (KINETICS)
add_dir_if_exists(odepack)
add_dir_if_exists(Multistrand_src)
install_include_tree(Multistrand_src)
endif (KINETICS)
install_include_tree(tclap)
install_include_tree(jsoncpp)

add_library(jsoncpp jsoncpp/jsoncpp.cpp)
