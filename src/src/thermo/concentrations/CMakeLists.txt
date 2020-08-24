
add_library(nupackconc CalcConc.c)

np_add_executable(concentrations concentrations.c ReadCommandLine.c InputFileReader.c
    OutputWriter.c CalcConc.c FracPair.c)

install(TARGETS nupackconc DESTINATION ${LIBRARY_INSTALL_LOCATION})
