# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 

# Makefile for concentrations for use in NUPACK
# src/thermo/concentrations

NUPACKHOME=../../..

ifndef CC
  CC=icc
endif
ifndef CXX
  CXX=icpc
endif
ifndef NUPACK_CFLAGS
  NUPACK_CFLAGS=-std=c99 -O3 -Wall -Wmissing-prototypes -Wmissing-declarations
endif
ifndef NUPACK_CXXFLAGS
  NUPACK_CXXFLAGS=-Wall -Wmissing-prototypes -Wmissing-declarations
endif

BIN_DIR=$(NUPACKHOME)/bin
CONSTANTS_DIR=$(NUPACKHOME)/src/shared/constants
LIBS=-L../../../lib $(NUPACK_LIBS)
PFUNC_LIB=$(NUPACKHOME)/lib/libpfunc.a
RNG_LIB=$(NUPACKHOME)/lib/libmt19937.a
SHARED_DIR=$(NUPACKHOME)/src/shared
UTILS_LIB=$(NUPACKHOME)/lib/libutils.a

INCLUDES=-I$(CONSTANTS_DIR) -I$(SHARED_DIR)

CFLAGS=$(NUPACK_CFLAGS) $(INCLUDES)


all: concentrations library

concentrations: concentrations.o ReadCommandLine.o InputFileReader.o \
                OutputWriter.o CalcConc.o FracPair.o \
                ConcentrationsHeaderFile.h \
                $(CONSTANTS_DIR)/physical_constants.h \
                $(CONSTANTS_DIR)/runtime_constants.h \
                $(SHARED_DIR)/utils.o $(UTILS_LIB)
	mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) concentrations.o ReadCommandLine.o InputFileReader.o \
              OutputWriter.o CalcConc.o FracPair.o \
              -o $(BIN_DIR)/concentrations $(LIBS) -lm $(UTILS_LIB) $(RNG_LIB)

LIBRARY_OBJECT=CalcConc.o
.PHONY: library
library: CalcConc.o
	mkdir -p $(NUPACKHOME)/lib
	ar rcs $(PFUNC_LIB) $(LIBRARY_OBJECT)

clean: 
	rm -f *.o
	rm -f $(BIN_DIR)/concentrations
	rm -f *~
