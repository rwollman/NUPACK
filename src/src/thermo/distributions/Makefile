# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 

# Makefile for thermo/distributions, for use with NUPACK.
# Requires that src/shared/rng and src/shared have already been built
# (to provide $LIB_DIR/libutils.a and $LIB_DIR/libpfunc.a).

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
INCLUDES=-I$(CONSTANTS_DIR)
LIB_DIR=$(NUPACKHOME)/lib
SHARED_DIR=$(NUPACKHOME)/src/shared
UTILS_LIB=libutils.a

CFLAGS=$(NUPACK_CFLAGS) $(INCLUDES)


all: $(BIN_DIR)/.exist distributions

$(BIN_DIR)/.exist:
	mkdir -p $(BIN_DIR)
	touch $(BIN_DIR)/.exist

distributions: distributions.o ReadCommandLine.o InputFileReader.o \
               OutputWriter.o CalcDist.o DistributionsHeaderFile.h \
               $(CONSTANTS_DIR)/physical_constants.h  \
               $(CONSTANTS_DIR)/runtime_constants.h \
               $(SHARED_DIR)/utils.o $(LIB_DIR)/$(UTILS_LIB)
	$(CC) $(CFLAGS) distributions.o ReadCommandLine.o InputFileReader.o \
              OutputWriter.o CalcDist.o -o $(BIN_DIR)/distributions -lm \
              $(LIB_DIR)/$(UTILS_LIB)

clean: 
	rm -f *.o
	rm -f $(BIN_DIR)/distributions
	rm -f *~
