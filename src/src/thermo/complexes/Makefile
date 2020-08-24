# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 1/2007

# Makefile for thermo/complexes, for use with NUPACK.

NUPACKHOME=../../..

BIN_DIR=$(NUPACKHOME)/bin
CONSTANTS_DIR=$(NUPACKHOME)/src/shared/constants
LIB_DIR=$(NUPACKHOME)/lib
PFUNCUTILS_DIR=../utils
PFUNC_LIB=libpfunc.a
SHARED_DIR=$(NUPACKHOME)/src/shared
UTILS_LIB=libutils.a

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

INCLUDES=-I$(PFUNCUTILS_DIR) -I$(SHARED_DIR) -I$(CONSTANTS_DIR)
LIBS=$(NUPACK_LIBS) -L$(LIB_DIR) -lpfunc -lutils -lmt19937 -lm

CFLAGS=$(NUPACK_CFLAGS) $(INCLUDES)


COMPLEXES_DEPEND=$(SHARED_DIR)/utils.o \
 $(SHARED_DIR)/utilsHeader.h\
 $(PFUNCUTILS_DIR)/pfuncUtilsHeader.h\
 $(PFUNCUTILS_DIR)/backtrack.o \
 $(PFUNCUTILS_DIR)/CalculateEnergy.o \
 $(PFUNCUTILS_DIR)/ene.o \
 $(PFUNCUTILS_DIR)/GetEnergy.o \
 $(PFUNCUTILS_DIR)/init.o \
 $(PFUNCUTILS_DIR)/mfeUtils.o \
 $(PFUNCUTILS_DIR)/min.o \
 $(PFUNCUTILS_DIR)/pairsPr.o \
 $(PFUNCUTILS_DIR)/pf.o \
 $(PFUNCUTILS_DIR)/pfuncUtils.o \
 $(PFUNCUTILS_DIR)/pknots.o \
 $(PFUNCUTILS_DIR)/ReadCommandLineNPK.o \
 $(PFUNCUTILS_DIR)/sumexp.o \
 $(PFUNCUTILS_DIR)/sumexp_pk.o \
 $(PFUNCUTILS_DIR)/pfuncUtilsConstants.h\
 $(PFUNCUTILS_DIR)/DNAExternals.h\
 $(PFUNCUTILS_DIR)/DNAGlobals.o\
 $(CONSTANTS_DIR)/runtime_constants.h\
 $(CONSTANTS_DIR)/physical_constants.h

complexes: $(LIB_DIR)/$(UTILS_LIB) $(LIB_DIR)/$(PFUNC_LIB) complexesHeader.h \
           complexesStructs.h complexes.o complexesUtils.o permBG.o \
           ReadCommandLine.o $(COMPLEXES_DEPEND)
	-mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) complexes.o complexesUtils.o permBG.o \
	ReadCommandLine.o -o $(BIN_DIR)/complexes $(LIBS)

clean: 
	rm -f *.o
	rm -f $(BIN_DIR)/complexes
	rm -f *~
