# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 

# Makefile for utility functions used in calculations
# based off of partition functions for use in NUPACK.
# Builds the $NUPACKHOME/lib 'pfunc' library.
# src/thermo/utils


NUPACKHOME=../../..

LIB_DIR=$(NUPACKHOME)/lib
SHARED_DIR=$(NUPACKHOME)/src/shared

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
  NUPACK_CXXFLAGS = -Wall -Wmissing-prototypes -Wmissing-declarations
endif

export NUPACK_SAMPLE
SAMPLE_FLAGS=-DNUPACK_SAMPLE
ifndef NUPACK_SELECTION_FLAGS
  NUPACK_SELECTION_FLAGS=$(SAMPLE_FLAGS)
endif

MKL_LIBS=
ifdef USE_MKL
  ifndef MKL_HOME
    MKL_HOME=$(HOME)/intel/mkl/10.0.1.014
  endif
  MKL_LIBS=-L$(MKL_HOME)/lib -lmkl -lguide -lpthread
  MKL_INCLUDES=-I$(MKL_HOME)/include
endif


CONSTANTS_DIR=$(SHARED_DIR)/constants
INCLUDES= -I$(SHARED_DIR) -I$(CONSTANTS_DIR)
PFUNC_LIB=$(LIB_DIR)/libpfunc.a

CFLAGS=$(NUPACK_CFLAGS) $(NUPACK_SELECTION_FLAGS) $(INCLUDES)

export CC
export MKL_HOME
export MKL_LIBS
export NUPACK_CFLAGS
export NUPACK_LIBS


LIBRARY_OBJECT=backtrack.o CalculateEnergy.o ene.o GetEnergy.o \
 init.o mfeUtils.o min.o nsStar.o pairsPr.o pf.o pfuncUtils.o \
 pknots.o sumexp.o sumexp_pk.o ReadCommandLineNPK.o DNAGlobals.o \
 $(SHARED_DIR)/utils.o

all: $(PFUNC_LIB)

$(PFUNC_LIB): $(LIBRARY_OBJECT) pfuncUtilsHeader.h \
	$(CONSTANTS_DIR)/runtime_constants.h \
	$(CONSTANTS_DIR)/runtime_constants.h DNAExternals.h \
	$(LIB_DIR)/libutils.a
	-mkdir -p $(LIB_DIR)
	$(AR) rcs $(PFUNC_LIB) $(LIBRARY_OBJECT)

clean:
	rm -f *o
	rm -f $(PFUNC_LIB)
	rm -f *~ 
