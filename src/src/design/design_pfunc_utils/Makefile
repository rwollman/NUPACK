# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 1/2007

# Makefile for complexes for use with NUPACK
# src/design/design_pfunc_utils
# (Create the design_pfunc_utils library.)

NUPACKHOME=../../..
LIB_DIR=$(NUPACKHOME)/lib

ifndef CC
  CC=icc
endif
ifndef CXX
  CXX=icpc
endif
ifndef NUPACK_CFLAGS
  NUPACK_CFLAGS = -std=c99 -O3 -Wall -Wmissing-prototypes -Wmissing-declarations
endif
ifndef NUPACK_CXXFLAGS
  NUPACK_CXXFLAGS = -Wall -Wmissing-prototypes -Wmissing-declarations
endif

CONSTANTS_DIR=$(NUPACKHOME)/src/shared/constants
DESIGN_LIB=$(LIB_DIR)/libdesign_pfunc_utils.a
LIBRARY_OBJECT=nsStar_psStar.o pairPrStruct.o
PFUNC_LIB=$(LIB_DIR)/libpfunc.a
PFUNCUTILS_DIR=../../thermo/utils
SHARED_DIR=$(NUPACKHOME)/src/shared
UTILS_LIB=$(LIB_DIR)/libutils.a

INCLUDES=-I$(PFUNCUTILS_DIR) -I$(CONSTANTS_DIR) -I${SHARED_DIR}

CFLAGS=$(NUPACK_CFLAGS) $(INCLUDES)


DESIGN_PFUNC_UTILS_DEPEND=$(SHARED_DIR)/utils.o \
 $(SHARED_DIR)/utilsHeader.h \
 $(PFUNCUTILS_DIR)/pfuncUtilsHeader.h \
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
 $(PFUNCUTILS_DIR)/pfuncUtilsConstants.h \
 $(PFUNCUTILS_DIR)/DNAExternals.h \
 $(PFUNCUTILS_DIR)/DNAGlobals.o \
 $(CONSTANTS_DIR)/runtime_constants.h \
 $(CONSTANTS_DIR)/physical_constants.h
 

all: $(DESIGN_LIB)

$(DESIGN_LIB): $(LIBRARY_OBJECT) design_pfunc_utils_header.h \
                  $(DESIGN_PFUNC_UTILS_DEPEND) $(UTILS_LIB) $(PFUNC_LIB)
	mkdir -p $(LIB_DIR)
	$(AR) rcs $(DESIGN_LIB) $(LIBRARY_OBJECT)

clean: 
	rm -f *.o
	rm -f $(DESIGN_LIBR)
	rm -f *~
