# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Joseph Zadeh 

# Makefile for single-stranded design functions
# src/design/

NUPACKHOME = ../../..
BIN_DIR = $(NUPACKHOME)/bin
CONSTANTS_DIR=$(NUPACKHOME)/src/shared/constants
LIB_DIR=$(NUPACKHOME)/lib
PFUNCUTILS_DIR=../../thermo/utils
RNG_DIR=$(SHARED_DIR)/rng
SHARED_DIR=$(NUPACKHOME)/src/shared

# Set up some defaults

ifndef CC
  CC=icc
endif
ifndef NUPACK_CFLAGS
  NUPACK_CFLAGS = -std=c99 -O3 -Wall -Wmissing-prototypes -Wmissing-declarations
endif

export NUPACK_SAMPLE
SAMPLE_FLAGS=-DNUPACK_SAMPLE
export SAMPLE_FLAGS
ifndef NUPACK_SELECTION_FLAGS
  NUPACK_SELECTION_FLAGS=${SAMPLE_FLAGS} 
endif

DESIGN_LIB=$(LIB_DIR)/libdesign_pfunc_utils.a
DESIGN_UTILS_DIR = $(NUPACKHOME)/src/design/design_pfunc_utils
LIBRARY_OBJECT=nsStar_psStar.o pairPrStruct.o
PFUNC_LIB=$(LIB_DIR)/libpfunc.a
RNG_LIB = $(LIB_DIR)/libmt19937.a
UTILS_DIR = $(NUPACKHOME)/src/thermo/utils
UTILS_LIB=$(LIB_DIR)/libutils.a

INCLUDES=-I$(CONSTANTS_DIR) -I$(UTILS_DIR) -I$(DESIGN_UTILS_DIR) -I$(RNG_DIR)
LIBS=-L../../../lib $(NUPACK_LIBS)

CFLAGS = $(NUPACK_SELECTION_FLAGS) $(NUPACK_CFLAGS) $(INCLUDES)

all: 	design

design: design.o design_utils.o design_engine.o design_test.o read_command_line.o\
 $(CONSTANTS_DIR)/physical_constants.h  $(CONSTANTS_DIR)/runtime_constants.h \
 $(UTILS_DIR)/DNAExternals.h $(DESIGN_UTILS_DIR)/design_pfunc_utils_header.h \
 $(SHARED_DIR)/utils.o $(PFUNC_LIB) $(UTILS_LIB)
	mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) design.o design_utils.o design_engine.o design_test.o \
	  read_command_line.o -o $(BIN_DIR)/design -lm \
	  $(DESIGN_LIB) $(PFUNC_LIB) $(RNG_LIB) $(UTILS_LIB) $(NUPACK_LIBS) $(LIBS)

clean: 
	rm -f *.o
	rm -f $(BIN_DIR)/design
	rm -f *~
