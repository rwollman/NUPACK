# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 1/2007
#
# Makefile for generic utilities for use with NUPACK.


NUPACKHOME=../..


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

CONSTANTS_DIR=constants
LIB_DIR=$(NUPACKHOME)/lib
UTILS_LIB=$(LIB_DIR)/libutils.a

CFLAGS=$(NUPACK_CFLAGS) -I$(CONSTANTS_DIR)


all: $(LIB_DIR)/.exist library

$(LIB_DIR)/.exist:
	mkdir -p $(LIB_DIR)
	touch $(LIB_DIR)/.exist        

utils:  $(LIB_DIR)/.exist utils.o utilsHeader.h hash.h \
        $(CONSTANTS_DIR)/runtime_constants.h \
        $(CONSTANTS_DIR)/physical_constants.h hash.o
	ar rcs $(UTILS_LIB) utils.o

library: utils.o utilsHeader.h $(CONSTANTS_DIR)/runtime_constants.h \
	$(CONSTANTS_DIR)/physical_constants.h hash.o
	$(AR) rcs $(UTILS_LIB) utils.o hash.o

clean:
	rm -f *o
	rm -f *~
	rm -f $(CONSTANTS_DIR)/*~
	rm -f $(UTILS_LIB)
