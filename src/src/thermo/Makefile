# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 1/2007

# Makefile for the thermo directory of the NUPACK distribution.
# src/thermo


NUPACKHOME=../..
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
  NUPACK_CXXFLAGS=-Wall -Wmissing-prototypes -Wmissing-declarations
endif

export NUPACK_SAMPLE
SAMPLE_FLAGS=-DNUPACK_SAMPLE
ifndef NUPACK_SELECTION_FLAGS
  NUPACK_SELECTION_FLAGS=-DNUPACK_SAMPLE
endif

INCLUDES= -I$(SHARED_DIR) -I$(SHARED_DIR)/constants

CFLAGS=$(NUPACK_CFLAGS) $(NUPACK_SELECTION_FLAGS) $(INCLUDES)


# 13Jan2015 http://en.wikipedia.org/wiki/Math_Kernel_Library:
# Intel Math Kernel Library (Intel MKL) is a library of optimized math routines 
# for science, engineering, and financial applications. Core math functions
#  include BLAS, LAPACK, ScaLAPACK, sparse solvers, fast Fourier transforms, 
# and vector math. The routines in MKL are hand-optimized specifically for 
# Intel processors.
MKL_LIBS=
#ifdef USE_MKL
  ifndef MKL_HOME
    MKL_HOME=$(HOME)/intel/mkl/10.0.1.014
  endif
  MKL_LIBS=-L$(MKL_HOME)/lib -lmkl -lguide -lpthread
  MKL_INCLUDES=-I$(MKL_HOME)/include
#endif

export CC
export MKL_HOME
export MKL_LIBS
export NUPACK_CFLAGS
export NUPACK_LIBS


# Targets
all: $(MALLOC_LIB) subdirs

# S
SUBDIRS = ../shared/rng ../shared utils  complexes \
          distributions concentrations  \
          basics

#$(MALLOC_LIB):
#	$(MAKE) -C ../ptmalloc3
#	cp -f ../ptmalloc3/libptmalloc3.a ../../lib

.PHONY: subdirs $(SUBDIRS)
               
subdirs: $(SUBDIRS)
                         
$(SUBDIRS):
	$(MAKE) -C $@ 

clean: 
	$(MAKE) -C ../../lib clean
	$(MAKE) -C ../shared clean
	$(MAKE) -C utils clean
	$(MAKE) -C basics clean
	$(MAKE) -C complexes clean
	$(MAKE) -C distributions clean
	$(MAKE) -C concentrations clean
