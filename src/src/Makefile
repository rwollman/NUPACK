# This Makefile is part of the NUPACK software suite
# Copyright (c) 2007 Caltech. All rights reserved.
# Coded by: Justin Bois 1/2007
# Makefile for making all of NUPACK

# NOTES:
# *SAMPLE* refs: thermo/basics  thermo/utils  design/single-complex
#
# --cplusplus refs: shared/rng               shared  
#                   thermo/complexes         thermo/distributions
#                   thermo/concentrations    thermo/utils
#                   design/design_pfunc_utils  

NUPACKHOME=..
LIB_DIR=$(NUPACKHOME)/lib

# Set up some defaults

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

export NUPACK_SAMPLE
SAMPLE_FLAGS=-DNUPACK_SAMPLE
export SAMPLE_FLAGS
ifndef NUPACK_SELECTION_FLAGS
  NUPACK_SELECTION_FLAGS=${SAMPLE_FLAGS} 
endif


export CC
export CXX
export NUPACK_CFLAGS
export NUPACK_CXXFLAGS
export NUPACK_FLAGS
export NUPACK_LIBS
export NUPACK_SELECTION_FLAGS

CFLAGS=$(NUPACK_CFLAGS) $(NUPACK_SELECTION_FLAGS)


# Targets
.PHONY: all
all: subdirs

# Subdirectories
SUBDIRS=shared/rng shared thermo/utils  thermo/complexes \
	thermo/distributions thermo/concentrations  \
	thermo/basics \
	design/design_pfunc_utils  \
	design/single-complex 

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	$(MAKE) -C ../lib clean
	$(MAKE) -C shared clean
	$(MAKE) -C shared/rng clean
	$(MAKE) -C thermo/basics clean
	$(MAKE) -C thermo/complexes clean
	$(MAKE) -C thermo/concentrations clean
	$(MAKE) -C thermo/distributions clean
	$(MAKE) -C thermo/utils clean
	$(MAKE) -C design/design_pfunc_utils clean
	$(MAKE) -C design/single-complex clean
