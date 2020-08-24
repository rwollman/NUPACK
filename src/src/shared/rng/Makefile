# This Makefile is part of the NUPACK software suite.
# Make random number generator library from src/shared/rng for use with NUPACK.


NUPACKHOME=../../..

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

CFLAGS=$(NUPACK_CFLAGS)
LIB_DIR=$(NUPACKHOME)/lib
LIBRARY_OBJECT=mt19937ar.o
PFUNC_LIB=$(LIB_DIR)/libmt19937.a


all: $(PFUNC_LIB)

$(PFUNC_LIB): $(LIBRARY_OBJECT)
	mkdir -p $(LIB_DIR)
	$(AR) rcs $(PFUNC_LIB) $(LIBRARY_OBJECT)

clean:
	rm -f *.o rng *.a *~ $(PFUNC_LIB)

#TEST_LIBS=-L$(LIB_DIR) -lmt19937 -L. -lgsl -lgslcblas -lm
#test: 	mtTest.o
#	$(CC) $(CFLAGS) mtTest.o -o rng $(TEST_LIBS)
