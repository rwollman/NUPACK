/* 
  pfuncUtilsConstants.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2001, Justin Bois 1/2007

  Useful constants for running partition function applications.
*/

#ifdef __cplusplus
extern "C" {
#endif 

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "physical_constants.h"
#include "runtime_constants.h"
#ifdef GC_DEBUG

#include "gc.h"
#define malloc(n) GC_MALLOC(n)
#define calloc(m,n) GC_MALLOC((m)*(n))
#define free(p) GC_FREE(p)
#define realloc(p,n) GC_REALLOC((p),(n))
#define CHECK_LEAKS() GC_gcollect()
#endif

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#ifdef NEDMALLOC
#include "nedmalloc.h"
#endif



//sets the type of floating point variables
#ifdef USE_DOUBLE

#define DBL_TYPE double
#define EXP_FUNC exp
#define LOG_FUNC log

#else

#define DBL_TYPE long double
#define EXP_FUNC expl
#define LOG_FUNC logl

#endif

//Minimum difference between energies before two are considered identical
#define ENERGY_TOLERANCE 0.0001


//max error in the bits of precision.  Used during pair probability
//calculations (where subtraction occurs) Setting this to zero can
//significantly slow down pair probability calculations.
#define MAXPRECERR 24 //max error in bits of precision

//Maximum seqeuence length
#define MAXSEQLENGTH 10000

//maximum # of strands in a complex
#define MAXSTRANDS 2000
 
//MATCH_PF will make the energy model used in energy calculations
//match the one used in mfe and partition function calculations.
//Otherwise, the energy of multiloops scales with the log of the size,
//rather than linearly.
//Other refinements, such as coaxial stacking, could also be included.
#define MATCH_PF

//Including NOGU will disallow all wobble pairs
//#define NOGU

//Including STRUCTURE_WARNINGS will produce error messages whenever
//a disconnected or illegal structure is evaluated, rather than just
//returning an "infinite" energy.
//#define STRUCTURE_WARNINGS

//Including DEBUG will cause various intermediate values to be printed during backtracking (backtrack.c)
//#define DEBUG

//Including FILEOUTPUT will cause pair probabilities to be spit out to a file named Pb_N#.txt, where
//the # is replaced with the time-complexity of the algorithm.

//#define FILEOUTPUT

//Including PRINTRESULTSONLY will limit output to the screen.  This needs to be updated.
//#define PRINTRESULTSONLY

//Including NODANGLES will set all dangle energies to zero.  Good for debugging purposes.
//#define NODANGLES

/* No changes below this line! */

//some physical constants
#define BASE_N 0 //ACTG
#define BASE_A 1
#define BASE_C 2
#define BASE_G 3
#define BASE_T 4
#define BASE_U 4
#define BASE_R 5 // AG
#define BASE_M 6 // AC
#define BASE_S 7 // CG
#define BASE_W 8 // AU
#define BASE_K 9 // GU
#define BASE_Y 10 // CU
#define BASE_V 11 // ACG
#define BASE_H 12 // ACU
#define BASE_D 13 // AGU
#define BASE_B 14 // CGU
#define STRAND_PLUS 15 // Strand break

#endif

#ifdef __cplusplus
}
#endif 
