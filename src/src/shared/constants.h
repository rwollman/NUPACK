#ifndef NUPACK_SHARED_CONSTANTS_H__
#define NUPACK_SHARED_CONSTANTS_H__
/*
  physical_constants.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007
*/

// Fix after vetting against test suite:
// #define kB 0.0019872041 // Boltzmann constant in kcal/mol/K
#define kB 0.00198717 // Boltzmann constant in kcal/mol/K
#define ZERO_C_IN_KELVIN 273.15 // Zero degrees C in Kelvin
#define AVOGADRO 6.022e23 // Avogadro's number

/*
  runtime_constants.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois, 1/2007, except where noted
*/

#define MAXLINE 10000 // Maximum characters in a line
#define NAD_INFINITY 100000 //artificial value for positive infinity
#define INF_CUTOFF 0.9 // fabs(1 - val/NAD_INFINITY) > INF_CUTOFF 
                      // if val is to be considered finite
#define NUM_PRECISION 1e-12 // A small number that's basically zero

//the character to use for comments.  This only affects the output,
//not the input files.
#define COMMENT_STRING "%"

// Constants used in random number generation
// These come from Numerical Recipes in C, 2nd edition, by Press, et al.
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// Error codes
#define ERR_FACTORIAL 65 // Error code for factorial overflow

// Nucleotides
enum BASES {
  BASE_NONE = -1, // No base is possible
  BASE_N = 0, //ACTG
  BASE_A = 1,
  BASE_C = 2,
  BASE_G = 3,
  BASE_T = 4,
  BASE_U = 4,
  BASE_R = 5, // AG
  BASE_M = 6, // AC
  BASE_S = 7, // CG
  BASE_W = 8, // AU
  BASE_K = 9, // GU
  BASE_Y = 10, // CU
  BASE_V = 11, // ACG
  BASE_H = 12, // ACU
  BASE_D = 13, // AGU
  BASE_B = 14, // CGU
  STRAND_PLUS = 15, // Strand break
};


//sets the type of floating point variables
#ifdef USE_DOUBLE

#define DBL_TYPE double
#define EXP_FUNC exp
#define LOG_FUNC log
#define ABS_FUNC fabs

#else

#define DBL_TYPE long double
#define EXP_FUNC expl
#define LOG_FUNC logl
#define ABS_FUNC fabsl

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
#define MAX_FILENAME_LEN 4096
 
//MATCH_PF will make the energy model used in energy calculations
//match the one used in mfe and partition function calculations.
//Otherwise, the energy of multiloops scales with the log of the size,
//rather than linearly.
//Other refinements, such as coaxial stacking, could also be included.
#define MATCH_PF

#define NUPACK_DEF_TEMPERATURE              310.15
#define NUPACK_DEF_SODIUM                   1.0
#define NUPACK_DEF_MAGNESIUM                0.00
#define NUPACK_DEF_MIN_PPAIR                0.00001
#define NUPACK_DEF_M_RESEED                 50
#define NUPACK_DEF_M_BAD                    300.0
#define NUPACK_DEF_M_REOPT                  3
#define NUPACK_DEF_F_SPLIT                  0.99
#define NUPACK_DEF_F_PASSIVE                0.01
#define NUPACK_DEF_F_STRINGENT              0.99
#define NUPACK_DEF_F_REDECOMP               0.03
#define NUPACK_DEF_F_REFOCUS                0.03
#define NUPACK_DEF_GC_INIT_PROB             0.50
#define NUPACK_DEF_MATERIAL                 RNA
#define NUPACK_DEF_DANGLE_TYPE              1
#define NUPACK_DEF_USE_LONG_HELIX           0
#define NUPACK_DEF_SEED                     0
#define NUPACK_DEF_H_SPLIT                  2
#define NUPACK_DEF_N_SPLIT                  12
#define NUPACK_DEF_POPULATION               1
#define NUPACK_DEF_ALLOW_MISMATCH           false
#define NUPACK_DEF_ALLOW_WOBBLE             false
#define NUPACK_DEF_MAX_PRINT_STEPS          1
#define NUPACK_DEF_DISABLE_DEFECT_WEIGHTS   false
#define NUPACK_DEF_DISABLE_FOCUS            false
#define NUPACK_DEF_FORBID_SPLITS            true
#define NUPACK_DEF_REDECOMPOSE              false
#define NUPACK_DEF_WORD_LEN                 4
#define NUPACK_DEF_STOP                     0.01
#define NUPACK_DEF_BONUS_PER_SPLIT          -25

#endif
