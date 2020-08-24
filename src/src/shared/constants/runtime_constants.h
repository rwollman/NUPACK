/*
  runtime_constants.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois, 1/2007, except where noted
*/

#define MAXLINE 10000 // Maximum characters in a line
#define NAD_INFINITY 100000 //artificial value for positive infinity
#define INF_CUTOFF 0.9 // fabs(1 - val/NAD_INTINITY) > INF_CUTOFF 
                      // if val is to be considered finite
#define NUM_PRECISION 1e-12 // A small number that's basically zero

//the character to use for comments.  This only affects the output,
//not the input files.
#define COMMENT_STRING "%"

// Version of NUPACK this is.
#define VERSION "3.0"

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

