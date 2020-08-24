/*
  DistributionsHeaderFile.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  Header file for use with distributions.c and related
  functions. Contains error messages and function prototypes, among
  other things.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <getopt.h>// Takes options from the command line

#include "shared.h"

// Where the help file to print is
#define DISTRIBUTIONS_HELP_FILE "src/thermo/distributions/distributions.help"

// Error codes
#define ERR_NOINPUT 2 // No prefix for a filename given
#define ERR_HELP 3 // User chose to display help
#define ERR_HELP_OPEN 4 // Error in opening README file to print help
#define ERR_COUNT 5 // Error opening .count input file
#define ERR_CX 6 // Error in opening .cx file
#define ERR_NONSEQUENTIAL 7 // CX file has nonsequential complex IDs
#define ERR_NOPERMS 8 // Free energies are too high to allow calcs with perms
#define ERR_BADROWINP 9 // Bad row in stoich. matrix in input file
#define ERR_DUPROWA 10 // Row of A is duplicated (indicates bug, not input err)
#define ERR_LOG 11 // Error opening log file
#define ERR_DIST 12 // Error opening .dist file
#define ERR_LAMBDATOOBIG 13 // Too many entries in lambda for method = 2
#define ERR_QBOXTOOBIG 14 // Overflow error in calculation of Q_{box}
#define ERR_LAMBDA 15 // Error opening .lam file
#define ERR_NONINTEGER 16 // Non integer input for initial single-strand count
#define ERR_NEGATIVEPROB 17 // Negative probability encountered

