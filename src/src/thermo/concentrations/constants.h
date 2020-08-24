/*
  concentrations.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.

  Header file with global variables and function prototypes for use
  with concentrations.c Contains trust region parameters, among
  others.
*/

#ifndef NUPACK_THERMO_CONCENTRATIONS_H
#define NUPACK_THERMO_CONCENTRATIONS_H

#include "shared.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <getopt.h>

// Constants used in trust region
#define TRUST_REGION_DELTABAR 1000.0 // Maximal size of trust region
#define TRUST_REGION_ETA 0.125 // Decision criterion for trust region

#define MAXLOGX 250 // Maximum logarithm of a concentration (prevents overflow)

// Where the help file to print is
#define NUPACK_THERMO_CONCENTRATIONS_HELP_FILE "src/thermo/concentrations/Concentrations.help"

// Error codes
#define ERR_NOCONVERGE 1 // Error code for failure to converge for method = 1
#define ERR_OVERFLOW 2 // Overflow in calculation of mole fractions
#define ERR_NOINPUT 3 // No prefix for a filename given
#define ERR_HELP 4 // User chose to display help
#define ERR_HELP_OPEN 5 // Error in opening README file to print help
#define ERR_CON 6 // Error opening .con input file
#define ERR_CX 7 // Error in opening .cx file
#define ERR_NONSEQUENTIAL 8 // CX file has nonsequential complex IDs
#define ERR_NOPERMS 9 // Free energies are too high to allow calcs with perms
#define ERR_BADROWINP 10 // Bad row in stoich. matrix in input file
#define ERR_DUPROWA 11 // Row of A is duplicated (indicates bug, not input err)
#define ERR_LOG 12 // Error opening log file
#define ERR_EQ 13 // Error opening .eq file
#define ERR_FPAIRS 14 // Error opening .fpairs file
#define ERR_NOSEQEQ 15 // No sequence information in comments of eq file
#define ERR_PAIRSFILE 16 // Error opening pairs file

#endif // NUPACK_THERMO_CONCENTRATIONS_H
