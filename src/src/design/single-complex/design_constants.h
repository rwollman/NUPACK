#ifndef DESIGN_CONSTANTS_H
#define DESIGN_CONSTANTS_H

#include "float.h"

/* **************************************************************** */
// Define salt parameters to be defaults.  THIS SHOULD CHANGE WHEN 
// SALT CORRECTIONS ARE IMPLEMENTED IN DESIGN.
//#define SODIUM_CONC 1.0
//#define MAGNESIUM_CONC 0.0
#define USE_LONG_SALT 0
/* **************************************************************** */

#define MAX_FILENAME_SIZE 1000

#define MAX_SEQ_LENGTH 4000
#define MIN_SEQ_LENGTH 0

#define RANDOM_INIT 0
#define SSM_INIT 1
#define CG_INIT 2
#define AU_INIT 3

#define N_OPTIMIZATION 0
#define MFE_OPTIMIZATION 1
#define P_OPTIMIZATION 2

#define H_MIN_RNA 2
#define H_MIN_DNA 3

#define U_MIN 20

#define N_RATIO_DEFAULT 0.01  

#define M_REOPT_DEFAULT 10
#define M_LEAFOPT_DEFAULT 3

#define M_UNFAVORABLE 4.0

#define M_EVAL 5000
#define B_ACCEPT 0.2

#define MFE_EPS 0.1
#define SSM_RETRY 10

//#define DEBUG
//#define DEBUG_PRINT
#define LEAFCORRECTION
//#define FREEZEBASES
//#define EARLY_TERMINATE
//#define FINDBUG



//#ifdef LEAFCORRECTION
// setting to 4000 as not to change previous results
// can change after paper submission
#define MAX_NODES 4000
//#endif 

#endif
