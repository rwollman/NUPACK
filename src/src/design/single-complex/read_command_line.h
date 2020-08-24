#ifndef RCL_H
#define RCL_H


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <ctype.h>

#include "pfuncUtilsHeader.h"
#include "physical_constants.h"
#include "runtime_constants.h"
#include "DNAExternals.h"
#include "design_constants.h"

void displayHelp(void);
void makeUCase(char *line);
void print_design_help(void);
void readCommandLine(int nargs, char **args, char* inputPrefix, char *psFile, int *initMode, 
			int *bypassDesign, int *bypassHeirarchy, int *bypassGuidance, int *designMode, 
		     int *loadSeeds, int *mReopt, int * mLeafopt, DBL_TYPE *nRatio, int *loadInit,
                     int *hMin, DBL_TYPE *pcut, int *quick, int * output_init,
                     int * output_seed, int * output_ppairs, int * output_json);
#endif
