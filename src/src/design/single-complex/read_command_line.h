#ifndef RCL_H
#define RCL_H


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <ctype.h>

#include "design_utils.h"
#include <shared.h>
#include <thermo.h>


void displayHelp(void);
void makeUCase(char *line);
void readCommandLine(int nargs, char **args, char* inputPrefix, char *psFile, int *initMode, 
			int *bypassDesign, int *bypassHeirarchy, int *bypassGuidance, int *designMode, 
		     int *loadSeeds, int *mReopt, int * mLeafopt, DBL_TYPE *nRatio, int *loadInit,
                     int *hMin, DBL_TYPE *pcut, int *quick, int * output_init,
                     int * output_seed, int * output_ppairs, int * output_json);
#endif
