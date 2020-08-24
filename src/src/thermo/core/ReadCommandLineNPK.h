#ifndef __READCOMMANDLINENPK_H__
#define __READCOMMANDLINENPK_H__

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <ctype.h>

#include "shared.h"

#ifdef __cplusplus
extern "C" {
#endif

//get all command line options, including a possible input file
int ReadCommandLineNPK(int nargs, char **args, char *inputfile);

//print out the thermodynamic options
void PrintNupackThermoHelp(void);
//print out the standard options for NUPACK Utilities executables
void PrintNupackUtilitiesHelp(void);

//read in input file
int ReadInputFile( char*, char*, int*, float*, char*, int*);

//get input interactively
void getUserInput(char*, int*,  float*, char*);

//determine if a permutation has a cyclic symmetry
int calculateVPi( int *, int);

//print outputs to the screen
void printInputs( int, char**, const char *, int,  const float*, const char*, char*);

//Prints a header with copyright information
void header( int, char**, char*, char*);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
