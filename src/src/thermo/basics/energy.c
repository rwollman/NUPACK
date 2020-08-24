/*
  energy.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 3/2006 and Justin Bois 1/2007

  ENERGY.C

  This code computes the energy for a specified secondary structure
  for a given sequence or set of sequences.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <thermo/core.h>

/* ************************************************ */

int main( int argc, char *argv[] ) {
  
  char seq[ MAXSEQLENGTH], parens[ MAXSEQLENGTH];
  int seqNum[MAXSEQLENGTH+1];
  int thepairs[MAXSEQLENGTH];
  DBL_TYPE ene;
  int vs, tmpLength;

  char inputFile[ MAXLINE];
  int inputFileSpecified;

  strcpy( inputFile, "");
  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  if(NupackShowHelp) {
    printf("Usage: energy [OPTIONS] PREFIX\n");
    printf("Calculate the free energy of the input sequence and structure\n");
    printf("Example: energy -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }
  if( !inputFileSpecified || 
      !ReadInputFile( inputFile, seq, &vs, NULL, parens, thepairs) ) {
       if (inputFileSpecified == 0) getUserInput( seq, &vs, NULL, parens);
       else abort();
  }

  header( argc, argv, "energy", "screen");

  tmpLength = strlen( seq);
  convertSeq(seq, seqNum, tmpLength);

  
  if (parens[0] == '\0') {
    ene = naEnergyPairsOrParensFullWithSym( thepairs, NULL, seqNum, DNARNACOUNT, DANGLETYPE, 
					    TEMP_K - ZERO_C_IN_KELVIN, vs,
					    SODIUM_CONC, MAGNESIUM_CONC, 
					    USE_LONG_HELIX_FOR_SALT_CORRECTION);
  }
  else {
    ene = naEnergyPairsOrParensFullWithSym( NULL, parens, seqNum, DNARNACOUNT, DANGLETYPE, 
					    TEMP_K - ZERO_C_IN_KELVIN, vs,
					    SODIUM_CONC, MAGNESIUM_CONC, 
					    USE_LONG_HELIX_FOR_SALT_CORRECTION);
  }

  // Check to see if the results is close to NAD_INFINITY and report
  // error if it is
  if (ABS_FUNC(1.0 - ene/NAD_INFINITY) < INF_CUTOFF) {
    printf("\n\n*** Error: invalid base pair(s) or disconnected complex. Check your inputs. ***\n\n");
    return 0;
  }

  // Print header and inputs to screen
  header( argc, argv, "energy", "screen");	
  printInputs( argc, argv, seq, vs, NULL, parens, "screen");
  
  printf("%s\n%s Energy (kcal/mol):\n",COMMENT_STRING,COMMENT_STRING);
  printf("%16.14Lf\n", (long double )ene);
	
  return 0;
}
/* ****** */




