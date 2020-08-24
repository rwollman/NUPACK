/*
  prob.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  This program computes the equilibrium probability that an inputted
  ordered complex will have a speficied secondary structure at
  equilibrium.  For use with NUPACK.
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
  int thepairs[MAXSEQLENGTH];
  int seqNum[ MAXSEQLENGTH+1];
  DBL_TYPE ene;
  DBL_TYPE pf;
  int vs;
  int tmpLength;
  int complexity;
  char inputFile[ MAXLINE];
  int inputFileSpecified;

  strcpy( inputFile, "");

  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  if(NupackShowHelp) {
    printf("Usage: prob [OPTIONS] PREFIX\n");
    printf("Calculate the equilibrium probability of the input structure\n");
    printf("Example: prob -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  header( argc, argv, "prob", "screen");

  if( inputFileSpecified == 0 ||
     !ReadInputFile( inputFile, seq, &vs, NULL, parens, thepairs) ) {
    if (inputFileSpecified == 0) getUserInput( seq, &vs, NULL, parens);
    else abort();
  }

  printInputs( argc, argv, seq, vs, NULL, parens, "screen");

  tmpLength = strlen( seq);
  convertSeq(seq, seqNum, tmpLength);

  if (parens[0] == '\0') {
    ene = naEnergyPairsOrParensFullWithSym( thepairs, NULL, seqNum, DNARNACOUNT, DANGLETYPE,
                                            TEMP_K - ZERO_C_IN_KELVIN, vs, SODIUM_CONC,
                                            MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
  }
  else {
    ene = naEnergyPairsOrParensFullWithSym( NULL, parens, seqNum, DNARNACOUNT, DANGLETYPE,
                                            TEMP_K - ZERO_C_IN_KELVIN, vs,  SODIUM_CONC,
                                            MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
  }

  if( !DO_PSEUDOKNOTS ) {
    complexity = 3;
  }
  else {
    complexity = 5;
  }

  //calculate partition function, without pairs info
  pf = pfuncFullWithSym(seqNum, complexity, DNARNACOUNT, DANGLETYPE,
                        TEMP_K - ZERO_C_IN_KELVIN, 0, vs, SODIUM_CONC,
                        MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);


  printf("%s\n%s Probability:\n", COMMENT_STRING,COMMENT_STRING);
  if(!NUPACK_VALIDATE) {
    printf("%4.3Le\n", (long double) (EXP_FUNC(-ene/(kB*TEMP_K))/pf) );
  } else {
    printf("%.14Le\n", (long double) (EXP_FUNC(-ene/(kB*TEMP_K))/pf) );
  }

  return 0;
}
/* ****** */




