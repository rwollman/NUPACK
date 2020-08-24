/*
  count.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 8/2006 and Justin Bois 1/2007
  
  The purpose of this program is to count the number of legal secondary
  structures for a given sequence.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <thermo/core.h>

/* ************************************************ */

int main( int argc, char *argv[] ) {
  
  char seq[ MAXSEQLENGTH];
  int seqNum[ MAXSEQLENGTH];

  DBL_TYPE pf;
  int vs;

  int complexity = 5;
  
  char inputFile[ MAXLINE];
  int inputFileSpecified;
  
  strcpy( inputFile, "");
  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  if(NupackShowHelp) {
    printf("Usage: count [OPTIONS] PREFIX\n");
    printf("Count the number of possible secondary structures of the input sequence\n");
    printf("Example: count -multi pnas04_hcr_basic\n");
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  header( argc, argv, "count", "screen");
  
  if( !inputFileSpecified || 
      !ReadInputFile( inputFile, seq, &vs, NULL, NULL, NULL) ) {
       if (inputFileSpecified == 0) getUserInput( seq, &vs, NULL, NULL);
       else abort();
  }

  
  printInputs( argc, argv, seq, vs, NULL, NULL, "screen");
  
  if( !DO_PSEUDOKNOTS ) {
    complexity = 3;
  }
  else {
    complexity = 5;
  }
  
  int tmpLength=strlen(seq);
  convertSeq(seq, seqNum, tmpLength);
  pf = pfuncFull(seqNum, complexity, COUNT, DANGLETYPE, TEMP_K-ZERO_C_IN_KELVIN,
        0, 1.0, 0.0, 0);
  
  printf("%s\n%s Total number of secondary structures:\n",COMMENT_STRING,
        COMMENT_STRING);
  printf( "%14.14Le\n", (long double) pf); 
 
  return 0;
}
/* ****** */













