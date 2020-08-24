/*
  subopt.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 6/2006 and Justin Bois 1/2007
  
  This is similar to mfe, but allows for the additional input of an
  energy range.  All structures whose energy is within the range of
  the mfe will be output.  Note that this algorithm slows down
  exponentially as the range increases.
  
  The complexity setting mimics that of pfunc.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <thermo/core.h>

/* ************************************************ */

int main( int argc, char *argv[] ) {
  
  //This function will calculate the all structures within a fixed
  //range of the "algorithmic" mfe.  The algorithmic mfe is the
  //minimum free energy of a structure, ignoring symmetry corrections.
  
  //The output is a list of the structures, including the difference
  //from the algorithmic mfe, followed by the symmetry factor and the
  //corrected energies, after adjusting for symmetries.  The list is
  //sorted by the corrected energies.
  
  char seq[ MAXSEQLENGTH];
  int seqNum[ MAXSEQLENGTH+1];
  int isNicked[ MAXSEQLENGTH];
  int nNicks = 0;
  
  int nicks[MAXSTRANDS];
  int nickIndex;
  int **etaN;	
  int complexity = 3;
  int length, tmpLength;
  float gap = -1;
  int i;
  int vs;
  char outFile[MAXLINE];
  int inputFileSpecified;
  FILE *fp;
  
  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY}; 
  //this struct will store
  //all the structures within the given range
  
  char inputFile[ MAXLINE];
  strcpy( inputFile, "");
  
  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  
  if(NupackShowHelp) {
    printf("Usage: subopt [OPTIONS] PREFIX\n");
    printf("Calculate and store all structures within the specified energy gap\n");
    printf("of the MFE structure.\n");
    printf("Example: subopt -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    exit(1);
  }

  if( !inputFileSpecified ) {
    printf("Enter output file prefix: ");
    scanf("%s", inputFile);
    strcat(inputFile,".in"); // Here, .in is just a placeholder
  }
  
  if( !inputFileSpecified ||
     !ReadInputFile( inputFile, seq, &vs, &gap, NULL, NULL) ) {
       if (inputFileSpecified==0) getUserInput( seq, &vs, &gap, NULL);
       else abort();
     }
  strncpy(outFile,inputFile,strlen(inputFile)-3);
  outFile[strlen(inputFile)-3] = '\0';
  strcat(outFile,".subopt");
  
  header( argc, argv, "subopt", outFile);
  printInputs( argc, argv, seq, vs, &gap, NULL, outFile);
  
  // Add newline for stylistic reasons
  fp = fopen(outFile,"a");
  fprintf(fp,"\n");
  fclose(fp);
  
  if( !DO_PSEUDOKNOTS ) {
    complexity = 3;
  }
  else {
    complexity = 5;
  }
  
  tmpLength = length = strlen( seq);
  convertSeq(seq, seqNum, tmpLength);

  mfeFullWithSym_SubOpt( seqNum, tmpLength, &mfeStructs, complexity, 
                        DNARNACOUNT, DANGLETYPE, 
                        TEMP_K - ZERO_C_IN_KELVIN,
			 vs, (DBL_TYPE) gap, 0, SODIUM_CONC, MAGNESIUM_CONC,
			 USE_LONG_HELIX_FOR_SALT_CORRECTION);

  //the rest is for printing purposes
  
  for( i = 0; i < tmpLength; i++) {
    isNicked[i] = 0;
    if( seq[i] == '+') {
      length--;
      isNicked[ i - nNicks++ -1] = 1;
    } 
  }
  
  //initialize nicks
  for( i = 0; i < MAXSTRANDS; i++) {
    nicks[i] = -1;
  }
  
  nickIndex = 0;
  for( i = 0; i < length; i++) {
    if( isNicked[i])
      nicks[ nickIndex++] = i;
  }
  
  //overkill, but convenient
  etaN = (int**) malloc( (length*(length+1)/2 + (length+1))*sizeof( int*));
  InitEtaN( etaN, nicks, length);
  
  PrintDnaStructures( &mfeStructs, etaN, nicks, vs, outFile);
  
  clearDnaStructures( &mfeStructs);
  
  return 0;
}
/* ****** */
