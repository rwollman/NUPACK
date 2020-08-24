/*
    mfe.c is part of the NUPACK software suite
    Copyright (c) 2007 Caltech. All rights reserved.
    Coded by: Robert Dirks, 6/2006 and Justin Bois 1/2007


    This function will calculate and print all mfe structures (if the
    -degenerate flag is selected) or one mfe structure (if not),
    taking into account symmetry corrections.  Consequently, if
    -degenerate is chosen, this could scale as poorly as exponential
    with regard to space and time.
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
  int seqNum[ MAXSEQLENGTH+1]; 
  int isNicked[ MAXSEQLENGTH];
  int nNicks = 0;


  int nicks[MAXSTRANDS];
  int nickIndex;
  int **etaN;
  int j, pf_ij;

  int complexity = 3;
  int length, tmpLength;
  DBL_TYPE mfe;
  int i;
  int vs;
  char inputFile[MAXLINE];
  char outFile[MAXLINE];
  int inputFileSpecified;
  FILE *fp;

  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};


  strcpy( inputFile, "");

  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);
  if(NupackShowHelp) {
    printf("Usage: mfe [OPTIONS] PREFIX\n");
    printf("Compute and store the minimum free energy and the MFE\n");
    printf("secondary structure(s) of the input sequence.\n");
    printf("Example: mfe -multi -T 25 -material dna example\n");
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
      !ReadInputFile( inputFile, seq, &vs, NULL, NULL, NULL) ) {
    if (inputFileSpecified == 0) getUserInput( seq, &vs, NULL, NULL);
    else abort();
  }
  strncpy(outFile,inputFile,strlen(inputFile)-3);
  outFile[strlen(inputFile)-3] = '\0';
  strcat(outFile,".mfe");


  header( argc, argv, "mfe",outFile);
  printInputs( argc, argv, seq, vs, NULL, NULL, outFile);


  if( !DO_PSEUDOKNOTS ) {
    complexity = 3;
  }
  else {
    complexity = 5;
  }

  tmpLength = strlen( seq);
  convertSeq(seq, seqNum, tmpLength);

  mfe = mfeFullWithSym( seqNum, tmpLength, &mfeStructs, complexity, DNARNACOUNT,
                        DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, vs,
                        ONLY_ONE_MFE, SODIUM_CONC, MAGNESIUM_CONC, 
			USE_LONG_HELIX_FOR_SALT_CORRECTION);


  //the rest is for printing purposes
  tmpLength = length = strlen( seq);


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

  /* Old way of printing information, cumbersome
  printf("Sequence/Structure");
  for( i = 19; i < tmpLength; i++) {
    printf(" ");
  }
  printf("  Diff  Sym Energy\n");
  printf("%s\n", seq);
  */


  // Put a new line before first MFE structure
  fp = fopen(outFile,"a");
  fprintf(fp,"\n");
  fclose(fp);

  PrintDnaStructures( &mfeStructs, etaN, nicks, vs,outFile);

  clearDnaStructures( &mfeStructs);

 
  for( i = 0; i <= length-1; i++) {
    for( j = i-1; j <= length-1; j++) {
      pf_ij = pf_index(i,j,length);
      free( etaN[pf_ij]);
    }
  }
  free( etaN);

  return 0;
}
/* ****** */
