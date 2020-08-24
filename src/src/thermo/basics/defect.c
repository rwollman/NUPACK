/*
  defect.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 3/2009

  This program computes the difference between the equilibrium
  ensemble of structures and an inputted structure, i.e., n(s*).  For
  use with NUPACK.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <thermo/core.h>

/* ************************************************ */
extern DBL_TYPE *pairPrPbg;  //for pseudoknots
extern DBL_TYPE *pairPrPb;  //for pseudoknots

extern double CUTOFF;
extern int Multistranded;
extern int perm[MAXSTRANDS];
extern int seqlengthArray[MAXSTRANDS];
extern int nUniqueSequences; // Number of unique sequences entered

int main( int argc, char *argv[] ) {

  int i,j; // indices
  char seqChar[ MAXSEQLENGTH], parens[ MAXSEQLENGTH];
  int thepairs[MAXSEQLENGTH];
  int seqNum[ MAXSEQLENGTH+1];
  DBL_TYPE nsStar;
  DBL_TYPE mfe; // Minimum free energy (not really used)
  DBL_TYPE ene; // Free energy (used to check if the structure is legal)
  int vs;
  int complexity;
  int tmpLength;
  int seqlength;
  char inputFile[ MAXLINE];
  int nNicks;  // Number of strands
  int doCalc; // Whether we need to compute the pair probability matrix/mfe or not
  char precompFile[MAXLINE];  //File that has ppairs data
  char line[MAXLINE]; // A line from a file
  char tokseps[] = " \t\n"; // Token separators
  dnaStructures mfeStructs = {NULL, 0, 0, 0, NAD_INFINITY};
  FILE *fp;
  int inputFileSpecified;

  // -degenerate flag not used for defect calcs, force ONLY_ONE_MFE
  ONLY_ONE_MFE = 1;

  strcpy( inputFile, "");

  doCalc = 0;

  inputFileSpecified = ReadCommandLineNPK(argc, argv, inputFile);
  if(NupackShowHelp) {
    printf("Usage: defect [OPTIONS] PREFIX\n");
    printf("Calculate the ensemble defect of the input sequence with respect\n");
    printf("to the input structure\n");
    printf("Example: defect -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    printf("Additional:\n");
    printf(" -mfe                         compute the MFE defect\n");
    printf("\n");
    exit(1);
  }

  header( argc, argv, "defect", "screen");

  // Read input
  if( !inputFileSpecified ||
      !ReadInputFile( inputFile, seqChar, &vs, NULL, parens, thepairs) ) {
    if(inputFileSpecified == 0) {
      getUserInput( seqChar, &vs, NULL, parens);
      doCalc = 1;
    } else {
      abort();
    }
  }

  // Get the sequence length and convert the sequence
  seqlength = tmpLength = strlen( seqChar);
  convertSeq(seqChar, seqNum, tmpLength);

  // Get the number of strand breaks
  nNicks = 0;
  for (i = 0; i < tmpLength; i++) {
    if (seqChar[i] == '+') {
      nNicks++;
    }
  }

  // New sequence length
  seqlength -= nNicks;


  // Make sure the target structure is in thepairs format
  if (parens[0] != '\0') {
    getStructureFromParens( parens, thepairs, seqlength);
  }


  // ******************************
  // THIS DOESN'T QUITE WORK YET.  WE NEED A LEGAL STRUCTURE CHECKER!
  // Compute the free energy to check if it's a legal structure
  ene = naEnergyPairsOrParensFullWithSym( thepairs, NULL, seqNum, DNARNACOUNT, DANGLETYPE,
            TEMP_K - ZERO_C_IN_KELVIN, vs,
            SODIUM_CONC, MAGNESIUM_CONC,
            USE_LONG_HELIX_FOR_SALT_CORRECTION);


  // Check to see if the results is close to NAD_INFINITY and report
  // error if it is
  // if (ABS_FUNC(1.0 - ene/NAD_INFINITY) < INF_CUTOFF) {
  //   printf("\n\n*** Error: target structure has invalid base pair(s) or disconnected complex. Check your inputs. ***\n\n");
  //   return 0;
  // }
  // ******************************


  // Allocate memory for storing pair probabilities
  pairPr = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), sizeof(DBL_TYPE));


  // Check to see if we can use the input file to get the pair probs
  // Get the filename prefix
  if (strlen(inputFile)>3)
    strncpy(precompFile, inputFile, strlen(inputFile)-3);
  if (USE_MFE) {
    strcat(precompFile,".mfe");

    // Open the mfe file
    if ((fp = fopen(precompFile,"r")) == NULL) { // mfe files does not exist
      doCalc = 1;
    }
    else {      // Parse the file for data

      // Initialize the structure to be all unpaired
      for (i = 0; i < seqlength; i++) {
        pairPr[i*(seqlength+1) + seqlength] = 1.0;
      }

      // Blow through comments and blank lines
      while (fgets(line,MAXLINE,fp) != NULL && isdigit(line[0]) == 0);

      // Read in the total number of bases
      seqlength = atoi(strtok(line,tokseps));

      fgets(line,MAXLINE,fp); // Free energy of MFE

      fgets(line,MAXLINE,fp); // Possibly secondary structure in dot-parens
      if (!isdigit(line[0])) { // It is secondary structure in dot-parens
        fgets(line,MAXLINE,fp);
      }

      // The current line is the first base pair in the MFE
      if (isdigit(line[0])) { // If the structure is not empty
        i = atoi(strtok(line,tokseps)) - 1;
        j = atoi(strtok(NULL,tokseps)) - 1;
        pairPr[i*(seqlength+1) + j] = 1.0;
        pairPr[j*(seqlength+1) + i] = 1.0;
        pairPr[i*(seqlength+1) + seqlength] = 0.0;
      }

      // Take the rest data out of the file
      while(fgets(line,MAXLINE,fp) != NULL && isdigit(line[0])) {
        i = atoi(strtok(line,tokseps)) - 1;
        j = atoi(strtok(NULL,tokseps)) - 1;
        pairPr[i*(seqlength+1) + j] = 1.0;
        pairPr[j*(seqlength+1) + i] = 1.0;
        pairPr[i*(seqlength+1) + seqlength] = 0.0;
      }

      fclose(fp);

      // Compute nsStar
      nsStar = seqlength - expectedCorrectBases(thepairs,seqlength);
    }
  }
  else {
    strcat(precompFile,".ppairs");

    // Open the ppairs file
    if ((fp = fopen(precompFile,"r")) == NULL) { // ppairs files does not exist
      doCalc = 1;
    }
    else {      // Parse the file for data

      // Blow through comments and blank lines
      while (fgets(line,MAXLINE,fp) != NULL && isdigit(line[0]) == 0);

      // Read in the total number of bases
      seqlength = atoi(strtok(line,tokseps));

      // Take the data out of the file
      while(fgets(line,MAXLINE,fp) != NULL && isdigit(line[0])) {
        i = atoi(strtok(line,tokseps)) - 1;
        j = atoi(strtok(NULL,tokseps)) - 1;
        pairPr[i*(seqlength+1) + j] = (DBL_TYPE) str2double(strtok(NULL,tokseps));
        pairPr[j*(seqlength+1) + i] = pairPr[i*(seqlength+1) + j];
      }

      fclose(fp);

      // Compute nsStar
      nsStar = seqlength - expectedCorrectBases(thepairs,seqlength);
    }

  }


  // Perform the calculation if need be
  if (doCalc) {
    printInputs( argc, argv, seqChar, vs, NULL, parens, "screen");

    if (USE_MFE) {
      if( !DO_PSEUDOKNOTS ) {
        complexity = 3;
      }
      else {
        complexity = 5;
      }

      // Compute MFE and MFE structure
      mfe = mfeFullWithSym( seqNum, tmpLength, &mfeStructs, complexity, DNARNACOUNT,
          DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, vs,
          ONLY_ONE_MFE, SODIUM_CONC, MAGNESIUM_CONC,
          USE_LONG_HELIX_FOR_SALT_CORRECTION);

      // Compute nsStar from output
      nsStar = 0.0;
      for (i = 0; i < seqlength; i++) {
        if (thepairs[i] != (mfeStructs.validStructs)[0].theStruct[i]) {
          nsStar += 1.0;
        }
      }
    }
    else{
      // Allocate memory for storing pair probabilities
      pairPrPbg = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), sizeof(DBL_TYPE));
      pairPrPb = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), sizeof(DBL_TYPE));

      if( !DO_PSEUDOKNOTS ) {
        complexity = 3;
      }
      else {
        complexity = 5;
      }

      nsStar = nsStarPairsOrParensFull(seqlength, seqNum, thepairs, NULL,
               complexity, DNARNACOUNT,DANGLETYPE,
               TEMP_K - ZERO_C_IN_KELVIN, SODIUM_CONC,
               MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);


      free( pairPrPbg);
      free( pairPrPb);
    }
  }
  else { // Tell the user we used an existing file for the probs
    if (USE_MFE) {
      printf("%s Results computed using the MFE in the file %s\n",
       COMMENT_STRING,precompFile);
    }
    else {
      printf("%s Results computed using the pair probabilities in the file %s\n",
       COMMENT_STRING,precompFile);
    }
    printf("%s You should check to ensure the correct parameters were used to generate %s\n",
     COMMENT_STRING,precompFile);
  }


  if (USE_MFE) {
    printf("%s\n%s Fraction of correct nucleotides vs. MFE:\n", COMMENT_STRING,COMMENT_STRING);
  }
  else {
    printf("%s\n%s Ensemble defect n(s,phi) and normalized ensemble defect n(s,phi)/N:\n", COMMENT_STRING,COMMENT_STRING);
  }
  if(!NUPACK_VALIDATE) {
    printf("%4.3Le\n", (long double) nsStar);
    printf("%4.3Le\n", (long double) nsStar / seqlength);
  } else {
    printf("%16.14Le\n",(long double) nsStar);
    printf("%16.14Le\n", (long double) nsStar / seqlength);
  }


  free( pairPr);

  return 0;
}
/* ****** */


