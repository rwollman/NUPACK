/*
  FracPair.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  For use with Concentrations.c.

  Computes the fraction of each strand that forms base pairs.
*/


#include "FracPair.h"
#include "constants.h"

/* ******************************************************************************** */
void FracPair(int numSS, int nTotal, int quiet, int NoPermID, int LargestCompID, 
              int *numPermsArray, char *eqFile, char *conFile, char *pairsFile, 
              char *fpairsFile, double cutoff,int NUPACK_VALIDATE) {
  /*
    Creates the matrix of base-pair fractions.

    This is basically a script, since it only deals with the output file from
    the program.  This could also be written in Perl, Python, or some scripting
    language, but since it's such a common calculation, we have included it in
    the functionality of concentrations itself.

    We don't bother messing with the units of x0 and x since they are the same
    and we always use their ratio.
  */

  // Go back through the input file and load in the permutation information.
  int i,j,k,n; // Counters
  char line[MAXLINE]; // A line from a file
  char InStr[MAXLINE]; // Input string from file
  char *tok; // Token
  char tokseps[] = " \t\n"; // Token separators  
  double *x; // Concentrations of complexes
  double PairProb; // Pair probability of entry
  double *x0; // The initial concentrations
  double **PairMat; // Base pairing matrix
  double PairFrac; // fraction of bases with pair type
  int *StrandLength;
  int N; // Number of bases total
  int *StrandKey; // Key for what index is in what strand
  int CompID, PermID; // Complex and permutation ID
  int basei, basej; // indices in pairing matrix
  int **ComplexKey; // ComplexKey[CompID][PermID] = corresponding index in x
  int ReachedEnd; // = 1 if we reached the end of the file
  FILE *fp; // file handle for con file and pairs file
  FILE *fpeq; // file handle for eq file
  FILE *fpfpairs; // file handle for fpairs file

  // Allocate memory for strand lengths
  StrandLength = (int *) malloc(numSS * sizeof(int));
  // Initialize
  for (i = 0; i < numSS; i++) {
    StrandLength[i] = 0;
  }

  // Allocate memory for concentrations
  x0 = (double *) malloc(numSS * sizeof(double));
  x = (double *) malloc(nTotal * sizeof(double));

  // Allocate memory for Complexes key
  ComplexKey = (int **) malloc(LargestCompID * sizeof(int *));
  for (j = 0; j < LargestCompID; j++) {
    ComplexKey[j] = (int *) malloc(numPermsArray[j] * sizeof(int));
  }

  /* ************ Read in intial concentrations ********************* */
  // Open the con file
  if ((fp = fopen(conFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",conFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_CON);
  }

  // Blow through comments and blank lines
  while (fgets(line,MAXLINE,fp) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n'));

  // Read in the initial concentrations
  if ( (tok = strtok(line,tokseps)) != NULL) {
    x0[0] = str2double(tok);
  }
  for (i = 1; i < numSS; i++) {
    fgets(line,MAXLINE,fp);
    if ( (tok = strtok(line,tokseps)) != NULL) {
      x0[i] = str2double(tok);
    }
  } 

  fclose(fp);
  /* *************************************************************** */


  /* *************** Get data from eq file ************************* */
  // Open the eq file
  if ((fpeq = fopen(eqFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",eqFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_EQ);
  }

 
  // Blow through comments and blank lines and pull out sequence lengths
  while (fgets(line,MAXLINE,fpeq) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n')) {

    if (strncmp(line,"% id sequence",10) == 0) { // This is the line before sequences
      for (i = 0; i < numSS; i++) {
        fgets(line,MAXLINE,fpeq);
        tok = strtok(line,tokseps); // tok = '%'
        tok = strtok(NULL,tokseps); // tok = sequence id
        tok = strtok(NULL,tokseps); // tok = sequence
        StrandLength[i] = strlen(tok);
      }
    }
  }

  // Make sure we got the sequence lengths
  for (i = 0; i < numSS; i++) {
    if (StrandLength[i] == 0) {
      if (quiet == 0) {
        printf("Error in getting sequence lengths from eq file.\n");
        printf("\nExiting....\n\n");
        exit(ERR_NOSEQEQ);
      }
    }
  }

  // Get N (number of bases)
  N = sumint(StrandLength,numSS);

  // Allocate memory for base pairing matrix and initialize
  PairMat = (double **) malloc(N * sizeof(double *));
  for (basei = 0; basei < N; basei++) {
    PairMat[basei] = (double *) malloc((N+1) * sizeof(double));
    for (basej = 0; basej < N+1; basej++) {
      PairMat[basei][basej] = 0.0;
    }
  }

  // Make the strand key
  StrandKey = (int *) malloc(N * sizeof(int));
  k = 0;
  for (i = 0; i < numSS; i++) {
    for (j = 0; j < StrandLength[i]; j++) {
      StrandKey[k] = i+1;
      k++;
    }
  }

  // Now go through rest of eq file and get concentration information
  PermID = 0; // Initialize to zero in case NoPermID = 1.
  k = 0;
  while (k < nTotal) {
    tok = strtok(line,tokseps);   // Complex ID
    CompID = atoi(tok) - 1;
    if (NoPermID == 0) {
      tok = strtok(NULL,tokseps);   // Permutation ID
      PermID = atoi(tok) - 1;
    }
    ComplexKey[CompID][PermID] = k;
    for (i = 0; i < numSS; i++) {  // stoichiometry 
      tok = strtok(NULL,tokseps);
    }
    tok = strtok(NULL,tokseps); // free energy
    tok = strtok(NULL,tokseps); // concentration
    x[k] = str2double(tok); 
    fgets(line,MAXLINE,fpeq);
    k++;
  }
  // close the eq file
  fclose(fpeq);
  /* *************************************************************** */


  /* *************** Get data from pairs file ********************** */
  if ((fp = fopen(pairsFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",pairsFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_PAIRSFILE);
  }


  // Blow through comments to first record
  while (fgets(line,MAXLINE,fp) != NULL && strncmp(line,"% complex",9) != 0
      && strncmp(line,"% composition",13) != 0);

  
  if (strncmp(line,"% complex",9) == 0 || strncmp(line,"% composition",13) == 0) { 
    // We have data in the file
    ReachedEnd = 0;
  }
  else { // We're already at the end of the file; no data
    ReachedEnd = 1;
  }

  // Go through records to get pairing data
  PermID = 0; // Initialize to zero in case NoPermID = 1.
  while (ReachedEnd == 0) {
    // Get the complex ID number
    i = 0;
    if (strncmp(line, "% complex", 9) == 0) {
      i = 9;  // Index in line string in v3.0.x
    }
    else if (strncmp(line, "% composition",13) == 0) {
      i = 13; // Index in line string in v3.x x > 0
    }
    n = i;
    while (isdigit(line[i])) {
      InStr[i-n] = line[i];
      i++;
    }
    InStr[i-n] = '\0';  // Add null character
    CompID = atoi(InStr) - 1;
    
    // Get the permutation ID number
    if (NoPermID == 0) {
      if (strncmp(line, "% complex", 9) == 0) {
        i += 6; // This skips us past the "-order" to the perm ID
      }
      else if (strncmp(line, "% composition",13) == 0) {
        i += 9; // This skips us past the "-ordering" to the perm ID
      }
      n = i;
      while (isdigit(line[i])) {
        InStr[i-n] = line[i];
        i++;
      }
      InStr[i-n] = '\0';  // Add null character
      PermID = atoi(InStr) - 1;
    }
    
    k = ComplexKey[CompID][PermID]; // index in x corresponding to this record

    // Number of bases total (already have it)
    fgets(line,MAXLINE,fp);

    
    fgets(line,MAXLINE,fp); // This is the first base-pair and pair prob
    // Go through pair potentials 
    while (line[0] != '%' && line[0] != '\n' && line[0] != '\0') {
      tok = strtok(line,tokseps);
      basei = atoi(tok) - 1; // Index i
      basej = atoi(strtok(NULL,tokseps)) - 1; // Index j
      PairProb = str2double(strtok(NULL,tokseps)); // pair prob
      PairMat[basei][basej] += PairProb*x[k];
      fgets(line,MAXLINE,fp);
    }

    // Scan ahead to the next record
    while (fgets(line,MAXLINE,fp) != NULL && strncmp(line,"% complex",9) != 0
        && strncmp(line,"% composition",13) != 0);

    if (strncmp(line,"% complex",9) != 0 && strncmp(line,"% composition",13) != 0) {
      // Hit end of file
      ReachedEnd = 1;
    }
  }
  fclose(fp);
  /* *************************************************************** */

  
  // Write out results
  if ((fpfpairs = fopen(fpairsFile,"a")) == NULL) {
    if (quiet == 0) {
      printf("Error opening %s.\n\nExiting....\n",fpairsFile);
    }
    exit(ERR_FPAIRS);
  }

  // Print total number of bases
  fprintf(fpfpairs,"%d\n",N);

  for (basei = 0; basei < N; basei++) {
    for (basej = 0; basej < N+1; basej++) {
      if (basej < N) {
        if (x0[StrandKey[basei]-1] > 0.0 && x0[StrandKey[basej]-1] > 0.0) {
          if(basei < basej) {
            PairFrac = PairMat[basei][basej]/x0[StrandKey[basei]-1];
            if (PairFrac >= cutoff) {
              fprintf(fpfpairs,"%d\t%d\t",basei + 1 , basej + 1);
              if(!NUPACK_VALIDATE) {
                fprintf(fpfpairs,"%8.6e\n",PairFrac);
              } else {
                fprintf(fpfpairs,"%.14e\n",PairFrac);
              }
            }
          } else {
            PairFrac = PairMat[basej][basei]/x0[StrandKey[basei]-1]; 
            if (PairFrac >= cutoff) { 
              fprintf(fpfpairs,"%d\t%d\t",basei + 1 , basej + 1); 
              if(!NUPACK_VALIDATE) { 
                fprintf(fpfpairs,"%8.6e\n",PairFrac); 
              } else { 
                fprintf(fpfpairs,"%.14e\n",PairFrac); 
              } 
            } 
          }
        }
      } else { // basej == N
        if (x0[StrandKey[basei]-1] > 0.0) {
          PairFrac = PairMat[basei][basej] / x0[StrandKey[basei]-1];
          if (PairFrac >= cutoff) {
            fprintf(fpfpairs,"%d\t%d\t",basei + 1 , basej + 1);
            if(!NUPACK_VALIDATE) {
              fprintf(fpfpairs,"%8.6e\n",PairFrac);
            } else {
              fprintf(fpfpairs,"%.14e\n",PairFrac);
            }
          }
        }
      }
    }
  }

  fclose(fpfpairs);

  free(StrandLength);
  free(StrandKey);
  free(x0);
  free(x);
  for (i = 0; i < N; i++) {
    free(PairMat[i]);
  }
  free(PairMat);
  for (i = 0; i < LargestCompID; i++) {
    free(ComplexKey[i]);
  }
  free(ComplexKey);

}
/* ******************************************************************************** */


