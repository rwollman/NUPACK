/*
  distributions.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  This program does thermodynamic analysis of interacting nucleic acid
  strands in a box where the partition functions for the complexes are
  known, as described in Dirks, Bois, Schaeffer, Winfree, and Pierce,
  "Thermodynamic Analysis of interacting nucleic acid strands", SIAM
  Review, 2007.  Variable names in the code should be referenced
  with variable names in that paper.

  The program computes the partition function, $Q_{box}$, and
  population probabilities for a small box containing only a few
  strands.  The program to do this calculation is CalcDist.c.

  For usage instructions, input and output formats, etc., see the
  associated manual.
*/


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "InputFileReader.h"
#include "OutputWriter.h"
#include "CalcDist.h"
#include "constants.h"
#include "ReadCommandLine.h"

/* ******************************************************************************** */
/*                                    BEGIN MAIN                                    */
/* ******************************************************************************** */

int main(int argc, char *argv[]) {
  
  int i,j; // Counter
  time_t StartTime; // Time calculation started
  time_t EndTime; // Time calculation ended
  char *StartTimeStr; // The string representing the start time
  char cxFile[MAXLINE]; // File containing complex ID's and free energies
  char countFile[MAXLINE]; // File containing initiatial monomer counts
  char logFile[MAXLINE]; // File containing data about the calculation
  char distFile[MAXLINE];  // Name of file for distributions
  char LambdaFile[MAXLINE]; // Name of file to which Lambda is written
  int numSS; // Number of single-strand (monomer) types. 
  int numSS0; // Number of monomer types including those with zero concentration
  int numTotal; // Total number of complexes
  int nTotal; // Total number of permutations
  int newnTotal; // Total number of permutations with unacheivable ones cut out
  int nComments; // Number of comment lines in .cx file before data start
  int LargestCompID; // Largest complex ID
  int SortOutput; // Sorting options for output
  int WriteLambda; // = 1 for writing out the entries in \Lambda
  int quiet; // = 1 for no displays of results to screen
  int NoPermID; // = 1 if there are no perumtation IDs in 2nd column of input file
  int WriteLogFile; // = 1 if information is to be written to a log file
  double kT; // The thermal energy in kcal/mol.
  double MaxSizeLambda; // Maximum number of populations to consider
  int Toverride; // = 1 if the user has enforced a temperature in the command line
  int **A; // A[i][j] is the number of monomers of type i in complex j
  double *G; // Free energies of complexes
  int *m0; // Total concentrations of single-species
  int maxm0; // Maximal entry in m0
  double **Pmn; // Pmn[j][n] = P(m_j = n)
  double *mEq; // mEq[j] = <m_j>
  int *numPermsArray; // Number of permutations of each species
  int *CompIDArray; // The complex ID's
  int *PermIDArray; // Permutation ID's
  double M; // Number of solvent molecules in the box
  FILE *fplog; // The logFile, which contains information about the run.
  FILE *fpdist; // The .dist file, which contains the output of the file.
  FILE *fpLam; // The file to which Lambda is written
  int NUPACK_VALIDATE ; // if set to 1 print out to 14 decimal places

  /* version 3 output */
  int v3;
  
  // Read command line arguments
  ReadCommandLine(argc,argv,cxFile,countFile,logFile,distFile,LambdaFile,&SortOutput,
                  &WriteLambda,&MaxSizeLambda,&kT,&quiet,&WriteLogFile,&Toverride,
                  &NoPermID,&NUPACK_VALIDATE, &v3);

  
  // Get the start time of the calculation
  StartTime = time(NULL);
  StartTimeStr = ctime(&StartTime);

  // Print run information to log file (.log file) and output file (.dist file)
  // and also .Lam file if necessary
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"w")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    //fclose(fplog);
    //fplog=stderr;
    fprintf(fplog,"*NUPACK %s\n", NUPACK_VERSION);
    fprintf(fplog,"*This is %s, a log file generated for a distributions\n",logFile);
    fprintf(fplog," calculation using input files %s and %s.\n",cxFile,countFile);
    fprintf(fplog,"*Command used: ");
    for (i = 0; i < argc; i++) {
      fprintf(fplog,"%s ",argv[i]);
    }
    fprintf(fplog,"\n");
    fprintf(fplog,"*Time calculation was begun: %s",StartTimeStr);
    fprintf(fplog,"*Pertinent run data:\n");
    fprintf(fplog,"   --Input files: %s\n",cxFile);
    fprintf(fplog,"                  %s\n",countFile);
    fprintf(fplog,"   --Output files: %s\n",distFile);
    if (WriteLambda == 1) {
      fprintf(fplog,"                   %s\n",LambdaFile);
    }
    else {
      fprintf(fplog,"                   Lambda is not written to a file.\n");
    }
    fprintf(fplog,"   --Single-species counts are: \n");
    fclose(fplog);
  }
  
  // Write information to dist file
  if ((fpdist = fopen(distFile,"w")) == NULL) {
    if (quiet == 0) {
      printf("Error opening %s.\n\nExiting....\n",distFile);
    }
    exit(ERR_DIST);
  }
  fprintf(fpdist,"%% NUPACK %s\n", NUPACK_VERSION);
  fprintf(fpdist,"%% This is %s, an output file generated for a \n",distFile);
  fprintf(fpdist,"%% calculation of equilibrium complex count distributions.\n");
  fprintf(fpdist,"%% Time calculation was begun: %s",StartTimeStr);
  fprintf(fpdist,"%% Inital monomer counts:\n");
  fclose(fpdist);
  
  // Write information to .Lam file
  if (WriteLambda) {
    if ((fpLam = fopen(LambdaFile,"w")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",LambdaFile);
      }
      exit(ERR_LAMBDA);
    }
    fprintf(fpLam,"%% NUPACK %s\n", NUPACK_VERSION);
    fprintf(fpLam,"%% This is %s, an output file generated for a \n",LambdaFile);
    fprintf(fpLam,"%% calculation of equilibrium complex count distributions.\n");
    fprintf(fpLam,"%% The contents of this file is the Lambda, the set of\n");
    fprintf(fpLam,"%% population vectors consistent with conservation of mass.\n"); 
    fprintf(fpLam,"%% Time calculation was begun: %s",StartTimeStr);
    fclose(fpLam);
  }
  
  // Get the size of the system.
  getSize(&numSS,&numTotal,&nTotal,&LargestCompID,&numPermsArray,&nComments,
          cxFile,countFile,quiet);
  
  // Read input files and sort if necessary.
  // Note: A, G, and m0 are all allocated in ReadInput
  if (NoPermID == 0) {
    ReadInputFilesPerm(&A,&G,&CompIDArray,&PermIDArray,&m0,&M,&numSS,&numSS0,
                       &newnTotal,nTotal,cxFile,countFile,&kT,Toverride,logFile,
                       distFile,quiet,WriteLogFile);
  }
  else {
    ReadInputFiles(&A,&G,&CompIDArray,&PermIDArray,&m0,&M,&numSS,&numSS0,&numTotal,
                   numPermsArray,cxFile,countFile,&kT,Toverride,logFile,distFile,
                   quiet,WriteLogFile);
  }
  
  // First write out information about the run
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    //fplog=stderr;
    fprintf(fplog,"   --Number of single-stranded species considered: %d\n",numSS);
    if (NoPermID == 0) {
      fprintf(fplog,"   --Total number of permutations considered: %d\n",newnTotal);
    }
    else {
      fprintf(fplog,"   --Total number of complexes considered: %d\n",numTotal);
    }
    fprintf(fplog,"   --Temperature (in deg. C): %g\n",kT/kB - 273.15);
    fclose(fplog);
  }
  
  maxm0 = maxint(m0,numSS); // Maximal entry in m0
  
  if (NoPermID == 0) {
    // Allocate memory for outputs of CalcDist
    mEq = (double *) malloc(newnTotal * sizeof(double));
    Pmn = (double **) malloc(newnTotal * sizeof(double *));
    for (j = 0; j < newnTotal; j++) {
      Pmn[j] = (double *) malloc((maxm0+1) * sizeof(double));
    }
    
    // Do the calculation
    CalcDist(mEq,Pmn,A,G,m0,M,numSS,newnTotal,MaxSizeLambda,LambdaFile,kT,
             WriteLambda,CompIDArray,PermIDArray,quiet,logFile,WriteLogFile,
             NoPermID);
    
    // Write output
    WriteOutput(mEq,Pmn,G,A,CompIDArray,PermIDArray,LargestCompID,numSS,newnTotal,
                nTotal,nComments,maxm0,kT,cxFile,SortOutput,distFile,quiet,NoPermID,NUPACK_VALIDATE);
  }
  else {
    // Allocate memory for outputs of CalcDist
    mEq = (double *) malloc(numTotal * sizeof(double));
    Pmn = (double **) malloc(numTotal * sizeof(double *));
    for (j = 0; j < numTotal; j++) {
      Pmn[j] = (double *) malloc((maxm0+1) * sizeof(double));
    }
    
    // Do the calculation
    CalcDist(mEq,Pmn,A,G,m0,M,numSS,numTotal,MaxSizeLambda,LambdaFile,kT,
             WriteLambda,CompIDArray,PermIDArray,quiet,logFile,WriteLogFile,
             NoPermID);
    
    // Write output
    WriteOutput(mEq,Pmn,G,A,CompIDArray,PermIDArray,LargestCompID,numSS,numTotal,
                nTotal,nComments,maxm0,kT,cxFile,SortOutput,distFile,quiet,NoPermID,NUPACK_VALIDATE);
  }
  
  
  EndTime = time(NULL);
  if (quiet == 0) {
    printf("Elapsed time: %g seconds.\n",difftime(EndTime,StartTime));
  }
  
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    fprintf(fplog,"   --Elapsed time of calculation: %g seconds\n",
            difftime(EndTime,StartTime));
    fclose(fplog);
  }
  
  // Free memory
  for (i = 0; i < numSS; i++) {
    free(A[i]); // Allocated in ReadInput
  }
  free(A); // Allocated in ReadInput
  free(G); // Allocated in ReadInput
  free(numPermsArray); // Allocated in getSize
  free(CompIDArray); // Allocated in ReadInput
  free(PermIDArray); // Allocated in ReadInput
  if (NoPermID == 0) {
    for (j = 0; j < newnTotal; j++) {
      free(Pmn[j]); // Allocated in main
    }
  }
  else{
    for (j = 0; j < numTotal; j++) {
      free(Pmn[j]); // Allocated in main
    }
  }
  free(Pmn); // Allocated in main
  free(mEq); // Allocated in main
  free(m0);  // Allocated in ReadInput
  
  return 0; // Return
  
}
/* ******************************************************************************** */
/*                                  END MAIN SCRIPT                                 */
/* ******************************************************************************** */
