/*
  CalcDist.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  For use with Distributions.c.

  This program solves the problem presented in Dirks, et al.,
  "Thermodynamic analysis of interacting nucleic acid strands", in
  progress.  All variable names are chosen to match those in that
  paper.

  The set $\Lambda$ is the set of all possible populations in a box
  containing a solution of complexes such that mass is conserved.
  This program enumerates all entries in $\Lambda$ and calculates
  their probabilty of occuring at equilibrium.  From this calculation,
  the partition function for the box, the probability distributions
  for the counts for each complex species, and the expectation value
  for the counts of each complex species.

  The inputs are as follows:
  A: A 2-D array; A[i][j] is the number of monomers of type i in complex j
  G: An array containing the corresponding complex free energies in units of kT.
     G[j] corresponds to the entries A[..][j].
  m0: Initial counts of the unit-size complexes.
  M: The number of solvent molecules in the box.
  numSS: The number of single-species (monomers) in the system.
  numTotal: The total number of complexes.
  MaxSizeLambda: The maximum number of elements in the set $\Lambda$.
  eqFile: The file to which the expectation values for the counts of the complexes
          is written.
  probFile: The file to which the probability distributions for each complex
            type is written.
  LambdaFile: The file to which the elements of the set Lambda and their
              respective probabilities are written.
  kT: kT in kcal/mol.
  quiet: = 1 for no printing of messages (except error messages) to screen.  
  fp: file for printing information about the run.

  For formats of the input and output files, etc., see the associated README file.
*/

#include "CalcDist.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // Must compile with -lm option

#include "constants.h"


/* ******************************************************************************** */
/*                                BEGIN CALCDIST                                    */
/* ******************************************************************************** */
void CalcDist(double *mEq, double **Pmn, int **A, double *G, int *m0, double M, 
          int numSS, int numTotal, double MaxSizeLambda, char *LambdaFile, 
          double kT, int WriteLambda, int *CompIDArray, int *PermIDArray,
          int quiet, char *logFile, int WriteLogFile, int NoPermID) {

  int i,j,k,n; // Counters, i: single-species, j: complexes, k: comb. of complexes
  int numComplex; // Number of complexes that are not single-species
  int LastInc; // Last entry in mComplex that was incremented
  int *mMax; // the maximal allowed count of each species
  int *mComplex; // Temporary array storing the counts of complexes excl. single-spec.
  int **Lambda; // Populations that satisfy conservation of mass
  double *logP; // Log of the probabilities for each element in Lambda
  double *P; // Probability for each element in Lambda
  long double Qbox; // The partition function (Q_box) 
  long double addQbox; // One term in the sum of Q_box 
  double FreeEnergy; // The free energy of the system.
  int SizeLambda; // The length of lambda
  double logM; // log of total numbers of particles
  double RefSum; // Term added to partition function to define reference state
  int maxm0; // Maximal entry in m0
  FILE *fpLam; // The file to which the elements of Lambda are written
  FILE *fplog; // file handle for log file

  // Preliminaries: get constants
  numComplex = numTotal - numSS;
  logM = log(M);
  RefSum = 0;
  for (i = 0; i < numSS; i++) {
    RefSum += m0[i]*(G[i] - logM)+ log(factorial(m0[i]));
  }
  maxm0 = maxint(m0,numSS);

  /* **************** Allocate memory for arrays ************* */
  Lambda = (int **) malloc(((int)(MaxSizeLambda)) * sizeof(int *));
  mMax = (int *) malloc(numComplex * sizeof(int));
  mComplex = (int *) malloc(numComplex * sizeof(int));
  /* ********************************************************* */


  // Calculate maximal counts based on exhausting limiting monomer
  for (j = 0; j < numComplex; j++) {
    mMax[j] = 0;
    for (i = 0; i < numSS; i++) {
      if (A[i][numSS+j] > 0) {
        mMax[j] = max2(floor(m0[i]/A[i][numSS+j]),mMax[j]);
      }
    }
  }

  // Initialialize mComplex
  for (j = 0; j < numComplex; j++) {
    mComplex[j] = 0;
  }

  // Get the list of populations  
  // First write the trivial population (all single-stranded)
  UpdateLambda(&Lambda,0,mComplex,m0,A,numSS,numTotal);

  LastInc = numComplex - 1;
  SizeLambda = 1;
  while (LastInc >= 0) {
    LastInc = next(mComplex,m0,A,numSS,numComplex,numTotal,mMax,LastInc);
    if (LastInc == numComplex-1) {
      UpdateLambda(&Lambda,SizeLambda,mComplex,m0,A,numSS,numTotal);
      SizeLambda++;
    }

    if (SizeLambda >= MaxSizeLambda) {
      if (quiet == 0) {
        printf("Exceeded maximum number of %g combinations!\n",MaxSizeLambda);
        printf("Try increasing the maximal size of Lambda or\n");
        printf("try using NUPACK's concentrations program to analyze\n");
        printf("larger systems.\n\nExiting...\n");
      }
      exit(ERR_LAMBDATOOBIG);
    }
  }

  // Write out  the number of populations in Lambda
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    fprintf(fplog,"   --Number of populations in Lambda: %d\n",SizeLambda);
    fclose(fplog);
  }

  if (quiet == 0) {
    printf("There are %d populations in Lambda.\n",SizeLambda);
  }

  // Allocate memory for probabilities
  logP = (double *) malloc(SizeLambda * sizeof(double));
  P = (double *) malloc(SizeLambda * sizeof(double));

  // Calculate terms in the partition function
  Qbox = 0.0;
  for (k = 0; k < SizeLambda; k++) {
    logP[k] = 0;
    
    j = 1;
    while (j < 2*Lambda[k][0]) {
      logP[k] += Lambda[k][j+1]*(logM - G[Lambda[k][j]]) 
      - log(factorial(Lambda[k][j+1]));
      j = j + 2;
    }

    logP[k] += RefSum;
    if ((addQbox = expl((long double) logP[k])) == HUGE_VAL) {
      if (quiet == 0) {
        printf("Overflow error in calculation of Q_{box}!\n");
        printf("Free energies of complexes is too high or the box is too big.\n");
        printf("Either adjust these or run the calculation on a large.\n");
        printf("system using the program concentrations.\n");
        printf("Exiting....\n\n");
      }
      exit(ERR_QBOXTOOBIG);
    }
    else {
      Qbox += addQbox;
    }
  }

  FreeEnergy = -log(Qbox);

  // Compute the probabilities of each element in Lambda
  for (k = 0; k < SizeLambda; k++) {
    logP[k] += FreeEnergy;
    P[k] = exp(logP[k]);
  }


  // Write out Lambda if necessary
  if (WriteLambda) {
    if ((fpLam = fopen(LambdaFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error in opening %s!\n\n",LambdaFile);
        printf("Exiting...\n");
      }
      exit(ERR_LAMBDA);
    }
    else {
      for (k = 0; k < SizeLambda; k++) {
        // The probability of the population occuring
        fprintf(fpLam,"%8.6e\t",P[k]);
        j = 1;
        if (NoPermID == 1) {
          while (j < 2*Lambda[k][0]) {
            fprintf(fpLam,"%d\t",CompIDArray[Lambda[k][j++]]);
            fprintf(fpLam,"%d\t",Lambda[k][j++]);
          }
        }
        else {
          while (j < 2*Lambda[k][0]) {
            fprintf(fpLam,"%d\t",CompIDArray[Lambda[k][j]]);
            fprintf(fpLam,"%d\t",PermIDArray[Lambda[k][j++]]);
            fprintf(fpLam,"%d\t",Lambda[k][j++]);
          }
        }
        fprintf(fpLam,"\n");
      }
      fclose(fpLam);
    }
  }


  // Initialize Pmn
  for (j = 0; j < numTotal; j++) {
    for (n = 0; n <= maxm0; n++) {
      Pmn[j][n] = 0.0;
    }
  }

  // Put the elements in Pmn
  for (k = 0; k < SizeLambda; k++) {
    j = 1;
    while (j < 2*Lambda[k][0]) {
      Pmn[Lambda[k][j]][Lambda[k][j+1]] += P[k];
      j = j + 2;
    }
  }

  // The probability that the count is zero
  for (j = 0; j < numTotal; j++) {
    Pmn[j][0] = 1.0 - sum(Pmn[j],maxm0+1);
    // Fix in case of precision error
    if (Pmn[j][0] < 0.0) {
      if (fabs(Pmn[j][0]) > NUM_PRECISION) {
	printf("Error: negative probability encountered.\n\nExiting....\n");
	exit(ERR_NEGATIVEPROB);
      }
      else {
	Pmn[j][0] = 0.0;
      }
    }
  }

  // Calculate the equilibrium counts
  for (j = 0; j < numTotal; j++) {
    mEq[j] = 0;
    for (n = 0; n <= maxm0; n++) {
      mEq[j] += n*Pmn[j][n];
    }
  }

  // Write out free energy
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    fprintf(fplog,"   --Free energy of the box: %8.6e kT, or %8.6e kcal\n",
            FreeEnergy,FreeEnergy*kT/AVOGADRO);
    fclose(fplog);
  }
  if (quiet == 0) {
    printf("Free energy of the box: %8.6e kT, or %8.6e kcal\n",
      FreeEnergy,FreeEnergy*kT/AVOGADRO);
  }


  // Free memory
  free(mMax);
  for (k = 0; k < SizeLambda; k++) {
    free(Lambda[k]);
  }
  free(Lambda);
  free(logP);
  free(P);
  free(mComplex);

}
/* ******************************************************************************** */
/*                                  END CALCDIST                                    */
/* ******************************************************************************** */


/* ******************************************************************************** */
int next(int *mComplex, int *m0, int **A, int numSS, int numComplex, int numTotal,
	 int *mMax, int LastInc) {

  /*
    Given the previous entries in Lambda, computes the next one.

    An example of what Lambda looks like for a system with two
    single-species and a maximal complex size of 3 is below:
    A     B         AA     AB    BB
    3     3         0      0     0
    3     1         0      0     1
    2     2         0      1     0
    2     0         0      1     1
    1     1         0      2     0
    0     0         0      3     0
    1     3         1      0     0
    1     2         1      0     1
    0     2         1      1     0
    0     0         1      1     1

    The single stranded counts are known from the multistranded by
    conservation of mass.  This function looks at the previous entry
    in Lambda and computes the next one.  For example, assume the
    previous entry is the top on in the table above.  The count of BB
    is advanced one and the resulting counts of the single species are
    checked.  There are both nonnegative, so the entry (0 0 1) for the
    complexes is accepted.  The last complex to be incremented is
    stored and returned as LastInc.  On the next call, BB is
    incremented again.  This time, the single strand count on strand B
    is -1, which is illegal, so the entry (0 0 2) is rejected.
    Therefore, BB is set to zero, and LastInc is incremented to
    correspond to AB.  The rejected entry is not inclulded in Lambda,
    and the next call to this function increments AB.  A check is
    again performed, and so on....
  */

  int i; // Counter

  (mComplex[numComplex-1])++;

  if (mComplex[numComplex-1] <= mMax[numComplex-1] 
      && NegCheck(mComplex,m0,A,numSS,numTotal) ) {
    LastInc = numComplex-1;
  }
  else {
    if (LastInc > 0) {
      (mComplex[LastInc-1])++;
      for (i = LastInc; i < numComplex-1; i++) {
	mComplex[i] = 0;
      }
      mComplex[numComplex-1] = -1;
      LastInc--;
    }
    else { // LastInc == 0
      LastInc = -1;
    }
  }

  return LastInc;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void UpdateLambda(int ***Lambda, int LambdaIndex, int *mComplex, int *m0, int **A, 
		  int numSS, int numTotal) {

  /*
    Adds the entry strored in mComplex to Lambda.
  */

  int i,j,k; // Counters
  int n; // number of nonzero entries in this row of lambda
  int *mss;
  int dotprod; // Dot product

  mss = (int *) malloc(numSS * sizeof(int));

  // Get the counts for a single-strand
  for (i = 0; i < numSS; i++) {
    dotprod = 0;
    for (j = numSS; j < numTotal; j++) {
      dotprod += A[i][j]*mComplex[j-numSS];
    }
    mss[i] = m0[i] - dotprod;
  }

  // Find out how many nonzero entries we have
  n = nnz(mss,numSS);
  n += nnz(mComplex,numTotal-numSS);

  // Allocate memory for this new row
  (*Lambda)[LambdaIndex] = (int *) malloc((2*n+1) * sizeof(int));

  // Enter the number of entries there are
  (*Lambda)[LambdaIndex][0] = n;

  // Enter the non-zero entries
  k = 1;
  for (i = 0; i < numSS; i++) {
    if (mss[i] > 0) {
      (*Lambda)[LambdaIndex][k++] = i;
      (*Lambda)[LambdaIndex][k++] = mss[i];
    }
  }
  for (j = numSS; j < numTotal; j++) {
    if (mComplex[j-numSS] > 0) {
      (*Lambda)[LambdaIndex][k++] = j;
      (*Lambda)[LambdaIndex][k++] = mComplex[j-numSS];
    }
  }

  free(mss);

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int NegCheck(int *mComplex, int *m0,int **A,int numSS,int numTotal) {

  /*
    Check to see if a given set of complexes in the solutions results in
    having a negative number of single-strand species.  This is illegal.

    Returns 0 if negative (illegal) and 1 if not (legal).
  */

  int i,j; // Counters
  double dotprod; // Dot product of row of A with m.

  for (i = 0; i < numSS; i++) {
    dotprod = 0;
    for (j = numSS; j < numTotal; j++) {
      dotprod += A[i][j]*mComplex[j-numSS];
    }
    if (m0[i] - dotprod < 0) {
      return 0;
    }
  }

  return 1;

}
/* ******************************************************************************** */


