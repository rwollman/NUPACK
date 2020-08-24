/*
  ConcentrationsHeaderFile.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.

  Header file with global variables and function prototypes for use
  with concentrations.c Contains trust region parameters, among
  others.
*/

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <getopt.h>// Takes options from the command line
#include "physical_constants.h"
#include "runtime_constants.h"

// Constants used in trust region
#define TRUST_REGION_DELTABAR 1000.0 // Maximal size of trust region
#define TRUST_REGION_ETA 0.125 // Decision criterion for trust region

#define MAXLOGX 250 // Maximum logarithm of a concentration (prevents overflow)

// Where the help file to print is
#define CONCENTRATIONS_HELP_FILE "src/thermo/concentrations/Concentrations.help"

// Error codes
#define ERR_NOCONVERGE 1 // Error code for failure to converge for method = 1
#define ERR_OVERFLOW 2 // Overflow in calculation of mole fractions
#define ERR_NOINPUT 3 // No prefix for a filename given
#define ERR_HELP 4 // User chose to display help
#define ERR_HELP_OPEN 5 // Error in opening README file to print help
#define ERR_CON 6 // Error opening .con input file
#define ERR_CX 7 // Error in opening .cx file
#define ERR_NONSEQUENTIAL 8 // CX file has nonsequential complex IDs
#define ERR_NOPERMS 9 // Free energies are too high to allow calcs with perms
#define ERR_BADROWINP 10 // Bad row in stoich. matrix in input file
#define ERR_DUPROWA 11 // Row of A is duplicated (indicates bug, not input err)
#define ERR_LOG 12 // Error opening log file
#define ERR_EQ 13 // Error opening .eq file
#define ERR_FPAIRS 14 // Error opening .fpairs file
#define ERR_NOSEQEQ 15 // No sequence information in comments of eq file
#define ERR_PAIRSFILE 16 // Error opening pairs file

// Function prototypes
/* ******************************* IN CALCCONC.C ********************************** */
void getInitialGuess(double *x0, double *lambda, double *G, int **AT, int **A,
		     int numSS, int numTotal, double PerturbScale, 
		     unsigned long rand_seed);
int getx(double *x, double *lambda, double *G, int **AT, int numSS, int numTotal);
void getGrad(double *Grad, double *x0, double *x, int **A, int numSS, int numTotal);
void getHes(double **Hes, double *x, int **A, int numSS, int numTotal);
int getSearchDir(double *p, double *Grad, double **Hes, double delta, int numSS);
double getRho(double *lambda, double *p, double *Grad, double *x, double **Hes, 
	      double *x0, double *G, int **AT, int numSS, int numTotal);
void getCauchyPoint(double *CauchyPoint, double **Hes, double *Grad, double delta, 
		    int numSS);
void PerturbLambda(double *lambda, double PerturbScale, double *G, int **AT, 
		   int numSS, int numTotal);
int CheckTol(double *Grad, double *AbsTol, int numSS);
int CalcConc(double *x, int **A, double *G, double *x0, int numSS, int numTotal, 
	     int MaxIters, double tol, double deltaBar, double eta, double kT, 
	     int MaxNoStep, int MaxTrial, double PerturbScale, int quiet, 
	     int WriteLogFile, char *logFile, double MolesWaterPerLiter, 
	     unsigned long seed);
/* ******************************************************************************** */

/* ************************ IN READCOMMANDLINE.C ********************************** */
void ReadCommandLine(int nargs, char **args, char *cxFile, char *conFile, 
		     char *logFile, char *eqFile, char *pairsFile, char *fpairsFile,
		     int *SortOutput, int *MaxIters, double *tol, double *kT,
		     int *MaxNoStep, int *MaxTrial, double *PerturbScale, int *quiet,
		     int *WriteLogFile, int *Toverride, int *NoPermID, 
		     int *DoBPfracs, unsigned long *seed, double *cutoff, int * NUPACK_VALIDATE);
void DisplayHelpConc(void);
/* ******************************************************************************** */

/* ************************ IN INPUTFILEREADER.C ********************************** */
void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
	     int **numPermsArray, char *cxFile, char *conFile, int quiet);
double ReadInputFiles(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
		      double **x0, int *numSS, int *numSS0, int *numTotal, 
		      int *numPermsArray, char *cxFile, char *conFile, double *kT, 
		      int Toverride, char  *logFile, char  *eqFile, 
		      char *fpairsFile, int quiet, int WriteLogFile, int DoBPfrac,
		      int NoPermID);
/* ******************************************************************************** */

/* *********************************** IN UTILS.C ********************************* */
double min2(double a, double b); // Returns minimum of two arguments
double max2(double a, double b); // Maximum of two arguments
double min(double *ar, int len);  // Minimum of array of doubles
double max(double *ar, int len); // Maximum entry in an array of doubles
int maxint(int *ar, int len); //  Maximum entry in array of ints
double maxabs(double *ar, int len); // Maximum absolute value of array of doubles
int nnz(int *ar, int len); // Number of non-zero elements in an array of ints
int FindNonZero(int *ar, int len); // Index of first nonzero entry
double sum(double *ar, int len);  // Sum of an array of doubles
int sumint(int *ar, int len); // Sum of an array of ints
double dot(double *v1, double *v2, int len);  // dot product of two double arrays
double didot(double *v1, int *v2, int len);  // dot product of int array with double
double norm(double *ar, int len);  // 2-norm of a vector of doubles
void IntTranspose(int **At, int **A, int nrowA, int ncolA);  // matrix transpose (int)
void SymMatrixMult(double **C, double **A, double **B, int n); // multiplication of 
                                                               // symmetric matrices
void MatrixVectorMult(double *c, double **A, double *b, int n); // Maxtix-vector prod.
int choleskyDecomposition(double **a, int n); // Cholesky decomposition
void choleskySolve(double **a, int n, double b[], double x[]);  // Cholesky solve
double factorial(int n);  // Computes n!
double str2double (char *str);  // Converts a string to a double
unsigned long GetRandSeed(unsigned long s); // seed for the random number generator
double WaterDensity(double T); // Computes density of water
/* ******************************************************************************** */

/* **************************** IN WRITEOUTPUT.C ********************************** */
void WriteOutput(double *x, double *G, int *CompIDArray, int LargestCompID, 
		 int numSS, int numTotal, int nTotal, double kT, char *cxFile, 
		 int SortOutput, char  *eqFile, double MolesWaterPerLiter, int quiet, 
		 int NoPermID,int NUPACK_VALIDATE);
int Compare11(const void *p1, const void *p2);  // Comparison function for sorting
int Compare12(const void *p1, const void *p2);  // Comparison function for sorting
int Compare13(const void *p1, const void *p2);  // Comparison function for sorting
int Compare14(const void *p1, const void *p2);  // Comparison function for sorting
/* ******************************************************************************** */

/* **************************** IN FRACPAIR.C ************************************* */
void FracPair(int numSS, int nTotal, int quiet, int NoPermID, int LargestCompID, 
	      int *numPermsArray, char *eqFile, char *conFile, char *pairsFile, 
	      char *fpairsFile, double cutoff,int NUPACK_VALIDATE);
/* ******************************************************************************** */

/* ***************************** IN MT19937AR.C *********************************** */
// Random number generation by Mersenne Twister
void init_genrand(unsigned long s);
double genrand_real1(void);
/* ******************************************************************************** */

#ifdef __cplusplus
}
#endif 
