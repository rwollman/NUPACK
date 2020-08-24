/*
  DistributionsHeaderFile.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  Header file for use with distributions.c and related
  functions. Contains error messages and function prototypes, among
  other things.
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

// Where the help file to print is
#define DISTRIBUTIONS_HELP_FILE "src/thermo/distributions/Distributions.help"

// Error codes
#define ERR_NOINPUT 2 // No prefix for a filename given
#define ERR_HELP 3 // User chose to display help
#define ERR_HELP_OPEN 4 // Error in opening README file to print help
#define ERR_COUNT 5 // Error opening .count input file
#define ERR_CX 6 // Error in opening .cx file
#define ERR_NONSEQUENTIAL 7 // CX file has nonsequential complex IDs
#define ERR_NOPERMS 8 // Free energies are too high to allow calcs with perms
#define ERR_BADROWINP 9 // Bad row in stoich. matrix in input file
#define ERR_DUPROWA 10 // Row of A is duplicated (indicates bug, not input err)
#define ERR_LOG 11 // Error opening log file
#define ERR_DIST 12 // Error opening .dist file
#define ERR_LAMBDATOOBIG 13 // Too many entries in lambda for method = 2
#define ERR_QBOXTOOBIG 14 // Overflow error in calculation of Q_{box}
#define ERR_LAMBDA 15 // Error opening .lam file
#define ERR_NONINTEGER 16 // Non integer input for initial single-strand count
#define ERR_NEGATIVEPROB 17 // Negative probability encountered

// Function prototypes
/* ******************************* IN CALCDIST.C ********************************** */
void CalcDist(double *mEq, double **Pmn, int **A, double *G, int *m0, double M, 
	      int numSS, int numTotal, double MaxSizeLambda, char *LambdaFile, 
	      double kT, int WriteLambda, int *CompIDArray, int *PermIDArray,
	      int quiet, char *logFile, int WriteLogFile, int NoPermID);
int next(int *mComplex, int *m0, int **A, int numSS, int numComplex, int numTotal,
	 int *mMax, int LastInc);
void UpdateLambda(int ***Lambda, int k, int *mComplex, int *m0, int **A, int numSS,
		  int numTotal);
int NegCheck(int *mComplex, int *m0,int **A,int numSS,int numTotal);
/* ******************************************************************************** */

/* ************************ IN READCOMMANDLINE.C ********************************** */
void ReadCommandLine(int nargs, char **args, char *cxFile, char *countFile, 
		     char *logFile, char *distFile, char *LambdaFile, 
		     int *SortOutput, int *WriteLambda, double *MaxSizeLambda,
		     double *kT, int *quiet, int *WriteLogFile, int *Toverride,
		     int *NoPermID, int * NUPACK_VALIDATE);
void DisplayHelp(int DummyArgument);
/* ******************************************************************************** */

/* ************************ IN INPUTFILEREADER.C ********************************** */
void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
	     int **numPermsArray, int *nComments, char *cxFile, char *countFile, 
	     int quiet);
void ReadInputFiles(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
		    int **m0, double *M, int *numSS, int *numSS0, int *numTotal,
		    int *numPermsArray, char *cxFile, char *countFile, double *kT,
		    int Toverride, char *logFile, char *eqFile, int quiet, 
		    int WriteLogFile);
void ReadInputFilesPerm(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
			int **m0, double *M, int *numSS, int *numSS0, 
			int *newnTotal, int nTotal, char *cxFile, char *countFile, 
			double *kT, int Toverride, char *logFile, char  *eqFile, 
			int quiet, int WriteLogFile);
int InputCompare(const void *p1, const void *p2);
int InputComparePerm(const void *p1, const void *p2);
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
int choldc(double **a, int n, double p[]); // Cholesky decomposition (Num. Rec.)
void cholsl(double **a, int n, double p[], double b[], double x[]);  // Cholesky solve
double factorial(int n);  // Computes n!
double str2double (char *str);  // Converts a string to a double
double WaterDensity(double T); // Computes density of water
/* ******************************************************************************** */

/* **************************** IN WRITEOUTPUT.C ********************************** */
void WriteOutput(double *mEq, double **Pmn, double *G, int **A, int *CompIDArray, 
		 int *PermIDArray, int LargestCompID, int numSS, int numTotal, 
		 int nTotal, int nComments, int maxm0, double kT, char *cxFile, 
		 int SortOutput, char *distFile, int quiet, int NoPermID,int NUPACK_VALIDATE);
int Compare21(const void *p1, const void *p2);
int Compare22(const void *p1, const void *p2);
int Compare23(const void *p1, const void *p2);
int Compare24(const void *p1, const void *p2);
int Compare25(const void *p1, const void *p2);
int Compare26(const void *p1, const void *p2);
int Compare27(const void *p1, const void *p2);
/* ******************************************************************************** */


#ifdef __cplusplus
}
#endif 
