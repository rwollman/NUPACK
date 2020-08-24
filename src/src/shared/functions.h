#ifndef NUPACK_SHARED_FUNCTIONS_H__
#define NUPACK_SHARED_FUNCTIONS_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "externals.h"
#include "constants.h"
#include "structs.h"


double WaterDensity(double T);
double str2double (char *str);
double min(double *ar, int len);
double max(double *ar, int len);
int maxint(int *ar, int len);
double maxabs(double *ar, int len);
int nnz(int *ar, int len);
int FindNonZero(int *ar, int len);
double sum(double *ar, int len);
int sumint(int *ar, int len);
double dot(double *v1, double *v2, int len);
double didot(double *v1, int *v2, int len);
double norm(double *ar, int len);
void IntTranspose(int **At, int **A, int nrowA, int ncolA);
void SymMatrixMult(double **C, double **A, double **B, int n);
void MatrixVectorMult(double *c, double **A, double *b, int n);
int choleskyDecomposition(double **A, int n);
void choleskySolve(double **A, int n, double *b, double *x);
void lowerTriSolve(double **L, int n, double *b, double *x);
void upperTriSolve(double **U, int n, double *b, double *x);
double min2(double a, double b);
double max2(double a, double b);
int gcd(int a, int b);
double factorial(int n);
long double factorial_long(int n);
int binomial_coefficient(int n, int k); //compute binomial coefficient. Returns -1 on overflow.
unsigned long GetRandSeed(unsigned long s);
int fileExists(char* file);
void SetExecutionPath(int nargs, char **args);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_SHARED_FUNCTIONS_H__ */
