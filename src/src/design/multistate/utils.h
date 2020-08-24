#pragma once

/*
  utilsHeader.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  UTILSHEADER.H

  Header file for use with utility functions included in NUPACK.  These
  utility functions are generate routines for handling vectors and
  matrices, random numbers, factorials, etc.
*/

#include "constants.h"

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

/* **************************** INDEXING MACRO ****************************** */
#define ij(i, j, n) ((i) * (n) + (j))   // Lex index for 2d array, n = # COLUMNS
/* ************************************************************************** */

#ifdef __cplusplus
extern "C" {
#endif

/* ************************* FUNCTION PROTOTYPES **************************** */
double water_density(double T); // Computes density of water
double str2double (char *str);  // Converts a string to a double
double min(double *ar, int len);  // Minimum of array of doubles
double max(double *ar, int len); // Maximum entry in an array of doubles
int sgn(double a); // Sign of a
void printf_matrix_double(double *a, int m, int n); // Print matrix to screen
int maxint(int *ar, int len); //  Maximum entry in array of ints
double maxabs(double *ar, int len); // Maximum abse value of array of doubles
int nnz(int *ar, int len); // Number of non-zero elements in an array of ints
int find_non_zero(int *ar, int len); // Index of first nonzero entry
double sum(double *ar, int len);  // Sum of an array of doubles
int sumint(int *ar, int len); // Sum of an array of ints
double dot(double *v1, double *v2, int len);  // dot prod of two double arrays
double didot(double *v1, int *v2, int len);  // dot prod of int with double
double norm(double *ar, int len);  // 2-norm of a vector of doubles
void transpose_int(int *At, int *A, int nrowA, int ncolA);  // transpose (int)
void sym_matrix_mult(double *C, double *A, double *B, int n); // mult sym mat
void matrix_vector_mult(double *c, double *A, double *b, int m, int n);
void matrix_matrix_mult(double *C, double *A, double *B, int m, int n, int p);
int cholesky_decomposition(double *a, int n); // Cholesky decomposition
int modified_cholesky(double *a, int * p, int n); // Modified Chol. decomp.
int lu_decomposition(double *A, int n); // LU decomposition
void lu_solve(double *a, int n, double b[], double x[]);  // LU solve
int get_approx_rank(double *A, int m, int n, double tol); 
// Get rank of matrix using Gaussian elimination
void lup_decomposition(double *A, int * p, int n); // LUP decomposition
int lup_solve(double *A, int * p, int n, double * x, double * b); // LUP solve
void cholesky_solve(double *a, int n, double b[], double x[]);  // Chol solve
int lower_tri_solve(double *L, int n, double *b, double *x);
int upper_tri_solve(double *U, int n, double *b, double *x);
int modified_cholesky_solve(double *a, int * p, int n, 
        double b[], double x[]);  // Modified Chol. solve
double min2(double a, double b); // Returns minimum of two arguments
double max2(double a, double b); // Maximum of two arguments
unsigned long get_random_seed(unsigned long s); // seed for random number gen
/* ************************************************************************** */
#ifdef __cplusplus
}
#endif
