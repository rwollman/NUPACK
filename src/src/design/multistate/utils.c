/*
  Set of utility functions, including string, array, and matrix
  processing.  Adapted from utils.c from NUPACK.

  utils.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois, 1/2007, except where noted.

  UTILS.C 

  General utility functions for processing arrays (e.g., dot products,
  matrix multiplication, etc.), numbers (e.g., factorials, random
  number generation), strings, etc.  There is one physical function in
  here, which is to compute the density of water at a given
  temperature.

  For use with NUPACK.
*/

#include "utils.h"


/* ************************************************************************** */
double water_density(double T) {
  /* 
     Calculates the number of moles of water per liter at temperature
     (T in degrees C).

     Density of water calculated using data from:
     Tanaka M., Girard, G., Davis, R., Peuto A.,
     Bignell, N.   Recommended table for the denisty
     of water..., Metrologia, 2001, 38, 301-309
  */
   double a1 = -3.983035;
   double a2 = 301.797;
   double a3 = 522528.9;
   double a4 = 69.34881;
   double a5 = 999.974950;

   return a5 * (1 - (T + a1) * (T + a1) * (T + a2) / a3 / (T + a4)) / 18.0152;
}
/* ************************************************************************** */

/*                          BEGIN ARRAY PROCESSING FUNCTIONS                  */
/* ************************************************************************** */

/* ************************************************************************** */
void printf_matrix_double(double *a, int m, int n) {
  /*
    Prints an mxn matrix of double to the screen.
  */

  int i, j; // Indices
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      printf("%g\t", a[ij(i, j, n)]);
    }
    printf("\n");
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
int find_non_zero(int *ar, int len) {
  /*
    Returns the index of the first nonzero entry in an array of ints
    of length len.

    Returns -1 if all entries are zero.
  */

  int i; // Counter

  for (i = 0; i < len; i++) {
    if (ar[i] != 0) {
      return i;
    }
  }

  return -1;

}
/* ************************************************************************** */

/* ************************************************************************** */
void transpose_int(int *At, int *A, int nrowA, int ncolA) {
  /* 
     Puts the transpose of matrix A into matrix At.
     A have nrowA rows and ncolA columns.
     The matrix is full of integers.
     At must be preallocated to have ncolA rows and nrowA columns.
  */

  int i,j; // Counters
  
  for (i = 0; i < nrowA; i++) {
    for (j = 0; j < ncolA; j++) {
      At[ij(j, i, nrowA)] = A[ij(i, j, ncolA)];
    }
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
void sym_matrix_mult(double *C, double *A, double *B, int n) {
  /*
    Performs matrix multiplication of A*B where A and B are symmetric
    n x n matrices.  All entries must be present in A and B, not just
    the upper or lower triangle.  The matrices are arrays of doubles.
    Matrix C holds the result and must be pre-allocated.
  */

  int i,j,k; // Counters
  double *v;  // Column in B that we dot with a row in A to get entry in product

  v = (double *) malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      C[ij(i, j, n)] = 0.0;
      for (k = 0; k < n; k++) {
        C[ij(i, j, n)] += A[ij(i, k, n)] * B[ij(k, j, n)];
      }
    }
  }

  // Fill out the lower triangle
  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      C[ij(i, j, n)] = C[ij(j, i, n)];
    }
  }

  free(v);
}
/* ************************************************************************** */


/* ************************************************************************** */
void matrix_vector_mult(double *c, double *A, double *b, int m, int n) {
  /*
    Performs the multiplication of the m x n matrix A with n-vector b
    and stores the result in vector c, which must be pre-allocated to
    have length m.

    All entries are doubles.
  */

  int i, j; // Indices

  for (i = 0; i < m; i++) {
    c[i] = 0.0;
    for (j = 0; j < n; j++) {
      c[i] += A[ij(i, j, n)] * b[j];
    }
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
void matrix_matrix_mult(double *C, double *A, double *B, int m, int n, int p) {
  /*
    Performs the multiplication of the m x n matrix A with n x p
    matrix B and stores the result in m x p matrix C, which must be
    pre-allocated to have m x p entries.

    All entries are doubles.
  */

  int i, j, k; // Indices

  for (i = 0; i < m; i++) {
    for (j = 0; j < p; j++) {
      C[ij(i, j, p)] = 0.0;
      for (k = 0; k < n; k++) {
        C[ij(i, j, p)] += A[ij(i, k, n)] * B[ij(k, j, p)];
      }
    }
  }
}
/* ************************************************************************** */

/* ************************************************************************** */
int modified_cholesky(double *A, int * p, int n) {
  /* 
   * Performs Modified Cholesky decomposition based on the algorithm 
   * GMW81 in Fang, O'Leary, 2006. Modified Cholesky Algorithms: A Catalog 
   * with New Approaches. We can replace this later with a more complicated
   * version if we want, but this should work since all our matrices
   * are extremely close to positive definite even after numerical error
   *
   * Only the lower triangular of L is computed and stored in A.
   *
   * This is almost a drop-in replacement of cholesky decomposition with the
   * addition of the permutation vector stored in p
   */
  int i, j, k;
  double beta;
  double c_sum;
  double mu_val;
  int mu;
  double temp;
  double eta = 0;
  double xi = 0;

  for (i = 0; i < n; i++) {
    p[i] = i;
  }

  if (n <= 0) {
    return 0;
  } 

  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      temp = fabs(A[ij(i, j, n)]);
      xi = max2(xi, temp);
    }
  }
  for (i = 0; i < n; i++) {
    temp = fabs(A[ij(i, i, n)]);
    eta = max2(eta, temp);
  }

  if (n > 1) {
    beta = sqrt(max2(eta, xi / sqrt(n*n - 1)));
  } else {
    beta = sqrt(eta);
  }
  beta = max2(beta, FLOAT_EPS);
  
  for (k = 0; k < n; k++) {
    // Pick a pivot
    mu_val = A[ij(k, k, n)];
    mu = k;
    for (i = k + 1; i < n; i++) {
      temp = A[ij(i, i, n)];
      if (mu_val < temp) {
        mu = i;
        mu_val = temp;
      }
    }
    // diagonal pivot k <=> mu
    i = p[mu];
    p[mu] = p[k];
    p[k] = i;

    for (i = 0; i < k; i++) {
      temp = A[ij(k, i, n)];
      A[ij(k, i, n)] = A[ij(mu, i, n)];
      A[ij(mu, i, n)] = temp;
    }
    temp = A[ij(k, k, n)];
    A[ij(k, k, n)] = A[ij(mu, mu, n)];
    A[ij(mu, mu, n)] = temp;
    for (i = k + 1; i < mu; i++) {
      temp = A[ij(i, k, n)];
      A[ij(i, k, n)] = A[ij(mu, i, n)];
      A[ij(mu, i, n)] = temp;
    }
    for (i = mu + 1; i < n; i++) {
      temp = A[ij(i, k, n)];
      A[ij(i, k, n)] = A[ij(i, mu, n)];
      A[ij(i, mu, n)] = temp;
    }

    // Compute c_sum
    c_sum = 0;
    for (i = k + 1; i < n; i++) {
      c_sum = max2(c_sum, fabs(A[ij(i, k, n)]));
    }
    c_sum /= beta;
    c_sum = c_sum * c_sum;

    temp = fabs(A[ij(k, k, n)]);
    temp = max2(temp, NUM_PRECISION * eta);
    temp = max2(temp, c_sum);
    A[ij(k, k, n)] = sqrt(temp);

    // Compute the current column of L
    for (i = k + 1; i < n; i++) {
      A[ij(i, k, n)] /= A[ij(k, k, n)];
    }

    // Adjust the \bar{A} 
    for (j = k+1; j < n; j++) {
      for (i = j; i < n; i++) {
        A[ij(i, j, n)] -= A[ij(i, k, n)] * A[ij(j, k, n)];
      } 
    }
  }
  return 0;
}
/* ************************************************************************** */

/* ************************************************************************** */
int cholesky_decomposition(double *A, int n) {
  /*
    Performs Cholesky decomposition on the positive-definite symmetric
    n by n matrix A.  Only the lower triangle of A need be supplied.
    Because only the lower triangle is used, A is assumed to be
    symmetric.  Symmetry is NOT checked.  

    The Cholesky decomposition is A = L L^T.  The lower triangular
    matrix (including its diagonal) is written in A.  I.e., A[i][j] =
    L[i][j] for i >= j.

    This is algorithm is the outer product Cholesky decomposition
    algorithm without pivoting.  We use algorithm 4.2.2 in Golub and
    van Loan.

    Returns 0 if Cholesky decompostion was successful and 1 if it fails.
  */

  int i,j,k;
  
  for (k = 0; k < n; k++) {
    if (A[ij(k, k, n)] <= 0.0) { // Cholesky decomposition failed, not pos def.
      return 1;
    }
    A[ij(k, k, n)] = sqrt(A[ij(k, k, n)]);
    for (i = k+1; i < n; i++) {
      A[ij(i, k, n)] /= A[ij(k ,k, n)];
    }
    for (j = k+1; j < n; j++) {
      for (i = j; i < n; i++) {
        A[ij(i, j, n)] -= A[ij(i, k, n)] * A[ij(j, k, n)];
      } 
    }
  }

  // Cholesky decomposition was successful, return 1
  return 0;
}
/* ************************************************************************** */

/* ************************************************************************** */
int lu_decomposition(double *A, int n) {
  /*
    Performs LU decomposition on a nonsingular square matrix A.

    The LU decomposition is A = L U.  The lower triangular matrix is
    written in A.  I.e., A[i][j] = L[i][j] for i > j.  The diagonal of
    L is all ones.  The upper triangular matrix is also written in A,
    including its diagonal.  I.e., A[i][j] = U[i][j] for i <= j.

    This is algorithm is the outer product LU decomposition by
    Gaussian eliminaryion without pivoting.  We use algorithm 3.2.1 in
    Golub and van Loan.

    Returns 0 if LU decompostion was successful and 1 if the matrix is
    found to be singular.
   */

  int i, j, k;
  
  for (k = 0; k < n - 1; k++) {
    // Check for singularity and return 0 if failed
    if (fabs(A[ij(k, k, n)]) < NUM_PRECISION) {
      return 1;
    }
    for (i = k + 1; i < n; i++) {
      A[ij(i, k, n)] /= A[ij(k, k, n)];
      for (j = k + 1; j < n; j++) {
        A[ij(i, j, n)] -= A[ij(i, k, n)] * A[ij(k, j, n)];
      }
    }
  }

  // LU decomposition was successful
  return 0;
}
/* ************************************************************************** */

/* ************************************************************************** */
void lu_solve(double *A, int n, double *b, double *x) {
  /*
   Solves the system of n linear equations Ax = b by Cholesky
   decomposition, where A = L U, where L and U have already been found
   by LU decomposition.  L is stored in the lower triangle of A (its
   diagonal is all ones).  U (with diagonal) is stored in the upper
   triangle of A. I.e., A HAS ALREADY UNDERGONE LU DECOMPOSITION using
   the function lu_decomposition.

   This function is to be used with the lu_decomposition function.
   The solution is returned in x, which must be preallocated of length
   n.

   Note that solving Ly = b and then U x = y amounts to solving Ax =
   b. To solve the first lower-triangular matrix equation, we use
   column-based forward substitution, outlined in algorithm 3.1.3 of
   Golub and Van Loan.  To solve the second upper-triangular matrix
   equation (U x = y), we use colum-based back substitution,
   algorithm 3.1.4 in Golub and Van Loan.
  */

  int i,j;
  double *L; // Lower triangular matrix made from A and diagonal of ones

  // Allocate and compute L
  L = (double *) malloc(n * n * sizeof(double));
  for (i = 0; i < n; i++) {
    L[ij(i, i, n)] = 1.0;
    for (j = 0; j < i; j++) {
      L[ij(i, j, n)] = A[ij(i, j, n)];
    }
  }

  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Ly = b, storing y in x.
  lower_tri_solve(L, n, b, x);

  // Solve Ux = y by back substitution
  upper_tri_solve(A, n, x, x);

  // Free L
  free(L);

}
/* ************************************************************************** */


/* ************************************************************************** */
int get_approx_rank(double * A, int m, int n, double tol) {
  /* 
   * Finds the rank of m x n matrix A using Gaussian elimination
   * I just need a fast way of calculating the approximate rank of A.
   * returns rank on success, -1 on error
   */
  
  int i, j, k;
  int ret_val = 0;
  int * p = NULL;
  double * A_copy = NULL;
  double tmp;
  int mu = 0;
  double mu_val = 0;

  int new_m;
  int new_n;
  if (m > n) {
    new_m = m;
    new_n = n;
  } else {
    new_m = n;
    new_n = m;
  }

  p = malloc(new_m * sizeof(int));
  A_copy = malloc(m * n * sizeof(double));
  if (p == NULL || A_copy == NULL) {
    ret_val = -1;
    goto end_approx_rank;
  }

  if (m > n) {
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        A_copy[ij(i, j, n)] = A[ij(i, j, n)];
      }
    }
  } else {
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        A_copy[ij(j, i, m)] = A[ij(i, j, n)];
      }
    }
  }

  for (i = 0; i < m; i++) {
    p[i] = i;
  }

  for (k = 0; k < new_n; k++) {
    // Find pivot
    mu = k;
    mu_val = fabs(A_copy[ij(k, k, new_n)]);
    for(j = k+1; j < new_m; j++) {
      if(fabs(A_copy[ij(j, k, new_n)]) > mu_val) {
        mu = j;
        mu_val = fabs(A_copy[ij(j, k, new_n)]);
      }
    }
    // Pivot
    for(j = 0; j < new_n; j++) {
      tmp = A_copy[ij(k, j, new_n)];
      A_copy[ij(k, j, new_n)] = A_copy[ij(mu, j, new_n)];
      A_copy[ij(mu, j, new_n)] = tmp;
    }

    i = p[k];
    p[k] = p[mu];
    p[mu] = i;
    if (mu_val > tol) {
      ret_val++;
      for(i = k+1; i < new_m; i++) {
        A_copy[ij(i, k, new_n)] = A_copy[ij(i, k, new_n)] / A_copy[ij(k, k, new_n)];
      }
      for(i = k+1; i < new_m; i++) {
        for(j = k+1; j < new_n; j++) {
          A_copy[ij(i, j, new_n)] = A_copy[ij(i, j, new_n)] - 
            A_copy[ij(i, k, new_n)]*A_copy[ij(k, j, new_n)];
        }
      }
    }
  }

end_approx_rank:
  if (A_copy) free(A_copy);
  if (p) free(p);

  return ret_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
void lup_decomposition(double *A, int *p, int n) {
  /* 
    Performs LUP decomposition of square matrix A.

    The LUP decomposition is P.A = L.U  
    The lower triangle A[i][j]: i > j is set to L[i][j]
    The diagonal of L is always unity
    The upper triangle including the diagonal: A[i][j]: i <= j is set to U[i][j]
    p[i] = j: P[i][j] = 1 in the above expression

    A should be n x n, p should be n elements.

    This algorithm is based on Algorithm 3.4.1 in Golub and van Loan with 
    extension from the end of 3.4.4

    Does not return any errors.
  */
  int i, j, k, mu;

  double m_val; // used to store max value
  double tmp;   // used to store temporary value during pivoting
  for(k = 0; k < n; k++) {
    p[k] = k;
  }

  for(k = 0; k < n-1; k++) {
    mu = k;
    m_val = fabs(A[ij(k, k, n)]);
    // Find pivot
    for(j = k+1; j < n; j++) {
      if(fabs(A[ij(j, k, n)]) > m_val) {
        mu = j;
        m_val = fabs(A[ij(j, k, n)]);
      }
    }
    // Pivot
    for(j = 0; j < n; j++) {
      tmp = A[ij(k, j, n)];
      A[ij(k, j, n)] = A[ij(mu, j, n)];
      A[ij(mu, j, n)] = tmp;
    }
    i = p[k];
    p[k] = p[mu];
    p[mu] = i;
    
    if(m_val > NUM_PRECISION) {
      // If we have a non-zero pivot
      for(i = k+1; i < n; i++) {
        A[ij(i, k, n)] = A[ij(i, k, n)] / A[ij(k, k, n)];
      }
      for(i = k+1; i < n; i++) {
        for(j = k+1; j < n; j++) {
          A[ij(i, j, n)] = A[ij(i, j, n)] - A[ij(i, k, n)]*A[ij(k, j, n)];
        }
      }
    }
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
int lup_solve(double *A, int *p, int n, double *b, double *x) {
  /*
    Performs LUP solve of square matrix A which has already undergone LUP
    decomposition (see lup_decomposition).
    
    This follows from Golub and van Loan (3.4.9)

    returns 0 on success
    returns 1 on memory error
   */
  int i,j;
  double *L = NULL;      // temporary L for lower_tri_solve
  double *b_perm = NULL; // permuted b
  int ret_val = 0;// return value

  if(NULL == (L = (double *) malloc(n * n * sizeof(double)))) {
    ret_val = 1;
    goto end_lup_solve;
  }
  if(NULL == (b_perm = (double *) malloc(n*sizeof(double)))) {
    free(L);
    ret_val = 1;
    goto end_lup_solve;
  }

  for(i = 0; i < n; i++) {
    L[ij(i, i, n)] = 1.0;
    for(j = 0; j < i; j++) {
      L[ij(i, j, n)] = A[ij(i, j, n)];
    }
  }

  for(i = 0; i < n; i++) {
    b_perm[i] = b[p[i]];
  }

  // Solve Ly = b_perm, storing y in x.
  ret_val = lower_tri_solve(L, n, b_perm, x);

  // Solve Ux = y by back substitution
  ret_val += upper_tri_solve(A, n, x, x);

end_lup_solve:
  if(L != NULL) free(L);
  if(b_perm != NULL) free(b_perm);

  return ret_val;
}
/* ************************************************************************** */

/* ************************************************************************** */
int modified_cholesky_solve(double *A, int * p, int n, double *b, double *x) {
  /*
  */

  int i,j;
  double *U = NULL; // Upper triangular matrix made from L^T
  double *xp = NULL;
  int ret_val = 0;

  // Allocate and compute U
  U = (double *) malloc(n * n * sizeof(double));
  xp = (double *) malloc(n * sizeof(double));

  if (U == NULL || xp == NULL) {
    ret_val = 1;
    goto end_mod_cholesky;
  }

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      U[ij(i, j, n)] = A[ij(j, i, n)];
    }
  }

  // Copy b to x
  for (i = 0; i < n; i++) {
    xp[i] = b[p[i]];
  }

  // Solve Ly = b, storing y in xp.
  lower_tri_solve(A, n, xp, x);

  // Solve L^T x = y by back substitution
  upper_tri_solve(U, n, x, xp);

  for (i = 0; i < n; i++) {
    x[p[i]] = xp[i];
  }

end_mod_cholesky:
  // Free U
  if (U != NULL) free(U);
  if (xp != NULL) free(xp);
  return ret_val;
}
/* ************************************************************************** */



/* ************************************************************************** */
void cholesky_solve(double *A, int n, double *b, double *x) {
  /*
   Solves the system of n linear equations Ax = b by Cholesky
   decomposition, where A = L L^T, where L has already been found by
   Cholesky decomposition.  L is stored in the lower triangle of A.
   I.e., A HAS ALREADY UNDERGONE CHOLESKY DECOMPOSITION.

   This function is to be used with the cholesky_decomposition
   function.  The solution is returned in x, which must be
   preallocated of length n.

   Note that solving Ly = b and then L^T x = y amounts to solving Ax =
   b. To solve the first lower-triangular matrix equation, we use
   column-based forward substitution, outlined in algorithm 3.1.3 of
   Golub and Van Loan.  To solve the second upper-triangular matrix
   equation (L^T x = y), we use colum-based back substitution,
   algorithm 3.1.4 in Golub and Van Loan.
  */

  int i,j;
  double *U; // Upper triangular matrix made from L^T

  // Allocate and compute U
  U = (double *) malloc(n * n * sizeof(double));
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      U[ij(i, j, n)] = A[ij(j, i, n)];
    }
  }

  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Ly = b, storing y in x.
  lower_tri_solve(A, n, b, x);

  // Solve L^T x = y by back substitution
  upper_tri_solve(U, n, x, x);

  // Free U
  free(U);

}
/* ************************************************************************** */

/* ************************************************************************** */
int lower_tri_solve(double *L, int n, double *b, double *x) {
  /*
    Solves the lower triangular system Lx = b.  L is a
    lower-triangular matrix (including diagonal).  We use column-based
    forward substitution, outlined in algorithm 3.1.3 of Golub and van
    Loan.
  */

  int i,j;
  int ret_val = 0;
  
  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Lx = b
  for (j = 0; j < n-1; j++) {
    if (fabs(L[ij(j, j, n)]) > FLOAT_EPS) {
      x[j] /= L[ij(j, j, n)];
      for (i = j+1; i < n; i++) {
        x[i] -= x[j] * L[ij(i, j, n)];
      }
    } else {
      x[j] = 0.0;
      ret_val = 1;
    }
  }
  if(n > 0) {
    if (fabs(L[ij(n-1, n-1, n)]) > FLOAT_EPS) {
      x[n-1] /= L[ij(n-1, n-1, n)];
    } else {
      x[n-1] = 0.0;
      ret_val = 1;
    }
  }
  return ret_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
int upper_tri_solve(double *U, int n, double *b, double *x) {
  /*
    Solves the lower triangular system Ux = b.  U is an
    upper-triangular matrix (including diagonal).  We use column-based
    forward substitution, outlined in algorithm 3.1.4 of Golub and van
    Loan.
  */

  int i,j;
  int ret_val = 0;
  
  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Ux = b by back substitution
  for (j = n-1; j > 0; j--) {
    if (fabs(U[ij(j, j, n)]) > FLOAT_EPS) {
      x[j] /= U[ij(j, j, n)];
      for (i = 0; i < j; i++) {
        x[i] -= x[j] * U[ij(i, j, n)];
      }
    } else {
      x[j] = 0.0;
      ret_val = 1;
    }
  }
  if(n > 0) {
    if (fabs(U[ij(0, 0, n)]) > FLOAT_EPS) {
      x[0] /= U[ij(0, 0, n)];
    } else {
      x[0] = 0.0;
      ret_val = 1;
    }
  }
  return ret_val;
}  
/* ************************************************************************** */


/* ************************************************************************** */
/*                          END ARRAY PROCESSING FUNCTIONS                    */
/* ************************************************************************** */



/* ************************************************************************** */
/*                         BEGIN NUMBER PROCESSING FUNCTIONS                  */
/* ************************************************************************** */

/* ************************************************************************** */
int sgn(double a) {
  /*
    Returns 1 if a is positive, -1 if negative, 0 if zero
  */

  if (a < 0.0) {
    return -1;
  }
  else if (a > 0.0) {
    return 1;
  }
  else {
    return 0;
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
unsigned long get_random_seed(unsigned long s){ 
  /*
    Generates seed for random number generator off of clock.
    If s = 0, uses clock seed, otherwise s is the seed.
  */
  time_t tseed;
  unsigned long seed;

  if(s == 0) { 
    time(&tseed);
    seed = (unsigned long)(tseed);
  }
  else {
    seed = s;
  }

  return seed;

}
/* ************************************************************************** */

/* ************************************************************************** */
/*                         END NUMBER PROCESSING FUNCTIONS                    */
/* ************************************************************************** */
