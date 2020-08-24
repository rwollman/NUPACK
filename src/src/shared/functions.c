/*
  shared.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois, 1/2007, except where noted.

  functions.c

  General utility functions for processing arrays (e.g., dot products,
  matrix multiplication, etc.), numbers (e.g., factorials, random
  number generation), strings, etc.  There is one physical function in
  here, which is to compute the density of water at a given
  temperature.

  For use with NUPACK.
*/

#include "functions.h"


/* ******************************************************************************** */
double WaterDensity(double T) {
  /* 
     Calculates the number of moles of water per liter at temperature
     (T in degrees C).

     Density of water calculated using data from:
     Tanaka M., Girard, G., Davis, R., Peuto A.,
     Bignell, N.   Recommended table for the density
     of water..., Metrologia, 2001, 38, 301-309
  */
   double a1 = -3.983035;
   double a2 = 301.797;
   double a3 = 522528.9;
   double a4 = 69.34881;
   double a5 = 999.974950;


   return a5 * (1 - (T+a1)*(T+a1)*(T+a2)/a3/(T+a4)) / 18.0152;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double str2double (char *str) {
  /* 
     Converts a string to a double.  The string may either be
     an integer or float.  E.g., 56 or 0.45 or 654.234.  It may 
     also be in scientific notation, e.g., 1.5e-12, 43.54e5,
     3.32E-8, 4E4, 2.45e+15, or 8.55E+12.  No spaces are allowed.
     This is easier than having a bunch of sscanf's with if state-
     ments in the main code.
   */

  int i,k; // counters
  int noE; // Haven't encountered an e or E yet.
  char *MantissaStr; // string storing the mantissa
  char *ExpStr; // string storing the exponent
  int Len; // length of string
  double mantissa; // number is mantissa * 10^exponent 
  int exponent; 

  Len = strlen(str);

  noE = 1;
  k = 0;
  while (k < Len && noE) {
    if (str[k] == 'e' || str[k] == 'E') {
      noE = 0;
    }
    k++;
  }

  if (k == Len) { // Not in scientific notation
    return atof(str);
  }

  // k is now the index of the start of the exponent
  ExpStr = (char *) malloc((Len-k+1) * sizeof(char));
  MantissaStr = (char *) malloc(k * sizeof(char));
  strncpy(MantissaStr,str,k-1);
  MantissaStr[k-1] = '\0';

  for (i = 0; i < Len-k; i++) {
    ExpStr[i] = str[k+i];
  }
  ExpStr[Len-k] = '\0';

  mantissa = atof(MantissaStr);
  exponent = atoi(ExpStr);

  free(MantissaStr);
  free(ExpStr);
  
  return (mantissa * pow(10,exponent));

}
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                          BEGIN ARRAY PROCESSING FUNCTIONS                        */
/* ******************************************************************************** */

/* ******************************************************************************** */
double min(double *ar, int len) {
  /* 
     Returns the minimum entry in an array of doubles of length len
  */

  int i; // Counter
  double MinEntry; // The entry we'll return

  MinEntry = ar[0];

  for (i = 1; i < len; i++) {
    if (ar[i] < MinEntry) {
      MinEntry = ar[i];
    }
  }

  return MinEntry;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double max(double *ar, int len) {
  /*
    Returns the maximum entry in an array of doubles of length len.
  */

  double MaxEntry; // The max entry
  int i; // Counter

  MaxEntry = ar[0];
  for (i = 0; i < len; i++) {
    if (ar[i] > MaxEntry) {
      MaxEntry = ar[i];
    }
  }

  return MaxEntry;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int maxint(int *ar, int len) {
  /*
    Returns the maximum entry in an array of integers of length len.
  */

  int MaxEntry; // The max entry
  int i; // Counter

  MaxEntry = ar[0];
  for (i = 0; i < len; i++) {
    if (ar[i] > MaxEntry) {
      MaxEntry = ar[i];
    }
  }

  return MaxEntry;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double maxabs(double *ar, int len) {
  /*
    Returns the maximum absolute value in an array of doubles of length len.
  */

  double MaxEntry; // The max entry
  int i; // Counter

  MaxEntry = fabs(ar[0]);
  for (i = 0; i < len; i++) {
    if (fabs(ar[i]) > MaxEntry) {
      MaxEntry = fabs(ar[i]);
    }
  }

  return MaxEntry;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int nnz(int *ar, int len) {
  /*
    Returns the number of nonzero elements in an array of integers of length len.
  */

  int i; // Counter
  int NumNonZero; // Number of non zere elements.... What we return.

  NumNonZero = 0;

  for (i = 0; i < len; i++) {
    if (ar[i] != 0) {
      NumNonZero++;
    }
  }

  return NumNonZero;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int FindNonZero(int *ar, int len) {
  /*
    Returns the index of the first nonzero entry in an array of ints of length len.

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
/* ******************************************************************************** */


/* ******************************************************************************** */
double sum(double *ar, int len) {
  /*
    Sums the elements in the array of doubles length len.
  */

  int i;
  double sumar;

  sumar = 0.0;
  for (i = 0; i < len; i++) {
    sumar += ar[i];
  }

  return sumar;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int sumint(int *ar, int len) {
  /*
    Sums the elements in the array of integers of length len.
  */

  int i;
  int sumar;

  sumar = 0;
  for (i = 0; i < len; i++) {
    sumar += ar[i];
  }

  return sumar;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double dot(double *v1, double *v2, int len) {
  /*
    Computes dot product of v1 and v2 (v1^T v2) and
    returns the result.

    v1 and v2 must be doubles.  They must have len entries in them.
  */
  
  int i; // Counter
  double dotprod; // The dot product
  
  dotprod = 0.0;
  for (i = 0; i < len; i++) {
    dotprod += v1[i]*v2[i];
  }

  return dotprod;
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double didot(double *v1, int *v2, int len) {
  /*
    Computes dot product of v1 and v2 (v1^T v2) and
    returns the result.

    v1 must be a double and v2 an array of ints.  They
    must have len entries in them.
  */
  
  int i; // Counter
  double dotprod; // The dot product

  dotprod = 0.0;
  for (i = 0; i < len; i++) {
    dotprod += v1[i]* ((double) v2[i]);
  }

  return dotprod;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
double norm(double *ar, int len) {
  /*
    Computes the norm of an array of double of length len.
  */

  return sqrt(dot(ar,ar,len));

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void IntTranspose(int **At, int **A, int nrowA, int ncolA) {
  /* 
     Puts the transpose of matrix A into matrix At.
     A have nrowA rows and ncolA columns.
     The matrix is full of integers.
     At must be preallocated to have ncolA rows and nrowA columns.
  */

  int i,j; // Counters
  
  for (i = 0; i < nrowA; i++) {
    for (j = 0; j < ncolA; j++) {
      At[j][i] = A[i][j];
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void SymMatrixMult(double **C, double **A, double **B, int n) {
  /*
    Performs matrix multiplication of A*B where A and B are symmetric n x n
    matrices.  All entries must be present in A and B, not just the upper
    or lower triangle.  The matrices are arrays of doubles.
    Matrix C holds the result and must be pre-allocated.
  */

  int i,j,k; // Counters
  double *v;  // Column in B that we dot with a row in A to get entry in product

  v = (double *) malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      for (k = 0; k < n; k++) {
	v[k] = B[k][j];
      }
      C[i][j] = dot(A[i],v,n);
    }
  }

  // Fill out the lower triangle
  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      C[i][j] = C[j][i];
    }
  }

  free(v);
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void MatrixVectorMult(double *c, double **A, double *b, int n) {
  /*
    Performs the multiplication of the n x n matrix A with n-vector b
    and stores the result in vector c, which must be pre-allocated.

    All entries are doubles.
  */

  int i; // Counter

  for (i = 0; i < n; i++) {
      c[i] = dot(A[i],b,n);
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int choleskyDecomposition(double **A, int n) {
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


    Returns 1 if Cholesky decompostion was successful and 0 if it fails.
  */

  int i,j,k;
  
  for (k = 0; k < n; k++) {
    if (A[k][k] <= 0.0) { // Cholesky decomposition failed, not pos def.
      return 0;
    }
    A[k][k] = sqrt(A[k][k]);
    for (i = k+1; i < n; i++) {
      A[i][k] /= A[k][k];
    }
    for (j = k+1; j < n; j++) {
      for (i = j; i < n; i++) {
	A[i][j] -= A[i][k]*A[j][k];
      } 
    }
  }

  // Cholesky decomposition was successful, return 1
  return 1;
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void choleskySolve(double **A, int n, double *b, double *x) {
  /*
   Solves the system of n linear equations Ax = b by Cholesky
   decomposition, where A = L L^T, where L has already been found by
   Cholesky decomposition.  L is stored in the lower triangle of A.
   I.e., A HAS ALREADY UNDERGONE CHOLESKY DECOMPOSITION.

   This function is to be used with the choleskyDecomposition
   function.  The solution is returned in x, which must be
   preallocated of length n.

   Note that solving Ly = b and then L^T x = y amounts to solving Ax =
   b. To solve the first lower-triangular matrix equation, we use
   column-based forward substitution, outline in algorithm 3.1.3 of
   Golub and van Loan.  To solve the second upper-triangular matrix
   equation (L^T x = y), we use colum-based back substitution,
   algorithm 3.1.4 in Golub and Van Loan.
  */

  int i,j;
  double **U; // Upper triangular matrix made from L^T

  // Allocate and compute U
  U = (double **) malloc(n * sizeof(double *));
  for (i = 0; i < n; i++) {
    U[i] = (double *) malloc(n * sizeof(double));
    for (j = i; j < n; j++) {
      U[i][j] = A[j][i];
    }
  }

  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Ly = b, storing y in x.
  lowerTriSolve(A,n,b,x);

  // Solve L^T x = y by back substitution
  upperTriSolve(U,n,x,x);

  // Free U
  for (i = 0; i < n; i++) {
    free(U[i]);
  }
  free(U);

}
/* ******************************************************************************** */

/* ******************************************************************************** */
void lowerTriSolve(double **L, int n, double *b, double *x) {
  /*
    Solves the lower triangular system Lx = b.  L is a
    lower-triangular matrix (including diagonal).  We use column-based
    forward substitution, outlined in algorithm 3.1.3 of Golub and van
    Loan.
  */

  int i,j;
  
  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Lx = b
  for (j = 0; j < n-1; j++) {
    x[j] /= L[j][j];
    for (i = j+1; i < n; i++) {
      x[i] -= x[j]*L[i][j];
    }
  }
  x[n-1] /= L[n-1][n-1];

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void upperTriSolve(double **U, int n, double *b, double *x) {
  /*
    Solves the lower triangular system Ux = b.  U is an
    upper-triangular matrix (including diagonal).  We use column-based
    forward substitution, outlined in algorithm 3.1.4 of Golub and van
    Loan.
  */

  int i,j;
  
  // Copy b to x
  for (i = 0; i < n; i++) {
    x[i] = b[i];
  }

  // Solve Ux = b by back substitution
  for (j = n-1; j > 0; j--) {
    x[j] /= U[j][j];
    for (i = 0; i < j; i++) {
      x[i] -= x[j]*U[i][j];
    }
  }
  x[0] /= U[0][0];

}  
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                          END ARRAY PROCESSING FUNCTIONS                          */
/* ******************************************************************************** */



/* ******************************************************************************** */
/*                         BEGIN NUMBER PROCESSING FUNCTIONS                        */
/* ******************************************************************************** */

/* ******************************************************************************** */
double min2(double a, double b) {
  /*
    Returns the min of two doubles, a and b
  */

  if (a < b) {
    return a;
  }
  else {
    return b;
  }
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double max2(double a, double b) {
  /*
    Returns the max of two doubles, a and b
  */

  if (a > b) {
    return a;
  }
  else {
    return b;
  }
}
/* ******************************************************************************** */


/* ******************************************************************************** */
int gcd(int a, int b) {
  /*
    Computes the greatest common divisor of a and b, returned as an int.
  */

  int c;
	
  while(b != 0) {
    c = a % b;
    a = b;
    b = c;
  }
	
  return abs(a);
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double factorial(int n) {
  /*
    Computes n!.  Returns the result as a double.

    Factorials up to 33! are stored and recalled.  Bigger factorials
    are calculated.

    WARNING: This program does not check to see if the argument of the factorial 
    is bigger than the size of the maximum double.  It may therefore produce overflow
    errors.
  */

  if (n < 0) {
    printf("Error!  Negative factorial not defined.\n");
    exit(ERR_FACTORIAL);
  }

  // Factorials
  static const double Fact[] = {
    1.00000000000000000000E0,
    1.00000000000000000000E0,
    2.00000000000000000000E0,
    6.00000000000000000000E0,
    2.40000000000000000000E1,
    1.20000000000000000000E2,
    7.20000000000000000000E2,
    5.04000000000000000000E3,
    4.03200000000000000000E4,
    3.62880000000000000000E5,
    3.62880000000000000000E6,
    3.99168000000000000000E7,
    4.79001600000000000000E8,
    6.22702080000000000000E9,
    8.71782912000000000000E10,
    1.30767436800000000000E12,
    2.09227898880000000000E13,
    3.55687428096000000000E14,
    6.40237370572800000000E15,
    1.21645100408832000000E17,
    2.43290200817664000000E18,
    5.10909421717094400000E19,
    1.12400072777760768000E21,
    2.58520167388849766400E22,
    6.20448401733239439360E23,
    1.55112100433309859840E25,
    4.03291461126605635584E26,
    1.0888869450418352160768E28,
    3.04888344611713860501504E29,
    8.841761993739701954543616E30,
    2.6525285981219105863630848E32,
    8.22283865417792281772556288E33,
    2.6313083693369353016721801216E35,
    8.68331761881188649551819440128E36
  };

  double nFact = 1.0; // the result
  int c; // counter

  if (n <= 33) {
    return Fact[n];
  }
  else {
    nFact = Fact[33];
    for (c = 34; c < nFact; c++) {
      nFact *= c;
    }
    return nFact;
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int binomial_coefficient(int n, int k) {
// Compute n choose k using factorials.
// return an error if there is an overflow.
  int cur_n;
  int cur_k;
  int val = 1;
  int prev_val = 1;
  int lim1 = k;
  int lim2 = n - k;
  if(k > n) {
    return 0;
  } 

  if(lim1 < lim2) {
    int temp = lim1;
    lim1 = lim2;
    lim2 = temp;
  }
  // lim1 is the greater limit to minimize iterations
  // of the for loop.
  // Each iteration of the for loop multiplies by a term
  // from the numerator, then divides by the most terms from
  // the denominator that it can without error. It checks for 
  // overflow by checking that the multiplication can be reversed
  // by dividing.
  cur_k = 2;

  for(cur_n = lim1+1; cur_n <= n ; cur_n ++) {
    prev_val = val;
    val = val * cur_n;
    if(val / cur_n != prev_val) {
      return -1;
    }
    // While we can divide by the denominator
    while(cur_k <= lim2 && val % cur_k == 0) {
      val = val / cur_k;
      cur_k ++;
    }
  }

  return val;
}
/* ******************************************************************************** */


/* ******************************************************************************** */
long double factorial_long(int n) {
  /*
    Computes n!.  Returns the result as a long double.

    Factorials up to 33! are stored and recalled.  Bigger factorials
    are calculated.

    WARNING: This program does not check to see if the argument of the
    factorial is bigger than the size of the maximum long double.  It
    may therefore produce overflow errors.
  */

  if (n < 0) {
    printf("Error!  Negative factorial not defined.\n");
    exit(ERR_FACTORIAL);
  }

  // Factorials
  static const long double Fact[] = {
    1.00000000000000000000E0,
    1.00000000000000000000E0,
    2.00000000000000000000E0,
    6.00000000000000000000E0,
    2.40000000000000000000E1,
    1.20000000000000000000E2,
    7.20000000000000000000E2,
    5.04000000000000000000E3,
    4.03200000000000000000E4,
    3.62880000000000000000E5,
    3.62880000000000000000E6,
    3.99168000000000000000E7,
    4.79001600000000000000E8,
    6.22702080000000000000E9,
    8.71782912000000000000E10,
    1.30767436800000000000E12,
    2.09227898880000000000E13,
    3.55687428096000000000E14,
    6.40237370572800000000E15,
    1.21645100408832000000E17,
    2.43290200817664000000E18,
    5.10909421717094400000E19,
    1.12400072777760768000E21,
    2.58520167388849766400E22,
    6.20448401733239439360E23,
    1.55112100433309859840E25,
    4.03291461126605635584E26,
    1.0888869450418352160768E28,
    3.04888344611713860501504E29,
    8.841761993739701954543616E30,
    2.6525285981219105863630848E32,
    8.22283865417792281772556288E33,
    2.6313083693369353016721801216E35,
    8.68331761881188649551819440128E36
  };

  long double nFact = 1.0; // the result
  int c; // counter

  if (n <= 33) {
    return Fact[n];
  }
  else {
    nFact = Fact[33];
    for (c = 34; c < n; c++) {
      nFact *= c;
    }
    return nFact;
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
unsigned long GetRandSeed(unsigned long s){ 
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
/* ******************************************************************************** */


/* ****** */
int fileExists( char* file) {
  FILE *fp;
  int exists = TRUE;

  fp = fopen( file, "r");
  if( fp == NULL) exists = FALSE;
  else fclose(fp);

  return exists;
}

/* ******** */


void SetExecutionPath(int nargs, char **args) {
  strncpy(COMMAND_PATH, args[0], MAX_FILENAME_LEN);
  char * cur_brk = COMMAND_PATH;
  char * last_brk = COMMAND_PATH;
  const char * path = getenv("PATH");
  const char * dir_separators = "\\/";
  const char * path_separators = ":";
  int pathlen = strlen(path) + 1;
  int commandlen = strlen(COMMAND_PATH) + 3;
  int filefound = 0;
  char * path_copy = malloc((pathlen + 1) * sizeof(char));
  char * path_tok = NULL;
  char * curpath = malloc((pathlen + commandlen) * sizeof(char));
  strncpy(path_copy, path, pathlen);

  while(cur_brk) {
    last_brk = cur_brk;
    cur_brk = strpbrk(cur_brk + 1, dir_separators);
  }
  if (last_brk != COMMAND_PATH) {
    last_brk[1] = '\0';
  } else {
    // Using a resolved path, checking against PATH
    path_tok = strtok(path_copy, path_separators);
    while(path_tok) {
      snprintf(curpath, pathlen + commandlen, "%s/%s", 
          path_tok, args[0]);
      if (fileExists(curpath)) {
        filefound = 1;
        break;
      }
      path_tok = strtok(NULL, path_separators);
    }
    if (path_tok) {
      strcpy(COMMAND_PATH, path_tok);
    } else {
      fprintf(stderr, 
          "warning: failed to find executable path, may have"
          "difficulty finding default parameters\n");
      COMMAND_PATH[0] = '\0';
    }
  }

  free(curpath);
  free(path_copy);
}

/* ******************************************************************************** */
/*                         BEGIN NUMBER PROCESSING FUNCTIONS                        */
/* ******************************************************************************** */
