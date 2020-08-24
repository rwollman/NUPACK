/*
  complexesHeader.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 8/2006 and Justin Bois 1/2007

  Header file for use with complexes.c.
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "pfuncUtilsHeader.h"
#include "complexesStructs.h"
#include "physical_constants.h"
#include "runtime_constants.h"
#include "DNAExternals.h"

#ifndef FUNCTIONSCX_H
#define FUNCTIONSCX_H


// Functions from complexesUtils.c
int getPermutation( int, int, int*);
void resetNicks( int, int*);
void nextMultiset( int, int*, int*, int*, int); //generate a multiset
int isCyclicP(int, int*, int*); //check if a cyclic permutation
void symmetryCheck( multiset*, int, permutation*); //check if symmetry
void printPerms( FILE*, int, int, multiset*);
void printMfesToFile( const dnaStructures *ds, FILE *fp,
                       const int *nicks);
int compareMultisets( const void*, const void*);
int comparePermutations(const void *, const void *);


//Functions for permBG.c
void neg( int t, int n, int k);
void gen( int t, int n, int k);
void PrintIt( void);
void BigGen( int);
void swap( int, int, int);
void initializeMP( int, int*);
void freeMP( void);
int nextPerm( void);
void setPerm( int*);
// Count the number of sets in a list of permutations
int CountSets(permutation * perms, int nPerms, int nStrands);
// Get the maximum complex size
int GetMaxComplexSize(multiset * allSets, int totalSets);
// Fill in the set information using the permutations
int FillSets(multiset * allSets, permutation * allPermutations,
              int totalSets, int totalPerms,
              int nStrands, int * seqlength);
// Generate fixed content necklaces.
// There is a more efficient algorithm
// (see Sawada 2003) but we don't need this for the commonly used (and actually
// documented) use cases of the executable. The only time this is used is when a
// C line stands alone and we need to generate all necklaces for that complex.
int makeFCPermutations(permutation * loc, int * composition, int length, int nStrands);
// fill out circular permutations of length "length" and
// alphabet size "nStrands"
int makePermutations(permutation * loc, int length, int nStrands);

// Functions from ReadCommandLine.c
int ReadCommandLine( int, char**);
int ReadInputFileComplexes(  char *filePrefix, int *nStrands,
                             char ***seqs, int **seqlength,
                             int *maxLength, int *maxComplexSize);
void printHeader( int nStrands, char **seqs, int maxComplexSize,
                  int totalOrders, int nSets, int nNewComplexes,
                  FILE *F_cx, int nargs, char **argv, int isPairs);

// Functions from utils.c in the shared directory
int gcd(int a, int b);
double factorial(int n);
long double factorial_long(int n);
int binomial_coefficient(int n, int k);
//long GetRandSeed(long s);

#endif
#ifdef __cplusplus
}
#endif
