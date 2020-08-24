/*
  complexesStructs.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 8/2006 and Justin Bois 1/2007

  Header file containing structs for use with complexes.c.
*/
#ifndef NUPACK_COMPLEXES_STRUCTS_H__
#define NUPACK_COMPLEXES_STRUCTS_H__

#ifdef __cplusplus
extern "C" {
#endif

struct LList{
  int nSeqs;
  int *code; //permutation, with absolute strand ids
  int *strand_sums;
  int *baseCode; //pairs of data, strand id, base #
  long double pf;
  char *seq;
  //long double *Q, *Qb, *Qm, *Q1r, *Q1; //N^4
  int symmetryFactor;
  struct LList *next;
};

typedef struct LList permutation;

typedef struct{
  int *code; //membership function (# equals multiplicity)

  int nSeqs; //total sequences in this subset
  char **seqs; //seqs in the subset

  int totalLength; //sum of all included sequences lengths
  int *seqlength; //pointer to sequences lengths

  long double pf; //partition function (sum over all circular permutations)
  int nPerms; // the number of permutations of the multiset
  permutation *perms; //pointer to first permutation


  long double **avgBp; //average number of pairs for strand i, base j

  //the next three values are not currently used, although they are still calculated
  long double mfe;
  int *mfePerms;  //stores permutation ids that contain mfe structures
  int nMfePerms; //number of permutations that have an mfe

} multiset;

typedef struct {
  int quiet;
  int permsOn;
  int dopairs;
  long double T;
  int dangles;
  int parameters;
  char inputFilePrefix[1000];
  int out;
  int timeonly;
  int listonly;
  int debug;
  int echo;
  int mfe;
  long double cutoff;
  int progress;
  int onlyOneMFE;
  long double sodiumconc;
  long double magnesiumconc;
  int uselongsalt;
  int dodefect;
  int v3;
} globalArgs_t;

#ifdef __cplusplus
}
#endif

#endif /* __cplusplus */
