/** \file init.h
 * init.h is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Robert Dirks 3/2006, Justin Bois 1/2007
 *
 * Functions to be run once at the beginning of the
 * partition function algorithm
 */

#ifndef PFUNC_INIT_H__
#define PFUNC_INIT_H__

#include <shared.h>
#include "ene.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ReadSequence() reads in a nucleic acid sequence from a file
   and allocates *seq to store it */
void ReadSequence(int *seqlength, char **seq, char *filename );

/* getSequenceLength() counts the number of strands (stored as nStrands)  in a given sequence (seq)
   and returns the sequence length */
int getSequenceLength(char *seq, int *nStrands);

/* getSequenceLength() counts the number of strands (stored as nStrands)  in a given sequence (seq)
   and returns the sequence length */
int getSequenceLengthInt(int seq[], int *nStrands);

/* processMultiSeqence() copies input sequence (containing strandPLUS members) to seq, but without the
   strandPLUS members.  The location of the strand breaks, strandPLUS, are stored in the nicks array */
void processMultiSequence(int inputSeq[], int seqlength, int nStrands,
                           int seq[], int nicks[]);

//Allocates Q and sets the values to zero.
void InitLDoublesMatrix(DBL_TYPE **Q, int size, char name[]);

//Sets Q to all zero
void ClearLDoublesMatrix(DBL_TYPE **Q, int size, char name[]);

//Memory management for "fast" interior loops subroutine
void manageQx(DBL_TYPE **Qx, DBL_TYPE **Qx_1,
               DBL_TYPE **Qx_2, int len, int seqlength);
//manageFx is the mfe version
void manageFx(DBL_TYPE **Fx, DBL_TYPE **Fx_1,
               DBL_TYPE **Fx_2, int len, int seqlength);


//Memory management routines for QgIx ("fast" gap spanning array treatment)
void manageQgIx(DBL_TYPE **QgIx, DBL_TYPE **QgIx_1,
                 DBL_TYPE **QgIx_2, int d, int seqlength);
void manageFgIx(DBL_TYPE **FgIx, DBL_TYPE **FgIx_1,
                 DBL_TYPE **FgIx_2, int d, int seqlength);

// Salt correction
DBL_TYPE computeSaltCorrection(DBL_TYPE sodiumConc, DBL_TYPE magenesiumConc,
                               int useLongHelix);

//Load energy parameters.  Global variable DNARNACOUNT determines parameter set
void LoadEnergies(void);
void setParametersToZero(void);

//Set Q[ pf_index(i, i-1, seqlength)] = 1;
void nonZeroInit(DBL_TYPE Q[], int seq[], int seqlength);

/*InitEtaN() initializes the etaN array.  etaN[ EtaNIndex(i,j,seqlength)][0] is the
  number of nicks between i and  j (i,j in [0.5, 1.5, 2.5, ...]).
  etaN[ EtaNIndex(i,j,seqlength)][1] is the index of the leftmost nick
  in the region between i and j, i.e. nicks[ EtaNIndex...] is the position
  of the leftmost nick between i and j.
*/
void InitEtaN(int **etaN, const int *nicks, int seqlength);
int EtaNIndex_old(float i, float j, int seqlength);

/* These functions set the size of Qg and Fg matrices is correct
   between successive calls to pfunc or mfe.  They also call PrecomputeValuesN5(f)  */
void initPF(int seqlength);
void initMfe(int seqlength);

//precomputes the sizeTerm array, used in extending interior loops from
//size i to i+2
void PrecomputeValuesN5(int seqlength);
//PrecomputeValuesN5f is for the mfe algorithm
void PrecomputeValuesN5f(int seqlength);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
