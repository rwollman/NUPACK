/** \file mfeUtils.h
 * mfeUtils.h is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Robert Dirks 7/2006, Justin Bois 1/2007
 * MPI added: Asif Khan 8/2009, Brian Wolfe 10/2009

 * The purpose of this program is to calculate the energy of the most
 * stable fold over all secondary structures of a given strand of
 * DNA/RNA, allowing for the simplest kinds of pseudoknots.  The energy
 * algorithm will follow the general format of Zuker and later
 * Hofacker.

 * The inclusion of pseudoknots and their corresponding energies relies
 * heavily on the ideas presented by (Rivas and Eddy 1999, J Mol Bio,
 * 285, 2053-68) although their notion of Gap matrices are not used.
 * The reason for this is that gap matrices allow a given secondary
 * structure to be obtained in multiple ways via recursions, which
 * leads to multiple possible energies per fold.  The method used in
 * this program will not allow for as general structures as Rivas and
 * Eddy, but will have unique representations of each structure.  This
 * is accomplished by explicitly creating pseudoknots in the
 * recursions.

 * 05/04/2007: Bug fix for -degenerate flag in MFE calculation
 */

#ifndef __MFEUTILS_H__
#define __MFEUTILS_H__

#include "pknots.h"
#include "min.h"
#include "init.h"
#include "backtrack.h"
#include "CalculateEnergy.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* mfeFull() returns an integer array (user allocated) describing the
   minimum free energy structure, with
   base i paired to base j <=> pairs[i] = j && pairs[j] = i;
   base i is unpaired <=> pairs[i].

   The arguments are the same as nsStarFull */
DBL_TYPE mfeFull( int inputSeq[], int seqLen, int *thepairs, int complexity,
                  int naType, int dangles, DBL_TYPE temperature,
                  DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

// mfe makes the same assumptions as nsStar
DBL_TYPE mfe( int seq[], int seqLength, int *thepairs);

/* mfeFullWithSym() is used when a strand permutation
   has a specified symmetry factor e.g. ABAB has a symmetry of 2.
   The algorithm first calculates the mfe structure assuming no symmetry (mfeFull).
   If the reported mfe has a symmetric secondary structure, then its corrected
   energy should is the reported value + kT log( symmetry).  Consequently, an
   enumeration is performed of all structures within at most kT log(symmetry) of the
   calculated mfe.  The interval may be less if a structure deviating by
   exactly one pair from the mfe has a smaller energy gap.  The true mfe, taking
   into account symmetry is then calculated and reported.
   The struct mfeStructures (allocated by the user)
   will contain the enumerated structures, with the best at the front of the list */
DBL_TYPE mfeFullWithSym( int inputSeq[], int seqLen,
              dnaStructures *mfeStructures, int complexity, int naType,
              int dangles, DBL_TYPE temperature, int symmetry, int onlyOne,
              DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//mfeFullWithSym_Subopt is similar to mfeFullWithSym, but enumerates all structures
//within range of the algorithmic mfe (by algorithmic, I mean ignoring symmetry)
DBL_TYPE mfeFullWithSym_SubOpt( int inputSeq[], int seqLen,
                               dnaStructures *mfeStructures, int complexity, int naType,
                               int dangles, DBL_TYPE temperature, int symmetry,
                                DBL_TYPE fixedSubOptRange,
                                int onlyOne, DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
                                int uselongsalt);


/* getStructure() will create a dot-parens notation (stored in structure, used allocated)
   for the structure described by pairs (see above).  Will be misleading for
   pseudoknotted structures */
void getStructure( int seqlength, const int *thepairs, char *structure);

//initialize a structure of type dnaStructures
void initMfeStructures( dnaStructures*, int);

//used by qsort to arrange enumerated structures by corrected energy.
int compareDnaStructs( const void *, const void *);

//used during validation runs to sort the output based on structure
//instead of based on floating point numbers.
int compareDnaStructsOutput(const void *, const void *);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
