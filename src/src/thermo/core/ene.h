/** \file ene.h
 * ene.h is part of the NUPACK software suite
 * Copyright (c) 2007 Caltech. All rights reserved.
 * Coded by: Robert Dirks 3/2006, Justin Bois 1/2007
 *
 * This file contains energy functions used for determining energies
 */

#ifndef ENE_H
#define ENE_H

#include <string.h>
#include <time.h>
#include <float.h>
#include <shared.h>

#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

//Nearest neighbor energy of two consecutive base pairs
DBL_TYPE HelixEnergy( int i, int j, int h, int m);

//interior mismatch energy
DBL_TYPE InteriorMM( char a, char b, char x, char y);

//hairpin energy
DBL_TYPE HairpinEnergy( int i, int j, int seq[] );

//interior loop energy
DBL_TYPE InteriorEnergy(  int i, int j, int h, int m, int seq[]);

//interior loop energy, but allows for the exclusion of the i-j mismatch
//contribution (for use with fast i loop calculations)
DBL_TYPE InteriorEnergyFull( int i, int j, int h, int m, int seq[], int);

//Calculate dangle energies, assuming the entire structure is known
//(and can hence accurately calculate wobble pair dangles)
DBL_TYPE DangleEnergyWithPairs( int i, int j, fold *thefold);

//Calculates the dangle energy of a subsequence (i,j), assuming
//i-1, j+1 are paired (unless near a strand break).
//DangleEnergy miscalculates the dangles of nearby wobble pairs
DBL_TYPE DangleEnergy( int i, int j, int seq[], int seqlength);

//Calculates exp(-(dangle energy)/RT) and
//exp( -(interior loop energy)/RT), respectively
extern unsigned int seqHash;
DBL_TYPE ExplDangle(int i, int j, int seq[], int seqlength);
DBL_TYPE ExplInternal(int i, int j, int h, int m, int seq[]);

//NickDangle calculates the dangle energy, taking into account the effects
//of strand breaks (nicks).  If hairpin == TRUE, then this region is a nicked hairpin
//and may be closed by a wobble pair
DBL_TYPE NickDangle(int i, int j, const int *nicks, int **etaN, int hairpin,
                    int seq[], int seqlength);

/* Computes the energy of an exterior loop with no secondary structure,
   and returns either exp( -energy/RT) or simply energy*/
DBL_TYPE NickedEmptyQ( int i, int j, const int nicks[], int seq[],
                       int seqlength, int **etaN);
DBL_TYPE NickedEmptyF( int i, int j, const int nicks[], int seq[],
                       int seqlength, int **etaN);

// Lookup table for the function 1.75*kB*TEMP_K*LOG_FUNC( size/30.0)
DBL_TYPE sizeLog(int size);

// Lookup table for contraction part of prFastILoops
DBL_TYPE sizeEnergyLog(int size);

//Computes the contribution of asymmetry and size to a large interior loop
DBL_TYPE asymmetryEfn( int L1, int L2, int size);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // ENE_H
