/** \file CalculateEnergy.h
 *CalculateEnergy.h is part of the NUPACK software suite
 *Copyright (c) 2007 Caltech. All rights reserved.
 *Coded by: Robert Dirks 4/2004, Justin Bois 1/2007
 */

#ifndef __CALCULATEENERGY_H__
#define __CALCULATEENERGY_H__

#include "shared.h"
#include "GetEnergy.h"
#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

/* naEnergyPairsOrParensFull() Computes the energy of a structure
   given either an int array of pairs (thepairs), or a character array
   (parens) in dot-parens notation.  See nsStarPairsOrParensFull().
   Assumes no symmetry.  Include the bimolecular association penalty
*/
DBL_TYPE naEnergyPairsOrParensFull( int *thepairs, char *parens,
                                   int inputSeq[], int naType,
                                    int dangles, DBL_TYPE temperature,
                                    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
                                    int uselongsalt);

//same as above, but allows a value for possible symmetry (v(pi) = symmetry of sequence permutation)
DBL_TYPE naEnergyPairsOrParensFullWithSym( int *thepairs, char *parens,
                                          int inputSeq[], int naType,
                                          int dangles, DBL_TYPE temperature,
                                           int possibleSymmetry,
                                           DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
                                           int uselongsalt);

//naEnergyFull() is the same as above, but obtains the secondary structure from prefix.fold
DBL_TYPE naEnergyFull( char prefix[], int inputSeq[], int naType,
                       int dangles, DBL_TYPE temperature,
                       DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//Same as naEnergyFull(), but allows possible symmetry (see above)
DBL_TYPE naEnergyFullWithSym( char prefix[], int inputSeq[], int naType,
                              int dangles, DBL_TYPE temperature, int possibleSymmetry,
                              DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//computes energy of secondary structure in file, using same assumptions as pfunc
DBL_TYPE naEnergy( char *prefix, int seq[]);
//same as above, but structure input in thepairs or parens
DBL_TYPE naEnergyPairsOrParens( int *thepairs, char *parens,
                                int inputSeq[]);

//Makefold() creates a structure of type fold from parens or thepairs.
void MakeFold( fold *thefold, int seqlength, int seq[], char *parens, int *thepairs);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
