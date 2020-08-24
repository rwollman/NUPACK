#ifndef PFUNC_PF_H__
#define PFUNC_PF_H__


#include "pairsPr.h"
#include "sumexp.h"
#include "sumexp_pk.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
   pfuncFull: Calculates the partition function.
   Arguments:
   InputSeq is the sequence, including '+' characters to demarcate strand breaks.

   Complexity = 3 is the multi-stranded, pseudoknot-free method  [ O(N^3)]
   Complexity = 5 is the single-stranded, pseudoknot algorithm   [ O(N^5)]

   naType is chosen from the enum list given above, and chooses the parameter set

   dangles = 0, no dangle energies
   dangles = 1, mfold treatment of dangles
   dangles = 2, dangle energies are always summed, regardless of nearby structures
   (same as the Vienna package -D2 option)
   temperature = the temperature in celsius
   calcPairs = 1 indicates that pair probabilities are calculated and stored in global pairPr
   calcPairs = 0 skips pair probability calculations

   Ignores the possibility of symmetry
*/

DBL_TYPE pfuncFull( int inputSeq[], int complexity, int naType, int dangles,
                    DBL_TYPE temperature, int calcPairs,
                    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

DBL_TYPE pfuncFullWithBonuses(int inputSeq[], int complexity, int naType,
    int dangles, DBL_TYPE temperature, int calc_pairs, int perm_sym, 
    DBL_TYPE sodium_conc, DBL_TYPE magnesium_conc, int use_long_salt, 
    DBL_TYPE * bonuses);

//pfuncFullWithSym is the Same as pfuncFull, but divides
//the result by permSym to account for symmetries
DBL_TYPE pfuncFullWithSym( int inputSeq[], int complexity, int naType,
                           int dangles, DBL_TYPE temperature, int calcPairs, int permSymmetry,
                           DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

DBL_TYPE pfuncFullWithSymHelper( int inputSeq[], int seqlength, int nStrands,
                                 int complexity, int naType,
                                 int dangles, DBL_TYPE temperature, int calcPairs,
                                 int permSymmetry, DBL_TYPE sodiumconc,
                                 DBL_TYPE magnesiumconc, int uselongsalt);


/* pfunc
   Calls pfuncFull, and assumes complexity = 3, DNA parameters, T = 37, dangles = 1,
   calcPairs = 1, [Na+] = 1.0, [Mg++] = 0.0, and short helix model for salt correction
*/
DBL_TYPE pfunc( int seq[]);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
