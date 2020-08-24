/* 
   design_pfunc_utils.h by Robert Dirks 07/24/2001
   updated: 2006/07/16
   updated: 2007/01/13 by Justin Bois
   
   Note that DBL_TYPE is defined in pfuncUtilsConstants.h, and
   indicates what type of floating point variables should be used
   (float, double, long double)
   
   See below for descriptions of each function
*/


#ifndef DESIGN_PFUNC_UTILS_H
#define DESIGN_PFUNC_UTILS_H

#include <shared.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/* ******************************************************************************** */
// The following functions are in pairPrStruct.c
/* ******************************************************************************** */
//The pairInfo structure is used to calculate a secondary structure based on
//the pair probability matrix.
typedef struct {
  int i;
  int j;
  double value;
  //base i paired to base j with probability to value
  //0 <= i <= seqlength -1; 0 <= j <= seqlength,
  //with j == seqlength indicating base i is unpaired
} pairInfo;

//used to sort pairInfo structs based on value
int compare_PI( const void*, const void*); //compare two pairInfo objects

//use a greedy algorithm to construct a structure
void makePairStruct( char*, DBL_TYPE*, int); 
/* ******************************************************************************** */


/* ******************************************************************************** */
//  The following functions are for computing design objective functions p(s*) and
//  n(s*).  See Dirks, et al, Nucleic Acids Res., 32. 1392, 2004.
/* ******************************************************************************** */

/* nsStarFull
   calculates the value of n(s*) for the given sequence (seq), and the structures stored in 
   prefix.fold.  Complexity, dangles, Temperature are the same as for pfuncFull
*/
DBL_TYPE nsStarFull( char prefix[], int seq[], int complexity, int naType, 
                     int dangles, DBL_TYPE temperature, DBL_TYPE sodiumconc,
		     DBL_TYPE magnesiumconc, int uselongsalt);

// nsStar computes n(s*) with the same assumptions as pfunc
DBL_TYPE nsStar( char prefix[], int seq[]);

/* nsStarPairsOrParens is the same as nsStar,
   except for the way in which the secondary structure is denoted.  If pairs != NULL, then
   pairs is an int array with base i paired to base j <=> pairs[i] = j && pairs[j] = i.
   Base i is unpaired <=> pairs[i];  Otherwise, parens is a character array that gives the
   structure in dot-parens notation, where dots indicate unpaired bases, and (), {}, [], or <>
   indicate base pairs.  Pairs of the same type are assumed to be nested with respect to one another.
   For these functions to behave properly, either pairs == NULL or parens == NULL.
*/
DBL_TYPE nsStarPairsOrParens( int seqlength, int seq[], int *pairs, 
                              char *parens);

DBL_TYPE nsStarPairsOrParensFull( int seqlength, int seq[], int *pairs, 
				  char *parens, int complexity, int naType, int dangles, 
				  DBL_TYPE temperature, DBL_TYPE sodiumconc,
				 DBL_TYPE magnesiumconc, int uselongsalt);


//psStar functions are identical to nsStar functions, except p(s*) is calculated rather than n(s*)
DBL_TYPE psStar( char *prefix, int seq[]);

DBL_TYPE psStarFull( char prefix[], int seq[], int complexity, int naType,
		      int dangles, DBL_TYPE temperature, DBL_TYPE sodiumconc,
		      DBL_TYPE magnesiumconc, int uselongsalt);

//The WithSym versions allows for v(pi) to be entered (i.e. symmetry of sequence permutation)
DBL_TYPE psStarFullWithSym( char prefix[], int seq[], int complexity, int naType,
          int dangles, DBL_TYPE temperature, int possibleSymmetry, DBL_TYPE sodiumconc,
	  DBL_TYPE magnesiumconc, int uselongsalt);

DBL_TYPE psStarPairsOrParens(int *pairs, char *parens, int seq[] );

DBL_TYPE psStarPairsOrParensFull( int *pairs, char *parens,  int seq[], 
          int complexity, int naType, int dangles, DBL_TYPE temperature,
	  DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

DBL_TYPE psStarPairsOrParensFullWithSym( int *pairs, char *parens,  int seq[], 
           int complexity, int naType, int dangles, 
           DBL_TYPE temperature, int possibleSymmetry, DBL_TYPE sodiumconc,
	   DBL_TYPE magnesiumconc, int uselongsalt);

/* nsStar*_ms are identical to the nsStar functions, except they also return a character array
   array ps (which the user needs to allocate), that contains a consensus secondary structure 
   (in dot-parens notation) using a greedy algorithm that chooses the most likely base-pairs
   and unpaired bases first */

DBL_TYPE nsStarPairsOrParens_ms( int seqlength, int seq[], int *pairs, 
                                char *parens, char *ps);

DBL_TYPE nsStarPairsOrParensFull_ms( int seqlength, int seq[], int *pairs, 
                                     char *parens, int complexity,
                                     int naType, int dangles, 
                                     DBL_TYPE temperature, DBL_TYPE sodiumconc,
				     DBL_TYPE magnesiumconc, int uselongsalt, 
				     char *ps);

/* calculates nsStar corrected for use in heirarchical deocmposition */
DBL_TYPE nsStarPairsOrParensCorrected( int seqlength, int nStrands,int seq[], 
          int pairs[], char *parens, int *dummyBases, DBL_TYPE *pfVal,
          int complexity, int naType, int dangles, 
	  DBL_TYPE temperature, DBL_TYPE sodiumconc,
	  DBL_TYPE magnesiumconc, int uselongsalt);

DBL_TYPE nsStarPairsOrParensWithSym( int seqlength , int nStrands, int seq[],
        int pairs[], char * parens, int * dummyBases, DBL_TYPE *pfVal,
        int complexity, int naType, int dangles, DBL_TYPE temperature,
        int symmetry, DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
        int uselongsalt);

DBL_TYPE expectedCorrectBasesWithDummy( int *structPairs, int *dummyBases, 
                                       int seqlength);

#ifdef __cplusplus
}
#endif 

#endif
