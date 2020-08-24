/*
  nsStar.c

  Functions for computing n(s*) and p(s*) using partition function
  algorithms.  Written by Robert Dirks.

  Modified by Justin Bois, 13 January 2007.
  Modified to include salt correction by JSB Feb 2009.
*/

#include "nsStar.h"


/* ******************** */
DBL_TYPE nsStarPairsOrParensFull( int seqlength, int seq[], int *pairs, 
				  char *parens, int complexity, int naType, int dangles, 
				  DBL_TYPE temperature, DBL_TYPE sodiumconc,
				  DBL_TYPE magnesiumconc, int uselongsalt) {
  
  DBL_TYPE explIncorrect;
  DBL_TYPE pf;
  int *thepairs;

  pf =  pfuncFull( seq, complexity, naType, dangles, temperature, 1, sodiumconc,
		   magnesiumconc, uselongsalt);
  
  if( pairs == NULL) {
    thepairs = (int*) malloc( (seqlength+1)*sizeof( int));
    getStructureFromParens( parens, thepairs, seqlength);
  }
  else
    thepairs = pairs;

  explIncorrect = seqlength - expectedCorrectBases( thepairs, seqlength);
  
  if( pairs == NULL) 
    free( thepairs);
  
  return explIncorrect;
}



