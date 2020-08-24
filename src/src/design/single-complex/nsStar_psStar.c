/*
  nsStar_psStar.c

  Functions for computing n(s*) and p(s*) using partition function
  algorithms.  Written by Robert Dirks.

  Modified by Justin Bois, 13 January 2007.
  Modified to include salt correction by JSB Feb 2009.
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "shared/structs.h" 
#include <thermo/core/pfuncUtils.h>
#include <thermo/core/pf.h>
#include "design_pfunc_utils.h"

/* ******************** */
DBL_TYPE nsStar( char prefix[], int seq[]) {
  // n(s*) with everyting set to defaults
  return nsStarFull( prefix, seq, 3, 0, 1, 37.0, 1.0, 0.0, 0);
}


/* ******************** */
DBL_TYPE nsStarFull( char prefix[], int seq[], int complexity, int naType, 
                     int dangles, DBL_TYPE temperature, DBL_TYPE sodiumconc,
		     DBL_TYPE magnesiumconc, int uselongsalt) {
  
  int size;
  char *foldFile;
  
  DBL_TYPE explIncorrect;
  //int seqlength;
  DBL_TYPE pf;
  
  /* DNA specific stuff */
  fold thefold;
  
  size = strlen(prefix) + 6;
  foldFile = (char*) malloc( size*sizeof( char));
  strcpy( foldFile, prefix); 
  strcat( foldFile, ".fold");
  
  LoadFold( &thefold, foldFile); 
  pf =  pfuncFull( seq, complexity, naType, dangles, temperature, 1 , sodiumconc,
  		   magnesiumconc, uselongsalt);
  
  explIncorrect = thefold.seqlength - 
                  expectedCorrectBases( thefold.pairs, thefold.seqlength);
  
  free( foldFile);
  free( thefold.pairs);
  free( thefold.pknots);
  free( thefold.fixedBases);
  free( thefold.isNicked);
  
  return explIncorrect;
}


/* ******************** */
DBL_TYPE nsStarPairsOrParens( int seqlength, int seq[], int *pairs, char *parens) {

  // Everything is set to defaults

  return nsStarPairsOrParensFull( seqlength, seq, pairs, parens, 3, DNA, 1, 37.0,
				  1.0, 0.0, 0);
}


/* ******************** */
DBL_TYPE nsStarPairsOrParensFull( int seqlength, int seq[], int *pairs, 
				  char *parens, int complexity, int naType, int dangles, 
				  DBL_TYPE temperature, DBL_TYPE sodiumconc,
				  DBL_TYPE magnesiumconc, int uselongsalt) {
  
  DBL_TYPE explIncorrect;
  DBL_TYPE pf;
  int *thepairs;

  /* DNA specific stuff */
  pf =  pfuncFull( seq, complexity, naType, dangles, temperature, 1, sodiumconc,
		   magnesiumconc, uselongsalt);
  //printf("%Le %d\n", pf, seqlength);
  
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



/* ******************** */

DBL_TYPE expectedCorrectBasesWithDummy( int *structPairs, int *dummyBases, 
                                       int seqlength) {
 
  int i;
  DBL_TYPE value = 0;
  int pair;
  extern DBL_TYPE *pairPr;
  int dummyCount = 0;

  value = 0;
  for( i = 0; i< seqlength; i++) {
    if (!dummyBases[i]) {
      pair = structPairs[i];
      if( pair == -1) {
        value += pairPr[ i*(seqlength+1) + seqlength]; //unpaired probability
      }
      else {
        value += pairPr[ i*(seqlength+1) + pair];
      }
    }
    else {
      dummyCount++;
    }
  }

return value + (DBL_TYPE)dummyCount;
}


/* ******************** */
DBL_TYPE nsStarPairsOrParensWithSym( int seqlength, int nStrands,int seq[], 
          int pairs[], char *parens, int *dummyBases, DBL_TYPE *pfVal,
          int complexity, int naType, int dangles, 
	  DBL_TYPE temperature, int symmetry, DBL_TYPE sodiumconc,
	  DBL_TYPE magnesiumconc, int uselongsalt) {
  
  DBL_TYPE explIncorrect;
  DBL_TYPE pf;
  int *thepairs;
  
  /* DNA specific stuff */
  pf =  pfuncFullWithSymHelper( seq, seqlength, nStrands, complexity, naType, 
				dangles, temperature, 1, symmetry, sodiumconc, 
				magnesiumconc, uselongsalt);


  *pfVal = pf;
  
  if( pairs == NULL) {
    thepairs = (int*) malloc( (seqlength+1)*sizeof( int));
    getStructureFromParens( parens, thepairs, seqlength);
  }
  else
    thepairs = pairs;

  explIncorrect = seqlength - 
              expectedCorrectBasesWithDummy( thepairs, dummyBases, seqlength);

  if( pairs == NULL) 
    free( thepairs);

  return explIncorrect;
}

/* ******************** */
DBL_TYPE nsStarPairsOrParensCorrected( int seqlength, int nStrands,int seq[], 
          int pairs[], char *parens, int *dummyBases, DBL_TYPE *pfVal,
          int complexity, int naType, int dangles, 
	  DBL_TYPE temperature, DBL_TYPE sodiumconc,
	  DBL_TYPE magnesiumconc, int uselongsalt) {
  
  DBL_TYPE explIncorrect;
  DBL_TYPE pf;
  int *thepairs;
  
  /* DNA specific stuff */
  pf =  pfuncFullWithSymHelper( seq, seqlength, nStrands, complexity, naType, 
				dangles, temperature, 1, 1, sodiumconc, 
				magnesiumconc, uselongsalt);


  *pfVal = pf;
  
  if( pairs == NULL) {
    thepairs = (int*) malloc( (seqlength+1)*sizeof( int));
    getStructureFromParens( parens, thepairs, seqlength);
  }
  else
    thepairs = pairs;

  explIncorrect = seqlength - 
              expectedCorrectBasesWithDummy( thepairs, dummyBases, seqlength);

  if( pairs == NULL) 
    free( thepairs);

  return explIncorrect;
}


/* ********** */
DBL_TYPE psStar( char prefix[], int seq[]) {
  return psStarFull( prefix, seq, 3, 0, 1, 37.0, 1.0, 0.0, 0);
}


/* ********** */
DBL_TYPE psStarFull( char prefix[], int seq[], int complexity, int naType,
		     int dangles, DBL_TYPE temperature, DBL_TYPE sodiumconc,
		     DBL_TYPE magnesiumconc, int uselongsalt) {
  
  return psStarFullWithSym( prefix, seq, complexity, naType,
                            dangles, temperature, 1, sodiumconc, magnesiumconc,
			    uselongsalt);
}


/* ********** */
DBL_TYPE psStarFullWithSym( char prefix[], int seq[], int complexity, int naType,
          int dangles, DBL_TYPE temperature, int possibleSymmetry, DBL_TYPE sodiumconc,
	  DBL_TYPE magnesiumconc, int uselongsalt) {

  DBL_TYPE pf, value;
  DBL_TYPE energy;

  int nStrands;

  /* DNA specific stuff */
  getSequenceLengthInt( seq, &nStrands); //used solely to set nStrands
  pf =  pfuncFullWithSym( seq, complexity, naType, dangles, temperature, 0,
                          possibleSymmetry, sodiumconc, magnesiumconc, uselongsalt);

  energy = naEnergyFullWithSym( prefix, seq, naType, dangles, temperature, 
                                possibleSymmetry, sodiumconc, magnesiumconc,
				uselongsalt);

  value = EXP_FUNC( -1.0*energy/(kB*TEMP_K) );
  value /= pf;

  return value;
}


/* ********** */
DBL_TYPE psStarPairsOrParens( int *pairs, char *parens, int seq[]) {
  return psStarPairsOrParensFull( pairs, parens, seq, 3, DNA, 1, 37.0, 1.0, 0.0, 0);
}

/* ********** */
DBL_TYPE psStarPairsOrParensFull( int *pairs, char *parens,  int seq[], 
          int complexity, int naType, int dangles, DBL_TYPE temperature,
          DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {
  return psStarPairsOrParensFullWithSym( pairs, parens, seq, complexity, naType,
					 dangles, temperature, 1, sodiumconc,
					 magnesiumconc, uselongsalt);
}

/* ********** */
DBL_TYPE psStarPairsOrParensFullWithSym( int *pairs, char *parens,  int seq[], 
           int complexity, int naType, int dangles, 
           DBL_TYPE temperature, int possibleSymmetry, DBL_TYPE sodiumconc,
 	   DBL_TYPE magnesiumconc, int uselongsalt) {
  
  DBL_TYPE pf, value;
  DBL_TYPE energy;
  int nStrands;
  int tmpLength;
  /* DNA specific stuff */
  
  tmpLength=getSequenceLengthInt( seq, &nStrands); //used solely to set nStrands

  pf =  pfuncFullWithSym( seq, complexity, naType, dangles, temperature, 0,
                          possibleSymmetry, sodiumconc, magnesiumconc, uselongsalt);
  
  energy = naEnergyPairsOrParensFullWithSym( pairs, parens, seq, naType, dangles, 
					     temperature, possibleSymmetry, sodiumconc,
					     magnesiumconc, uselongsalt);
  
  value = EXP_FUNC( -1.0*energy/(kB*TEMP_K) );
  value /= pf;
  
  return (DBL_TYPE) value;
}


/* ************* */
DBL_TYPE nsStarPairsOrParens_ms( int seqlength, int seq[], int *pairs, 
                                char *parens, char *ps) {

  return nsStarPairsOrParensFull_ms( seqlength, seq, pairs, parens, 3, DNA, 
				     1, 37.0, 1.0, 0.0, 0, ps);
}


/* ************* */
DBL_TYPE nsStarPairsOrParensFull_ms( int seqlength, int seq[], int *pairs, 
                                     char *parens, int complexity,
                                     int naType, int dangles, 
                                     DBL_TYPE temperature, DBL_TYPE sodiumconc,
				     DBL_TYPE magnesiumconc, int uselongsalt, 
				     char *ps) {

    DBL_TYPE expIncorrect = 
    nsStarPairsOrParensFull( seqlength, seq, pairs, parens, 
                             complexity, naType, dangles,
                             temperature, sodiumconc, magnesiumconc, uselongsalt);
  
  strncpy( ps, parens, seqlength);
  makePairStruct( ps, pairPr, seqlength);
  
  return expIncorrect;
  
}

