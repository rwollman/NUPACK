#include "CalculateEnergy.h"

//***********************************************************

DBL_TYPE naEnergy( char *prefix, int seq[]) {
  // Give energy of DNA strand with all other params set to defaults
  return naEnergyFull( prefix, seq, DNA, 1, 37, 1.0, 0.0, 0);
}



DBL_TYPE naEnergyFull( char prefix[], int inputSeq[], int naType, 
		       int dangles, DBL_TYPE temperature,
		       DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {
  // Give energy with possibly symmetry set to 1 and gets secondary structure from prefix.fold

  return naEnergyFullWithSym( prefix, inputSeq, naType, dangles, temperature,
			      1, sodiumconc, magnesiumconc, uselongsalt);
}


DBL_TYPE naEnergyFullWithSym( char prefix[], int inputSeq[], int naType, 
			      int dangles, DBL_TYPE temperature, int possibleSymmetry,
			      DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {

  fold thefold;
  DBL_TYPE energy;
  char *foldFile;
  int *seq;
  int size;
  int nicks[ MAXSTRANDS];
  int seqlength, nStrands;

  seqlength = getSequenceLengthInt( inputSeq, &nStrands);
  seq = (int*) malloc( (seqlength+1)*sizeof( int) ); 
  processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);

  TEMP_K = temperature + ZERO_C_IN_KELVIN;
  DNARNACOUNT = naType;
  DANGLETYPE = dangles;
  SODIUM_CONC = sodiumconc;
  MAGNESIUM_CONC = magnesiumconc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = uselongsalt;

  size = strlen(prefix) + 6;
  foldFile = (char*) malloc( size*sizeof( char));
  strcpy( foldFile, prefix); 
  strcat( foldFile, ".fold");

  //LoadEnergies();  Energies are loaded in GetEnergy

  LoadFold( &thefold, foldFile); // get input file
  thefold.seq = seq;

  energy = GetEnergy( &thefold); //Calculates Energy  

  energy += LOG_FUNC( (DBL_TYPE) checkSymmetry( thefold.pairs, seqlength, nicks, 
                                          possibleSymmetry, nStrands))*kB*TEMP_K;
  energy += (BIMOLECULAR + SALT_CORRECTION)*(nStrands-1);


  free( foldFile);
  free( thefold.pairs);
  free( thefold.pknots);
  free( thefold.fixedBases);
  free( thefold.isNicked);
  free( seq);

  return energy;
}


/* ******************* */
DBL_TYPE naEnergyPairsOrParens( int *thepairs, char *parens, int seq[]) {
  // Returns energy of structure for DNA with all other parameters set to defaults

  return naEnergyPairsOrParensFull( thepairs, parens, seq, DNA, 1, 37, 1.0, 0.0, 0);
}

/* ******************* */
DBL_TYPE naEnergyPairsOrParensFull( int *thepairs, char *parens, 
                                   int inputSeq[], int naType,
				    int dangles, DBL_TYPE temperature,
				    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, 
				    int uselongsalt) {

  // Give energy with possibly symmetry set to 1
 
  return naEnergyPairsOrParensFullWithSym( thepairs, parens,
                                         inputSeq, naType, dangles, 
					   temperature, 1, sodiumconc, magnesiumconc,
					   uselongsalt);
}

/* ******************* */
DBL_TYPE naEnergyPairsOrParensFullWithSym( int *thepairs, char *parens, 
                                          int inputSeq[], int naType,
                                          int dangles, DBL_TYPE temperature,
                                          int possibleSymmetry,
                                          DBL_TYPE sodiumconc, 
                                          DBL_TYPE magnesiumconc,   
                                          int uselongsalt) {
                                            
  fold thefold;
  DBL_TYPE energy;
  int nStrands;
  int *seq; //without the strand breaks
  int seqlength;
  int nicks[ MAXSTRANDS];
  int i;
  
  for( i = 0; i < MAXSTRANDS; i++) { //initialize nicks array
    nicks[i] = -1;
  }
  
  seqlength = getSequenceLengthInt( inputSeq, &nStrands);
  seq = (int*) malloc( (seqlength+1)*sizeof( int) ); 
  processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);
  
  TEMP_K = temperature + ZERO_C_IN_KELVIN;
  DNARNACOUNT = naType;
  DANGLETYPE = dangles;
  SODIUM_CONC = sodiumconc;
  MAGNESIUM_CONC = magnesiumconc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = uselongsalt;
  
  MakeFold( &thefold, seqlength, seq, parens, thepairs);
  thefold.seq = seq;
  
  thefold.isNicked = (int *) calloc( seqlength, sizeof( int));
  i = 0;
  while( nicks[i] != -1) {
    thefold.isNicked[ nicks[i++]] = 1;
  }
  
  energy = GetEnergy( &thefold); //Calculates Energy  
  energy += LOG_FUNC( (DBL_TYPE) checkSymmetry( thefold.pairs, seqlength, nicks, 
                                           possibleSymmetry, nStrands))*kB*TEMP_K;
  energy += (BIMOLECULAR + SALT_CORRECTION)*(nStrands-1);
  
  free( thefold.pairs);
  free( thefold.pknots);
  free( thefold.isNicked);
  free( seq);
  
  return energy;
}



/* ***************************************************** */
void MakeFold( fold *thefold, int seqlength, int seq[], char *parens, int *thepairs) {
  
  
  int init, i; // loop indices for initializations
  
  int pairsFromParens[ MAXSEQLENGTH];
  
  
  thefold->seqlength = seqlength;
  thefold->seq = seq;
  
  pairsFromParens[0] = -5; 
  if( parens != NULL) 
    getStructureFromParens( parens, pairsFromParens,  
                           thefold->seqlength);
  
  thefold->pairs = 
    (int*) calloc( thefold->seqlength+1, sizeof(int));
  if( thefold->pairs == NULL) {
    printf("Unable to allocate fold file!\n");
    exit(1);
  }
  
  thefold->pknots = 
    (int*) calloc( thefold->seqlength+1, sizeof(int));
  if( thefold->pknots == NULL) {
    printf("Unable to allocate fold file!\n");
    exit(1);
  }
  
  for( init = 0; init <= thefold->seqlength; init++) {
    thefold->pairs[init] = -1;
    thefold->pknots[init] = -1;
  }
  
  for( init = 0; init <= thefold->seqlength - 1; init++) {
    if( parens != NULL) {
      thefold->pairs[init] = pairsFromParens[init];
    } else {
      thefold->pairs[init] = thepairs[init];
    }
  }
  
  
  //the following pknot finding routine should be optimized
  for( init = 0; init <= thefold->seqlength-1; init++) {
    if( thefold->pairs[init] > init) {
      for( i = 0; i < init; i++) {
        if( thefold->pairs[i] > init && thefold->pairs[i] < thefold->pairs[init]) {
          if( thefold->pknots[i] == -1) {
            
            thefold->pknots[i] = thefold->pairs[init];
            thefold->pknots[ thefold->pairs[init]] = i;
          }
          break;
        }
      }
      
    }
  }
  
}


/* *************** */


