#include "mfeUtils.h"

int containsPk;

/* ****************************** */
// This is the main MFE calculator.  Actually finds all suboptimal folds
// with energy below fixedSubOptRange, which if < 0, does MFE
DBL_TYPE mfeFullWithSym_SubOpt( int inputSeq[], int seqLen, 
                               dnaStructures *mfeStructures, int complexity, int naType, 
                               int dangles, DBL_TYPE temperature, int symmetry, DBL_TYPE fixedSubOptRange,
                               int onlyOne, DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, 
                               int uselongsalt) {
                                 
  //if fixedSubOptRange > 0, then enumerate all structures with fsor, rather than find mfe

  DBL_TYPE result;
  int seqlength;

  DBL_TYPE *F = NULL;
  DBL_TYPE *Fb = NULL;
  DBL_TYPE *Fm = NULL;

  //N^3
  DBL_TYPE *Fx = NULL;
  DBL_TYPE *Fx_1 = NULL;
  DBL_TYPE *Fx_2 = NULL;
  DBL_TYPE *Fs = NULL;
  DBL_TYPE *Fms = NULL;

  //PKNOTS
  DBL_TYPE *Fp = NULL;
  DBL_TYPE *Fz = NULL;  //O(N^2)
  DBL_TYPE *Fg = NULL; //O(N^4)

  //N^5
  DBL_TYPE *FgIx = NULL;
  DBL_TYPE *FgIx_1 = NULL;
  DBL_TYPE *FgIx_2 = NULL;
  DBL_TYPE *Fgls = NULL;
  DBL_TYPE *Fgrs = NULL;

  DBL_TYPE *Fgl = NULL;
  DBL_TYPE *Fgr = NULL; //O(N^4) space


  /*  F-type matrices are dynamically allocated matrices that
  contain minimum energies restricted to a subsequence of the
  strand.  Each of the above should be accessed by the call
  F[ pf_index(i, j, seqlength)] to indicate the partition function between
  i and j, inclusive. 

  Descriptions of each are in the referenced paper (see pfunc.c)
  */

  int i, j, k; // the beginning and end bases for F
  long int maxIndex;
  int L; //This the length of the current subsequence
  DBL_TYPE min_energy;
  int pf_ij;
  DBL_TYPE tempMin;

  extern long int maxGapIndex;
  short *possiblePairs;

  int nicks[ MAXSTRANDS];  //the entries must be strictly increasing
  //nicks[i] = N means a strand ends with base N, and a new one starts at N+1

  int **etaN;
  int arraySize;
  int nStrands;
  int *seq;
  int *foldparens;

  DBL_TYPE mfeEpsilon;
  DBL_TYPE *minILoopEnergyBySize;
  int *maxILoopSize;

  DBL_TYPE localEnergy;
  int oldp1;
  int symmetryOfStruct = 1; // Must be initialized
  int *thepairs;  


  //assign global variables
  TEMP_K = temperature + ZERO_C_IN_KELVIN;
  DNARNACOUNT = naType;
  DANGLETYPE = dangles;
  SODIUM_CONC = sodiumconc;
  MAGNESIUM_CONC = magnesiumconc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = uselongsalt;

  // Get dG_salt. It will be used to calculate 
  DBL_TYPE salt_correction = computeSaltCorrection(sodiumconc,magnesiumconc,uselongsalt);

  seqlength = getSequenceLengthInt( inputSeq, &nStrands);




  mfeEpsilon = kB*TEMP_K*LOG_FUNC(symmetry); //max range to search when looking for mfe

  if( fixedSubOptRange > 0) {
   mfeEpsilon = fixedSubOptRange + mfeEpsilon;
  }

  for( i = 0; i < MAXSTRANDS; i++) { //initialize nicks array
   nicks[i] = -1;
  }

  seq = (int *) malloc( (seqLen+1)*sizeof( int) );
  processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);
  foldparens = (int*) malloc( (seqLen+nStrands)*sizeof(int));

  if( nStrands >= 2 && complexity >= 5) {
   printf("Warning, pseudoknots not allowed for multi-stranded complexes!");
   printf("  Pseudoknots disabled.\n");
   complexity = 3;
  }
  LoadEnergies();

  if( complexity >= 5) //pseudoknotted
   initMfe( seqlength);


  arraySize = seqlength*(seqlength+1)/2 + (seqlength+1);
  // Allocate and Initialize Matrices
  InitLDoublesMatrix( &F, arraySize, "F");
  InitLDoublesMatrix( &Fb, arraySize, "Fb");
  InitLDoublesMatrix( &Fm, arraySize, "Fm");

  etaN = (int**) malloc( arraySize*sizeof( int*));
  InitEtaN( etaN, nicks, seqlength);

  maxILoopSize = (int*) malloc( arraySize*sizeof( int));
  minILoopEnergyBySize = (DBL_TYPE*) malloc( seqlength*sizeof( DBL_TYPE));

  if( complexity == 3) {
   InitLDoublesMatrix( &Fs, arraySize, "Fs");
   InitLDoublesMatrix( &Fms, arraySize, "Fms");
  }

  if( complexity >= 5) {
   InitLDoublesMatrix( &Fp, arraySize, "Fp");
   InitLDoublesMatrix( &Fz, arraySize, "Fz");
   InitLDoublesMatrix( &Fg, maxGapIndex, "Fg");
   
   if( complexity == 5) {
     InitLDoublesMatrix( &Fgl, maxGapIndex, "Fgl");
     InitLDoublesMatrix( &Fgr, maxGapIndex, "Fgr");
     InitLDoublesMatrix( &Fgls, maxGapIndex, "Fgls");
     InitLDoublesMatrix( &Fgrs, maxGapIndex, "Fgrs");
     CheckPossiblePairs( &possiblePairs, seqlength, seq);
   }
  }

  //Initialization to NAD_INFINITY

  if( complexity >= 5) 
   maxIndex = maxGapIndex; //beware overflow
  else
   maxIndex = arraySize;

  for( i = 0; i < maxIndex; i++) {    
   if(  i < arraySize ) {
     F[i] = Fb[i] = Fm[i] = NAD_INFINITY; 
     if( complexity == 3)
       Fs[i] = Fms[i] = NAD_INFINITY;
     
     if( complexity >= 5) 
       Fp[i] = Fz[i] = NAD_INFINITY; 
   } 
   if( complexity >= 5) {
     Fg[i] = NAD_INFINITY;
     if( complexity == 5)
       Fgl[i] = Fgr[i] = Fgls[i] = Fgrs[i] = NAD_INFINITY;
   }
  }

  for( i = 0; i <= seqlength; i++) {
   pf_ij = pf_index( i, i-1, seqlength);
   F[ pf_ij] = NickDangle(i, i-1, nicks, etaN, FALSE, seq, seqlength);
   
   if( complexity >= 5)
     Fz[ pf_ij] = F[ pf_ij];
  }


  for( L = 1; L <= seqlength; L++) {
   /* Calculate all sub energies for
    length = 0, then 1, then 2.... */
    int iMin = 0;
    int iMax = seqlength - L;

   
    if( complexity == 3) 
      manageFx( &Fx, &Fx_1, &Fx_2, L-1, seqlength);   
   //allocate/deallocate memory
   
    if( complexity == 5) 
      manageFgIx( &FgIx, &FgIx_1, &FgIx_2, L-1, seqlength);
   //manageQgIx manages the temporary matrices needed for 
   //calculating Qg_closed in time n^5
   
    for( i = iMin; i <= iMax; i++) {
      j = i + L - 1;
      pf_ij = pf_index( i, j, seqlength);
     
     //store the maximum iloop size with mfeEpsilon of mfe
     for( k = 0; k < L; k++) minILoopEnergyBySize[k] = NAD_INFINITY; //initialize to zero;
     
     /* Recursions for Fb */
     /* bp = base pairs, pk = pseudoknots */
     
     min_energy = NAD_INFINITY;
     if( CanPair( seq[ i], seq[ j]) == FALSE) {
       Fb[ pf_ij] = NAD_INFINITY;
     }
     else {
       min_energy = MinHairpin( i, j, seq, seqlength, etaN);
       // Exactly 1 bp
       
       if( complexity == 3) {
         if( etaN[ EtaNIndex(i+0.5, i+0.5, seqlength)][0] == 0 &&
            etaN[ EtaNIndex(j-0.5, j-0.5, seqlength)][0] == 0) {
              //regular multiloop.  No top-level nicks
              
              tempMin = MinMultiloops(i, j, seq, Fms, Fm,
                                      seqlength, etaN);
              min_energy = MIN( tempMin, min_energy);
            }
         
         if( etaN[ EtaNIndex(i+0.5, j-0.5, seqlength)][0] >= 1) {
           //Exterior loop (created by nick)
           tempMin = MinExteriorLoop( i, j, seq, seqlength, 
                                     F, nicks, etaN); 
           min_energy = MIN( tempMin, min_energy);
         }
       }
       
       if( complexity != 3) {
         // Interior Loop and Multiloop Case
         tempMin = MinInterior_Multi( i, j, seq, seqlength, Fm, Fb, nicks, etaN);
         min_energy = MIN( tempMin, min_energy);
       }
       
       if( complexity >= 5) {
         
         tempMin =  MinFb_Pk( i, j, seq, seqlength, Fp, Fm );
         min_energy = MIN( tempMin, min_energy);
       }
       
     }
     
     if( complexity == 3) 
       MinFastILoops( i, j, L, seqlength, seq, etaN, Fb, Fx, Fx_2, minILoopEnergyBySize);
     
     Fb[pf_ij] = MIN( Fb[ pf_ij], min_energy);

     maxILoopSize[ pf_ij] = 0;
     if( CanPair( seq[i], seq[j]) == TRUE) { 
       
       for( k = 0; k < L; k++) {
        if( minILoopEnergyBySize[k] < Fb[ pf_ij] + mfeEpsilon + ENERGY_TOLERANCE ) {
         maxILoopSize[ pf_ij] = k;
        }
       }
     }
     
     
     // Recursions for Fg, Fgls, Fgrs, Fgl, Fgr
     if( complexity == 5) {
       MakeFg_N5(i, j, seq, seqlength, Fg, Fm, Fgls, Fgrs, FgIx, FgIx_2,
                 possiblePairs);
       
       MakeFgls( i, j, seq, seqlength, Fg, Fm, Fgls);
       MakeFgrs( i, j, seq, seqlength, Fg, Fm, Fgrs);
       MakeFgl(i, j, seq, seqlength, Fg, Fgl, Fz);
       MakeFgr(i, j, seq, seqlength, Fgr, Fgl, Fz);
       
       Fp[ pf_ij] = MinFp_N5( i, j, seq, seqlength, Fgl, Fgr, Fg, Fz);
     }
     else if( complexity == 8) {
       //MakeFg_N8( i, j, seq, seqlength, Fg, Fm);
       //Fp[ pf_ij] = MinFp_N8( i, j, seq, seqlength, Fg, Fz);
     }
     
     
     if( complexity == 3) {
       /* Recursions for Fms, Fs */
       MakeFs_Fms( i, j, seq, seqlength, Fs, Fms, Fb, nicks, etaN);
       
       /* Recursions for Q, Qm, Qz */
       MakeF_Fm_N3( i, j, seq, seqlength, F, Fs, Fms, Fm,
                   nicks,etaN);  
     }
     
#ifdef test
     if( complexity == 4) 
        MakeF_Fm_N4( i, j, seq, seqlength, F, Fm, Fb);
#endif
      if( complexity >= 5)
        MakeF_Fm_Fz(i, j, seq, seqlength, F, Fm, Fz, Fb, Fp);
     
    }
  }
    result = F[ pf_index(0,seqlength-1,seqlength)];  
    if( result < NAD_INFINITY/2.0) {

     initMfeStructures( mfeStructures, seqlength);
     if( complexity == 3) {
       
       if( fixedSubOptRange <= 0) {
         bktrF_Fm_N3( 0, seqlength - 1, seq, seqlength, F, Fb, Fm, Fs, Fms,
                     nicks, etaN, mfeStructures, "F", maxILoopSize, 0, onlyOne && !NUPACK_VALIDATE);
         thepairs = mfeStructures->validStructs[0].theStruct;
         
         symmetryOfStruct = checkSymmetry( thepairs, seqlength, nicks, symmetry,
                                          nStrands);
         
         // THIS IS WHERE WE KNOW WHETHER OR NOT WE HAVE TO DO THE ENUMERATION
         
         mfeEpsilon = kB*TEMP_K*LOG_FUNC( (DBL_TYPE) symmetryOfStruct);
         //default search space is within RT log( sym) of the mfe
         
         if( mfeEpsilon > ENERGY_TOLERANCE) {
           for( i = 0; i < seqlength; i++) { 
             //check structures that differ by one base pair before doing full enumeration
             oldp1 = thepairs[i];
             if( oldp1 >= 0) {
               thepairs[i] = -1;
               thepairs[ oldp1] = -1;
               
               //no symmetry is possible if the original structure was symmetric
               localEnergy = naEnergyPairsOrParensFull( thepairs, NULL, inputSeq, naType,
                    dangles, temperature, SODIUM_CONC,
                    MAGNESIUM_CONC, 
                    USE_LONG_HELIX_FOR_SALT_CORRECTION) - 
                 ( BIMOLECULAR + SALT_CORRECTION ) *(nStrands-1); //for comparison purposes, remove bimolecular term
               
               mfeEpsilon = MIN( mfeEpsilon, localEnergy - result);

               thepairs[i] = oldp1;
               thepairs[oldp1] = i;
             }
           }
         }
       }

       //find all structures within mfeEpsilon of the mfe
       if (fixedSubOptRange > 0 || symmetryOfStruct > 1) {
         clearDnaStructures( mfeStructures);
         initMfeStructures( mfeStructures, seqlength);
         bktrF_Fm_N3( 0, seqlength - 1, seq, seqlength, F, Fb, Fm, Fs, Fms,
                     nicks, etaN, mfeStructures, "F", maxILoopSize, mfeEpsilon, FALSE);
       }
    }
    else if( complexity == 5) {
      if( fixedSubOptRange < 0) mfeEpsilon = 0;
      
      bktrF_Fm_FzN5( 0, seqlength - 1, seq, seqlength, F, Fb, Fm, Fp,
                    Fz, Fg, Fgls, Fgrs, Fgl, Fgr, mfeStructures, nicks,
                    etaN, mfeEpsilon, "F");
    }
#ifdef test
    else if( complexity == 4)
       bktrF_Fm_N4( 0, seqlength - 1, seq, seqlength, result, F, Fb, Fm, 
                  , "F");
    else if( complexity == 8) 
       bktrF_Fm_FzN8( 0, seqlength - 1, seq, seqlength, result, F, Fb, Fm, Fp,
                    Fz, Fg, thepairs, "F");
#endif
    }

    if( mfeStructures->nStructs >= 1) {
      //correct energies for symmetries
      DBL_TYPE minimum_energy = NAD_INFINITY;
      for( i = 0; i < mfeStructures->nStructs; i++) {
        mfeStructures->validStructs[i].slength = mfeStructures->seqlength;
        mfeStructures->validStructs[i].correctedEnergy = result + mfeStructures->validStructs[i].error +
         LOG_FUNC( (DBL_TYPE) checkSymmetry( (mfeStructures->validStructs)[i].theStruct, seqlength, nicks, symmetry,
                                            nStrands))*kB*TEMP_K + (BIMOLECULAR + salt_correction) *(nStrands-1);
        if(minimum_energy > mfeStructures->validStructs[i].correctedEnergy) {
          minimum_energy = mfeStructures->validStructs[i].correctedEnergy;
        }
      }

      int offset = 0;
      DBL_TYPE max_energy = minimum_energy + fixedSubOptRange;
      int num_structs = mfeStructures->nStructs;

      if(fixedSubOptRange > 0) {
        for(i = 0 ; i < num_structs; i++) {
          if(mfeStructures->validStructs[i].correctedEnergy <= max_energy) {
            mfeStructures->validStructs[i-offset].theStruct = mfeStructures->validStructs[i].theStruct;
            mfeStructures->validStructs[i-offset].error = mfeStructures->validStructs[i].error;
            mfeStructures->validStructs[i-offset].correctedEnergy = mfeStructures->validStructs[i].correctedEnergy;
            mfeStructures->validStructs[i-offset].slength = mfeStructures->validStructs[i].slength;
          } else {
            free((mfeStructures->validStructs)[i].theStruct);
            (mfeStructures->validStructs)[i].theStruct = NULL;
            offset ++;
            mfeStructures->nStructs --;
          }
        }
      }

      // Commented out: no need to report this.
      //    printf("There are %d structs.\n",mfeStructures->nStructs);

      //sort results by corrected energies
      qsort( mfeStructures->validStructs, mfeStructures->nStructs,
           sizeof( oneDnaStruct), &compareDnaStructs);

      result = mfeStructures->validStructs[0].correctedEnergy; //correct the energy 
 
      // Eliminate nonunique structures (only for MFE calculation)
      if( fixedSubOptRange <= 0) {  
        // Eliminate duplicates, keep the right output permutation
        findUniqueMins(mfeStructures,nicks,symmetry,nStrands,0);

        if (onlyOne) { // mfeStructures has only the first in list
          for( i = 1; i < mfeStructures->nStructs; i++) {
            free( (mfeStructures->validStructs)[i].theStruct);
            (mfeStructures->validStructs)[i].theStruct = NULL;
          }
          mfeStructures->nStructs = 1;
          mfeStructures->nAlloc = 1;
          mfeStructures->minError = 0.0;  
        } else {
          findUniqueMins( mfeStructures, nicks, symmetry, nStrands, 0);
        }
      }

      if(mfe_sort_method == 1) {
        qsort(mfeStructures->validStructs,mfeStructures->nStructs,
          sizeof(oneDnaStruct),&compareDnaStructsOutput);
      }
    } 

  free( seq);
  free( foldparens);
  seq = foldparens = NULL;


  free( F);
  free( Fb);
  free( Fm);

  F = Fb = Fm = NULL;

  if(complexity == 3) {
   free( Fs);
   free( Fms);
   free( Fx);
   free( Fx_1);
   free( Fx_2);
   
   Fs = Fms = Fx = Fx_1 = Fx_2 = NULL;
  }

    if( complexity  >= 5) {
    free( Fp);
    free( Fz);
    free( Fg);

    Fp = Fz = Fg = NULL;

    if( complexity == 5) {
      free( Fgl);
      free( Fgr);
      free( Fgls);
      free( Fgrs);
      free( possiblePairs);
      free( FgIx);
      free( FgIx_1);
      free( FgIx_2);

      Fgl = Fgr = Fgls = Fgrs = FgIx = FgIx_1 = FgIx_2 = NULL;
      possiblePairs = NULL;

      free(sizeTerm);
      sizeTerm = NULL;
    }
  }

  for( i = 0; i <= seqlength-1; i++) {
    for( j = i-1; j <= seqlength-1; j++) {
      pf_ij = pf_index(i,j,seqlength);
      free( etaN[pf_ij]);
    }
  }
  free( etaN);

  free( maxILoopSize); maxILoopSize = NULL;
  free( minILoopEnergyBySize); minILoopEnergyBySize = NULL;


  return result;
}



/* ****************************** */
DBL_TYPE mfe( int seq[], int seqLen, int *thepairs) { 
  //ignores symmetry, single mfe, default parameters for DNA
  return mfeFull( seq, seqLen, thepairs, 3, DNA, 1, 37, 1.0, 0.0, 0);
}


/* ****************************** */
DBL_TYPE mfeFull( int inputSeq[], int seqLen, int *thepairs, int complexity, 
                  int naType, int dangles, DBL_TYPE temperature,
      DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) { 
//ignores symmetry, single mfe
  DBL_TYPE returnValue;
  int i;

  dnaStructures mfeStructures = {NULL, 0, 0, 0, 0}; //all structures withi epsilon of mfe
  returnValue = mfeFullWithSym( inputSeq, seqLen, &mfeStructures, complexity, 
        naType, dangles, temperature, 1, 1, sodiumconc, magnesiumconc,
        uselongsalt);

  for( i = 0; i < mfeStructures.seqlength; i++) {
    thepairs[i] = mfeStructures.validStructs[0].theStruct[i];
  }

  clearDnaStructures( &mfeStructures); 
  return returnValue;

}


/* ****************************** */
DBL_TYPE mfeFullWithSym( int inputSeq[], int seqLen, 
              dnaStructures *mfeStructures, int complexity, int naType, 
        int dangles, DBL_TYPE temperature, int symmetry, int onlyOne,
              DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {

  int fixedSubOptRange = -1; //When < 0, will find mfe structures

  return mfeFullWithSym_SubOpt( inputSeq, seqLen, mfeStructures, complexity, naType, 
        dangles, temperature, symmetry, fixedSubOptRange, onlyOne,
        sodiumconc, magnesiumconc, uselongsalt);
}

/* ****************************** */
//This function converts intpairs to parens notation 
//allocate char *structure to seqlength + 1 before passing in 

void getStructure( int seqlength, const int thepairs[], char *structure) {
  int i;

  for( i = 0; i < seqlength; i++) {
    if( thepairs[i] != -1) {
      if( thepairs[i] > i) {
        structure[i] = '(';
      }
      else {
        structure[i] = ')';
      }
    }
    else {
      structure[i] = '.';
    }
  }
  structure[ seqlength] = '\0';
}


/* ******************** */
void initMfeStructures( dnaStructures *mfeStructures, int seqlength) {

  int i;

  mfeStructures->minError = 0;
  mfeStructures->nStructs = 1;
  mfeStructures->nAlloc = 1;
  mfeStructures->seqlength = seqlength;
  mfeStructures->validStructs = (oneDnaStruct *) malloc( 1*sizeof( oneDnaStruct) ); 
  (mfeStructures->validStructs)[0].error = 0;
  (mfeStructures->validStructs)[0].correctedEnergy = 0;
  (mfeStructures->validStructs)[0].theStruct = (int *) malloc( seqlength*sizeof(int));
  
  for( i = 0; i < seqlength; i++) {
    (mfeStructures->validStructs)[0].theStruct[i] = -1;
  }

}


/* ******************** */
int compareDnaStructs( const void *p1, const void *p2) {
  const oneDnaStruct *s1 = (oneDnaStruct *)p1;
  const oneDnaStruct *s2 = (oneDnaStruct *)p2;

  if( s1->correctedEnergy < s2->correctedEnergy) return -1;
  if( s1->correctedEnergy > s2->correctedEnergy) return 1;
  int st_val = compareDnaStructsOutput(p1,p2);
  return st_val;
}

int compareDnaStructsOutput(const void *p1, const void * p2) {
  const oneDnaStruct *s1 = (oneDnaStruct *)p1;
  const oneDnaStruct *s2 = (oneDnaStruct *)p2;

  int index = 0;
  int length = s1->slength;
  for(index = 0 ; index < length ; index++) {
    if(s1->theStruct[index] < s2->theStruct[index]) {
      return -1;
    }
    if(s1->theStruct[index] > s2->theStruct[index]) {
      return 1;
    }
  }
  return 0;
}
