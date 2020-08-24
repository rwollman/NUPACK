/*
  pf.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007
            Asif Khan 8/2009 Brian Wolfe 10/2009
  This file moves the partition function algorithm to a function, so
  that it can be more readily used as in a library.  In addition,
  scaling is incorporated to allow for longer sequences.
*/

#include "pf.h"

/* ************************************************ */


// Just for a test. Please do not release this.
static DBL_TYPE * pairing_bonuses = NULL;
static int use_bonuses = 0;

/* ******************** */
DBL_TYPE pfuncFullWithBonuses( int inputSeq[], int complexity, int naType, int dangles, 
                    DBL_TYPE temperature, int calcPairs, int perm_symm, DBL_TYPE sodiumconc, 
                    DBL_TYPE magnesiumconc, int uselongsalt, DBL_TYPE * bonuses) {

  pairing_bonuses = bonuses;
  use_bonuses = 1;

  int nStrands;
  int seqlength=getSequenceLengthInt (inputSeq, &nStrands);
  DBL_TYPE res = pfuncFullWithSymHelper(inputSeq, seqlength, nStrands, complexity, naType, 
				dangles, temperature, calcPairs, perm_symm, sodiumconc,
				magnesiumconc, uselongsalt);

  use_bonuses = 0;
  pairing_bonuses = NULL;
  return res;
}


// This is the main function for computing partition functions.
DBL_TYPE pfuncFullWithSymHelper( int inputSeq[], int seqlength, int nStrands,
                                int complexity, int naType, 
				 int dangles, DBL_TYPE temperature, int calcPairs, int permSymmetry,
				 DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {

  //complexity: 3 = N^3, 4 = N^4, 5 = N^5, 8 = N^8
  //naType: DNA = 0, RNA = 1
  //dangles: 0 = none, 1 = normal, 2 = add both

  seqHash=0; // Invalidate ExplDangle cache every time
  int *seq = (int*) calloc( (seqlength+1),sizeof( int) );
  extern int use_cache;
  use_cache=1;

  DBL_TYPE *Q = NULL;
  DBL_TYPE *Qb = NULL;
  DBL_TYPE *Qm = NULL; //O(N^2)
  DBL_TYPE *Qb_bonus = NULL;

  //Multiple strand arrays

  int isPairPrExtern=FALSE;

  //N^3 arrays
  DBL_TYPE *Qx = NULL;
  DBL_TYPE *Qx_1 = NULL;
  DBL_TYPE *Qx_2 = NULL; 
  DBL_TYPE *Qs = NULL;
  DBL_TYPE *Qms = NULL;

  //Pseudoknot arrays
  DBL_TYPE *Qp = NULL;
  DBL_TYPE *Qz = NULL; //O(N^2)
  DBL_TYPE *Qg = NULL; //O(N^4) space

  DBL_TYPE *QgIx = NULL;
  DBL_TYPE *QgIx_1 = NULL;
  DBL_TYPE *QgIx_2 = NULL;
  DBL_TYPE *Qgls = NULL;
  DBL_TYPE *Qgrs = NULL; //O(N^4) 
  DBL_TYPE *Qgl = NULL;
  DBL_TYPE *Qgr = NULL; //O(N^4) space

  //extern DBL_TYPE *sizeTerm;

  //Pair probabilities
  DBL_TYPE *Pb = NULL;
  DBL_TYPE *P = NULL;
  DBL_TYPE *Pm = NULL;
  DBL_TYPE *Pms = NULL;
  DBL_TYPE *Ps = NULL;

  //pseudoknots
  DBL_TYPE *Pz = NULL;
  DBL_TYPE *Pp = NULL;
  DBL_TYPE *Pg = NULL; 
  DBL_TYPE *Pbg = NULL;

  //N^5
  DBL_TYPE *Pgl = NULL;
  DBL_TYPE *Pgr = NULL;
  DBL_TYPE *Pgls = NULL;
  DBL_TYPE *Pgrs = NULL;

  /*  
  The above matrices are dynamically allocated matrices that
  contain partition functions restricted to a subsequence of the
  strand.  Each of the above should be accessed by the call
  to Q[ pf_index(i, j)] indicate the partition function between
  i and j, inclusive. 

  They are described in the paper mentioned above.
  */

  int i, j; // the beginning and end bases for Q;
  int L; //This the length of the current subsequence 
  int pf_ij; //index for O(N^2) matrixes; used to reduce calls to pf_index;
#ifdef NUPACK_SAMPLE
  // variables for sampling
  int sample_count;
  int sample_offset;
  int nNicks;
  int pair_flag;
#endif // NUPACK_SAMPLE
  DBL_TYPE returnValue;


  int iMin;
  int iMax;
  
  //pseudoknots
  extern long int maxGapIndex; 
  //used to minimize memory allocation for fastiloops

  //pseudoknots
  short *possiblePairs; //a speedup for fastiloops (not in paper)

  int nicks[ MAXSTRANDS];  //the entries must be strictly increasing
  for (i=0;i<MAXSTRANDS;i++){
    nicks[i]=-1;
  }
  //nicks[i] = N means a strand ends with base N, and a new one starts at N+1
  // isNicked[n] is 0 if no nick at n, 1 otherwise

  int **etaN;
  int arraySize;

  //assign global variables
  TEMP_K = temperature + ZERO_C_IN_KELVIN;

  DNARNACOUNT = naType;
  DANGLETYPE = dangles;
  SODIUM_CONC = sodiumconc;
  MAGNESIUM_CONC = magnesiumconc;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = uselongsalt;

  for( i = 0; i < MAXSTRANDS; i++) { //initialize nicks array
    nicks[i] = -1;
  }


  processMultiSequence( inputSeq, seqlength, nStrands, seq, nicks);

  if( nStrands >= 2 && complexity >= 5) {
    printf("Warning, pseudoknots not allowed for multi-stranded complexes!");
    printf("  Pseudoknots disabled.\n");
    complexity = 3;
  }
#ifdef NUPACK_SAMPLE
  if(nupack_sample && complexity != 3) {
    printf("Warning, sampling only supported for complexity = 3\n");
    complexity = 3;
  }
#endif

  LoadEnergies();

  if( complexity >= 5) //pseudoknotted
    initPF( seqlength); //precompute values

  // Allocate and Initialize Matrices
  arraySize = seqlength*(seqlength+1)/2+(seqlength+1);
  InitLDoublesMatrix( &Q, arraySize, "Q");
  InitLDoublesMatrix( &Qb, arraySize, "Qb");
  InitLDoublesMatrix( &Qm, arraySize, "Qm");
  InitLDoublesMatrix( &Qb_bonus, arraySize, "Qb_bonus");
  //InitLDoublesMatrix( &Qn, arraySize, "Qn");
  //InitLDoublesMatrix( &Qsn, arraySize, "Qsn");

  etaN = (int**) malloc( arraySize*sizeof( int*));
  InitEtaN( etaN, nicks, seqlength);
  nonZeroInit( Q, seq, seqlength);

  if( complexity == 3) {
    InitLDoublesMatrix( &Qs, arraySize, "Qs");
    InitLDoublesMatrix( &Qms, arraySize, "Qms");
  }

  if( complexity >= 5) {
    InitLDoublesMatrix( &Qp, arraySize, "Qp");
    InitLDoublesMatrix( &Qz, arraySize, "Qz");
    nonZeroInit( Qz, seq, seqlength);

    InitLDoublesMatrix( &Qg, maxGapIndex, "Qg");

    if( complexity == 5) {
      InitLDoublesMatrix( &Qgl, maxGapIndex, "Qgl");
      InitLDoublesMatrix( &Qgr, maxGapIndex, "Qgr");
      InitLDoublesMatrix( &Qgls, maxGapIndex, "Qgls");
      InitLDoublesMatrix( &Qgrs, maxGapIndex, "Qgrs");
      CheckPossiblePairs( &possiblePairs, seqlength, seq);
    }
  }

  if (pairing_bonuses && use_bonuses) {
    for (i = 0; i < seqlength; i++) {
      for (j = i; j < seqlength; j++) {
        Qb_bonus[pf_index(i, j, seqlength)] = pairing_bonuses[i*seqlength + j];
        if (Qb_bonus[pf_index(i, j, seqlength)] > 1.1 || Qb_bonus[pf_index(i, j, seqlength)] < 0.9) {
          // fprintf(stderr, "%i %i: %Lf\n", i, j, - kB * TEMP_K * LOG_FUNC(Qb_bonus[pf_index(i, j, seqlength)]));
        }
      }
    }
  } else {
    for (i = 0; i < seqlength; i++) {
      for (j = i; j < seqlength; j++) {
        Qb_bonus[pf_index(i, j, seqlength)] = 1.0;
      }
    }
  }

  for( L = 1; L <= seqlength; L++) {
    /* Calculate all sub partition functions for
    distance = 0, then 1, then 2.... */

    if( complexity == 3) {
      manageQx( &Qx, &Qx_1, &Qx_2, L-1, seqlength);   
      //allocate/deallocate memory
    }

    if( complexity == 5) {
      manageQgIx( &QgIx, &QgIx_1, &QgIx_2, L-1, seqlength);
      //manageQgIx manages the temporary matrices needed for 
      //calculating Qg_closed in time n^5
    }
    //Default without parallelization
    iMin = 0;
    iMax = seqlength - L; 


    for( i = iMin; i <= iMax; i++) {
      j = i + L - 1;
      pf_ij = pf_index( i, j, seqlength);
      /* Recursions for Qb.  See figure 13 of paper */
      /* bp = base pairs, pk = pseudoknots */
      if( CanPair( seq[ i], seq[ j]) == FALSE) {
        Qb[ pf_ij] = 0.0; //scaling still gives 0
      }
      else {
        Qb[ pf_ij] += Qb_bonus[pf_ij] * ExplHairpin( i, j, seq, seqlength, etaN);

        //no nicked haripins allowed in previous function
        if( complexity == 3) {
          if( etaN[ EtaNIndex(i+0.5, i+0.5, seqlength)][0] == 0 &&
             etaN[ EtaNIndex(j-0.5, j-0.5, seqlength)][0] == 0) {
               //regular multiloop.  No top-level nicks
               
                Qb[ pf_ij] += Qb_bonus[pf_ij] * SumExpMultiloops(i, j, seq, Qms, Qm,
                                              seqlength, etaN);
          }

          if( etaN[ EtaNIndex(i+0.5, j-0.5, seqlength)][0] >= 1) {
            //Exterior loop (created by nick)
            Qb[ pf_ij] += Qb_bonus[pf_ij] * SumExpExteriorLoop( i, j, seq, seqlength, 
                                             Q, nicks, etaN); 
          }
        }

        if( complexity != 3) { //N^4
          // Interior Loop and Multiloop Case
          Qb[ pf_ij] += Qb_bonus[pf_ij] * SumExpInterior_Multi( i, j, seq, seqlength, Qm, Qb);
        }

        if( complexity >= 5) {
          Qb[ pf_ij] += Qb_bonus[pf_ij] * SumExpQb_Pk( i, j, seq, seqlength, Qp, Qm );
        }
      }

      if( complexity == 3) {
        fastILoops( i, j, L, seqlength, seq, etaN, Qb, Qx, Qx_2, Qb_bonus);
      }


#ifdef TEST
      /* Recursions for Qg.  Figures 16, 19 of paper */
      if( complexity == 8) {
        MakeQg_N8( i, j, seq, seqlength, Qg, Qm);
      }
#endif
      
      if( complexity == 5) {
        MakeQg_N5( i, j, seq, seqlength, Qg, Qm, Qgls, Qgrs, QgIx, QgIx_2,
                  possiblePairs);

        //figure 20
        MakeQgls( i, j, seq, seqlength, Qg, Qm, Qgls);
        MakeQgrs( i, j, seq, seqlength, Qg, Qm, Qgrs);

        //figure 18
        MakeQgl(i, j, seq, seqlength, Qg, Qgl, Qz);
        MakeQgr(i, j, seq, seqlength, Qgr, Qgl, Qz);
      }

#ifdef TEST
      /* Recursions for Qp.  figure 15, 17 */
      if( complexity == 8) {
        Qp[ pf_ij] += SumExpQp_N8( i, j, seq, seqlength, Qg, Qz);
      }
#endif

      if( complexity == 5) {
        Qp[ pf_ij] += SumExpQp_N5( i, j, seq, seqlength, Qgl, Qgr,
                                  Qg, Qz);
      }

      if( complexity == 3) {
        /* Recursions for Qms, Qs */
        MakeQs_Qms( i, j, seq, seqlength, Qs, Qms, Qb, nicks, etaN);
        
        /* Recursions for Q, Qm, Qz */
        MakeQ_Qm_N3( i, j, seq, seqlength, Q, Qs, Qms, Qm,
                    nicks,etaN);
      }

      /*if( complexity == 4) {
      MakeQ_Qm_N4( i, j, seq, seqlength, Q, Qm, Qb);
      } */

      if( complexity >= 5) {
        //figures 12, 14
        MakeQ_Qm_Qz(i, j, seq, seqlength, Q, Qm, Qz, Qb, Qp);
      }
    }
  }

  //adjust this for nStrands, symmetry at rank == 0 node
    returnValue = EXP_FUNC( -1*(BIMOLECULAR + SALT_CORRECTION)*(nStrands-1)/(kB*TEMP_K) )*
      Q[ pf_index(0,seqlength-1, seqlength)]/((DBL_TYPE) permSymmetry);




#ifdef NUPACK_SAMPLE
  if(nupack_sample) {
    if(nupack_sample_list == NULL) {
      printf("NULL pointer for structure storage, exiting\n");
      exit(1);
    }
    isPairPrExtern = TRUE;
    if( pairPr == NULL 
        ) {
      pairPr = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), 
                                  sizeof(DBL_TYPE));
      isPairPrExtern = FALSE;
    }


    InitLDoublesMatrix( &P, arraySize, "P");
    InitLDoublesMatrix( &Pb, arraySize, "Pb");
    InitLDoublesMatrix( &Pm, arraySize, "Pm");
    InitLDoublesMatrix( &Pms, arraySize, "Pms");
    InitLDoublesMatrix( &Ps, arraySize, "Ps");
      

    for(sample_count = 0 ; sample_count < nupack_num_samples ; sample_count++) {
      ClearLDoublesMatrix(&P, arraySize, "P");
      ClearLDoublesMatrix(&Pb, arraySize, "Pb");
      ClearLDoublesMatrix(&Pm, arraySize, "Pm");
      ClearLDoublesMatrix(&Pms, arraySize, "Pms");
      ClearLDoublesMatrix(&Ps, arraySize, "Ps");
    
      calculatePairsN3(Q, Qb, Qm, Qms, Qs,
                       &Qx, &Qx_1, &Qx_2, Qb_bonus, P, Pb, Pm, Pms,
                       Ps, seqlength, seq, nicks, etaN);
      // Copy the structure from the matrix to the dot-paren format
      sample_offset = 0;
      nNicks = etaN[EtaNIndex(0-0.5,seqlength-0.5,seqlength)][0];
      // First, fill in with unpaired bases
      for(i = 0 ; i < seqlength + nNicks ; i++) {
        nupack_sample_list[sample_count][i] = '.';
      }
      for(i = 0 ; i < seqlength ; i++) {
        pair_flag = FALSE;
        if(nNicks > sample_offset && nicks[sample_offset] < i) {
          nupack_sample_list[sample_count][i + etaN[EtaNIndex(-0.5,i-1.5,seqlength)][0]] 
          = '+';
          sample_offset ++;
        }
        for(j = i+1; j < seqlength ; j++) {
          if(Pb[pf_index(i,j,seqlength)] > 0.5) {
          // Really going to be either 0 or 1, but whatever
            pair_flag = TRUE;
            // The +etaN is just to make sure that I account for + in the sequence
            nupack_sample_list[sample_count][i+etaN[EtaNIndex(-0.5,i-0.5,seqlength)][0]] 
            = '(';
            nupack_sample_list[sample_count][j+etaN[EtaNIndex(-0.5,j-0.5,seqlength)][0]]  
            = ')';
          } 
        }
      }
      nupack_sample_list[sample_count][seqlength + nNicks] = '\0';
    }
    free(P);
    free(Pb);
    free(Pm);
    free(Pms);
    free(Ps);
  } 
#endif //NUPACK_SAMPLE
  //Calculate Pair Probabilities as needed 
  if(calcPairs) {
    /*
    if( complexity == 4) {
    InitLDoublesMatrix( &Pb, arraySize, "Pb");
    calculatePairsN4( Q, Qb, Qm, Pb, seqlength, seq);
    }
    */
    isPairPrExtern = TRUE;
    if( pairPr == NULL 
        ) {
      pairPr = (DBL_TYPE*) calloc( (seqlength+1)*(seqlength+1), 
                                  sizeof(DBL_TYPE));
      isPairPrExtern = FALSE;
    #ifdef NUPACK_SAMPLE
    } else if(nupack_sample) {
      isPairPrExtern = FALSE;
    #endif //NUPACK_SAMPLE
    }
    if( complexity == 3) {
      InitLDoublesMatrix(  &P, arraySize, "P");
      InitLDoublesMatrix(  &Pb, arraySize, "Pb");
      InitLDoublesMatrix(  &Pm, arraySize, "Pm");
      InitLDoublesMatrix(  &Pms, arraySize, "Pms");
      InitLDoublesMatrix(  &Ps, arraySize, "Ps");
      
      calculatePairsN3( Q, Qb, Qm, Qms, Qs, 
                       &Qx, &Qx_1, &Qx_2, Qb_bonus, P, Pb, Pm, Pms,
                       Ps, seqlength, seq, nicks, etaN);
    }

    /*
    if( complexity == 8) {
    InitLDoublesMatrix( &P, arraySize, "P");
    InitLDoublesMatrix( &Pb, arraySize, "Pb");
    InitLDoublesMatrix( &Pm, arraySize, "Pm");
    InitLDoublesMatrix( &Pp, arraySize, "Pp");
    InitLDoublesMatrix( &Pz, arraySize, "Pz");
    InitLDoublesMatrix( &Pg, 
    maxGapIndex,
    "Pg");
    InitLDoublesMatrix( &Pbg, arraySize, "Pbg");
    
    P[ pf_index( 0, seqlength-1, seqlength)] = 1.0;
    //calculatePairsN8( Q, Qb, Qm, Qp, Qz, Qg, P, Pb, Pp, Pz, Pg, Pbg, Pm,seqlength, seq);
    }
    */

    if( complexity == 5) {
      InitLDoublesMatrix( &P, arraySize, "P");
      InitLDoublesMatrix( &Pb, arraySize, "Pb");
      InitLDoublesMatrix( &Pm, arraySize, "Pm");
      InitLDoublesMatrix( &Pp, arraySize, "Pp");
      InitLDoublesMatrix( &Pz, arraySize, "Pz");
      InitLDoublesMatrix( &Pg, maxGapIndex, "Pg");
      InitLDoublesMatrix( &Pgl, maxGapIndex, "Pgl");
      InitLDoublesMatrix( &Pgr, maxGapIndex, "Pgr");
      InitLDoublesMatrix( &Pgls, maxGapIndex, "Pgls");
      InitLDoublesMatrix( &Pgrs, maxGapIndex,"Pgrs");
      InitLDoublesMatrix( &Pbg, seqlength*(seqlength+1)/2, "Pbg");

      P[ pf_index( 0, seqlength-1, seqlength)] = 1.0;
      
      calculatePairsN5( Q, Qb, Qm, Qp, Qz, Qg, Qgl, Qgr, Qgls, Qgrs, &QgIx,
                       &QgIx_1, &QgIx_2, P, Pb, Pp, Pz, Pg, Pbg, Pm, Pgl, Pgr, 
                       Pgls, Pgrs, seqlength, seq);
    }
  }

  if(EXTERN_Q != NULL) {
    // Fill out the external Q and Qb matrices.
    // This is used by the design code to do approximations to the energy function
    for(i = 0 ; i < seqlength ; i++) {
      for(j = i ; j < seqlength ; j++) {
        EXTERN_Q[pf_index(i,j,seqlength)] = Q[pf_index(i,j,seqlength)];
      }
    }
  } else {
    EXTERN_Q = NULL;
  }

  if(EXTERN_QB != NULL) {
    // Fill out the external Q and Qb matrices.
    // This is used by the design code to do approximations to the energy function
    for(i = 0 ; i < seqlength ; i++) {
      for(j = 0 ; j < seqlength + 1; j++) {
        EXTERN_QB[pf_index(i,j,seqlength)] = Qb[pf_index(i,j,seqlength)];
      }
    }
  } else {
    EXTERN_QB = NULL;
  }


  free( Q);
  free( Qb);
  free( Qm);
  free( Qb_bonus);

  Q = Qb = Qm = NULL;

  if( complexity == 3) {
    free( Qs);
    free( Qms);
    
    free( Qx);
    free( Qx_1);
    free( Qx_2);
    
    Qs = Qms = Qx = Qx_1 = Qx_2 = NULL;
  }

  if( complexity >= 5) {
    free( Qp);
    free( Qz);
    free( Qg);

    Qp = Qz = Qg = NULL;

    if( complexity == 5) {
      free(Qgl);
      free(Qgr);
      free(Qgls);
      free(Qgrs);
      free(QgIx);
      free(QgIx_1);
      free(QgIx_2);
      free(possiblePairs);
      
      Qgl = Qgr = Qgls = Qgrs = QgIx = QgIx_1 = QgIx_2 = NULL;
      possiblePairs = NULL;
      
      free(sizeTerm);
      sizeTerm = NULL;
    }
  }

  if( calcPairs) {
    if( isPairPrExtern == FALSE) {
      free( pairPr);
      pairPr = NULL;
    }

    /*
    if( complexity == 4) {
    free( Pb);
    Pb = NULL;
    }
    */

    if( complexity == 3) {
      free(   P);
      free(  Pb);
      free(  Pm);
      free( Pms);
      free(  Ps);
      
      P = Pb = Pm = Pms = Ps = NULL;
    }

    if( complexity == 8) {
      free(P);
      free(Pb);
      free(Pz);
      free(Pp);
      free(Pg);
      free(Pbg);
      free(Pm);
      P = Pb = Pz = Pp = Pg = Pbg = Pm = NULL;
    }

    if( complexity == 5) {
      free(P);
      free(Pb);
      free(Pz);
      free(Pp);
      free(Pg);
      free(Pbg);
      free(Pm);
      free(Pgl);
      free(Pgr);
      free(Pgls);
      free(Pgrs);
      P = Pb = Pz = Pp = Pg = Pbg = Pm = Pgl = Pgr = Pgls = Pgrs = NULL;
    }
  }

  free( seq);

  for( i = 0; i <= seqlength-1; i++) {
    for( j = i-1; j <= seqlength-1; j++) {
      pf_ij = pf_index(i,j,seqlength);
      free( etaN[pf_ij]);
    }
  }

  free( etaN);

  return returnValue;
}
/* ****** */


// This is the main function for computing partition functions.
DBL_TYPE pfuncFullWithSym( int inputSeq[], int complexity, int naType, 
			   int dangles, DBL_TYPE temperature, int calcPairs, int permSymmetry,
			   DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {
  int nStrands;
  int seqlength=getSequenceLengthInt (inputSeq, &nStrands);
  return pfuncFullWithSymHelper(inputSeq, seqlength, nStrands, complexity, naType, 
				dangles, temperature, calcPairs, permSymmetry, sodiumconc,
				magnesiumconc, uselongsalt);
  }

/* ******************** */
DBL_TYPE pfunc( int seq[]) {
  // Does a DNA pfunc calculation with everything else set to defaults
  return pfuncFull( seq, 3, DNA, 1, 37, 1, 1.0, 0.0, 0);
}

/* ******************** */
DBL_TYPE pfuncFull( int inputSeq[], int complexity, int naType, int dangles, 
                    DBL_TYPE temperature, int calcPairs,
		    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt) {
  return pfuncFullWithSym( inputSeq, complexity, naType,
                          dangles, temperature,
			   calcPairs, 1, sodiumconc, magnesiumconc, uselongsalt);
}


