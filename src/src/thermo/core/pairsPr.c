/*
  pairsPr.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 8/2002, Justin Bois 1/2007
    Asif Khan 8/2009 Brian Wolfe 10/2009

  Implementation of the base pairing probability matrix (McCaskill,
  1990)

  Originally, just for non-pseudoknotted case, but I have added the
  pseudoknot version as well.
*/

#include "pairsPr.h"

struct {
  // Current sum
  DBL_TYPE pq;
  DBL_TYPE pqb;
  DBL_TYPE pqm;
  DBL_TYPE pqs;
  DBL_TYPE pqms;
  DBL_TYPE pqx;

  // Threshold
  DBL_TYPE Zq;
  DBL_TYPE Zqb;
  DBL_TYPE Zqm;
  DBL_TYPE Zqs;
  DBL_TYPE Zqms;
  DBL_TYPE Zqx;

  int pqset;
  int pqbset;
  int pqmset;
  int pqsset;
  int pqmsset;
  int pqxset;
} samplingGlobals ;


/* **************** */
long double log2l(long double);
float subtractLongDouble( DBL_TYPE *ap, DBL_TYPE b) {
  DBL_TYPE a = *ap;
  float precisionLost;
  DBL_TYPE value;

  if( b >= a) {
    *ap = 0;
    return 100;
  }

  value = -b/a + 1.0;

  precisionLost = (float) -1*log2l(value);

#ifdef DEBUG
  if( a < 0) printf("%Le %Le %f\n", (long double) a, (long double) b, precisionLost);
#endif

  *ap = a-b;
#ifdef DEBUG
  if( *ap < 0) printf("! %Le %Le %Le %f\n", (long double) a, (long double) b,
                      (long double) a-b, precisionLost);
#endif
  return precisionLost;
}
/* ******************* */


void calculatePairsN3( DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE *Qm, DBL_TYPE *Qms,
                       DBL_TYPE *Qs, DBL_TYPE **Qx,
                       DBL_TYPE **Qx_1, DBL_TYPE **Qx_2,
                       DBL_TYPE *Qb_bonus,
                       DBL_TYPE *P, DBL_TYPE *Pb, DBL_TYPE *Pm, DBL_TYPE *Pms,
                       DBL_TYPE *Ps, int seqlength,
                       int seq[], int *nicks, int** etaN) {

  int L, i, j, indI;
  DBL_TYPE rowsum;
  int pf_ij;
  DBL_TYPE *Px, *Px_1, *Px_2;
  float *preX, *preX_1, *preX_2;
  int iMin, iMax;

#ifdef NUPACK_SAMPLE
  samplingGlobals.pq = samplingGlobals.pqb = samplingGlobals.pqm = 
    samplingGlobals.pqs = samplingGlobals.pqms = samplingGlobals.pqx = 0.0;
  samplingGlobals.Zq = samplingGlobals.Zqb = samplingGlobals.Zqm = 
    samplingGlobals.pqs = samplingGlobals.Zqms = samplingGlobals.Zqx = 0.0;

  samplingGlobals.pqset = samplingGlobals.pqbset = samplingGlobals.pqmsset = 0;
  samplingGlobals.pqmset = samplingGlobals.pqxset = samplingGlobals.pqsset = 0;
#endif // NUPACK_SAMPLE

  Px = Px_1 = Px_2 = NULL;
  preX = preX_1 = preX_2 = NULL;

  P[ pf_index(0, seqlength-1, seqlength)] = 1;
  //Pn[ pf_index(0, seqlength-1, seqlength)] = 1;

  for( L = seqlength; L >= 1; L--) {

    prManageQx( Qx, Qx_1, Qx_2, &Px, &Px_1, &Px_2,
                &preX, &preX_1, &preX_2, L-1, seqlength);
    // Break up the sequence for parallelization
    iMin = 0; iMax = seqlength - L;

    for( i = iMin; i <= iMax ; i++) {
#ifdef NUPACK_SAMPLE
      if(nupack_sample) {
        samplingGlobals.Zq =(DBL_TYPE) genrand_res53();
        samplingGlobals.Zqb =(DBL_TYPE) genrand_res53() ;
        samplingGlobals.Zqm = (DBL_TYPE) genrand_res53();
        samplingGlobals.Zqms = (DBL_TYPE) genrand_res53();
        samplingGlobals.Zqx = (DBL_TYPE) genrand_res53();
        samplingGlobals.Zqs = (DBL_TYPE) genrand_res53();


        samplingGlobals.pq = 0;
        samplingGlobals.pqb = 0;
        samplingGlobals.pqm = 0;
        samplingGlobals.pqms = 0;
        samplingGlobals.pqx = 0;
        samplingGlobals.pqs = 0;

        samplingGlobals.pqset = 0;
        samplingGlobals.pqbset = 0;
        samplingGlobals.pqmset = 0;
        samplingGlobals.pqmsset = 0;
        samplingGlobals.pqxset = 0;
        samplingGlobals.pqsset = 0;
      }
#endif // NUPACK_SAMPLE

      j = i + L - 1;
      pf_ij = pf_index(i, j, seqlength);

      MakeP_Pm_N3(i,j,seq,seqlength, Q, Qs, Qms, Qm,
                  P, Ps, Pms, Pm, etaN);
      MakePs_Pms(i,j,seq,seqlength, Qs, Qms, Qb,
                 Ps, Pms, Pb, nicks, etaN);
      pf_ij = pf_index(i,j,seqlength);
#ifdef USE_N4_INTLOOPS
      prInteriorLoopsN4MS( i, j, seq, seqlength, Qb, Pb, nicks);

#else

      prFastILoops( i, j, L, seqlength, seq, Qb, *Qx, *Qx_2, Qb_bonus, Pb, Px,
                    Px_2, nicks, etaN, preX, preX_2);
#endif


      if( etaN[ EtaNIndex(i+0.5, j-0.5, seqlength)][0] >= 1) {
        prExterior_N3(i, j, seq, seqlength,
                      Q, Qb, Qb_bonus, P, Pb,
                      nicks, etaN);
      }

      if( etaN[ EtaNIndex_same(i+0.5, seqlength)][0] == 0 &&
          etaN[ EtaNIndex_same(j-0.5, seqlength)][0] == 0) {
        prMultiBp_N3(i, j, seq, seqlength, Qb, Qms, Qm, Qb_bonus,
                     Pb, Pms, Pm, etaN);
      }
    }
  }

  //check rowsums:
  // only for rank 0
    for( i = 0; i <= seqlength - 1; i++) {
      rowsum = 0;
      for ( j = 0; j <= seqlength - 1; j++) {
        if( i <= j) {
          rowsum += Pb[  pf_index(i, j, seqlength)];
          if( Pb[ pf_index(i,j,seqlength)] > 1.0 + NUM_PRECISION) {
            printf("Numerical precision loss in pairsPr.c\n");
            printf("P(%d,%d) = %Le !!\n", i, j, (long double) Pb[ pf_index(i,j,seqlength)]);
          }
        }
        else {
          rowsum += Pb[  pf_index(j, i, seqlength)];
        }
      }
  
      if( rowsum > 1.0 + NUM_PRECISION) {
        printf("Numerical precision loss in pairsPr.c\n");
        printf( "rs: %d %.16Le\n", i, (long double) rowsum);
        printf("Error!!!! %d\n", seq[0]);
        exit(1);
      }
    }

    //store values in pairPr
    for( i = 0; i <= seqlength-1; i++) {
      rowsum = 0;
      indI = i*(seqlength+1);
      for( j = 0; j <= seqlength; j++) {
        if( j == seqlength) {
          pairPr[ indI + j] = 1 - rowsum;
        }
        else if ( i <= j) {
          pf_ij = pf_index(i,j, seqlength);
          rowsum += Pb[ pf_ij];
          pairPr[ indI + j] = Pb[ pf_ij];
        } else {
          pf_ij = pf_index(j,i, seqlength);
          rowsum += Pb[ pf_ij];
          pairPr[ indI + j] = Pb[ pf_ij];
        }
      }
    }

  free( Px); Px = NULL;
  free( Px_1); Px_1 = NULL;
  free( Px_2); Px_2 = NULL;

  free( preX);
  free( preX_1);
  free( preX_2);
  preX = preX_1 = preX_2 = NULL;

}

/* ************************ */

void prFastILoops( int i, int j, int L, int seqlength, int seq[],
                   DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2,
                   DBL_TYPE *Qb_bonus,
                   DBL_TYPE *Pb, DBL_TYPE *Px, DBL_TYPE *Px_2,
                   int *nicks, int **etaN, float *preX, float *preX_2) {

  int d, e, L1, L2, pf_de;
  int pf_ij = pf_index(i,j,seqlength);
  DBL_TYPE pr = 0;
  int size;
  //extern DBL_TYPE loop37[];
  DBL_TYPE oldSizeEnergy;
  DBL_TYPE newSizeEnergy;
  DBL_TYPE energy;
  int fbix, fbix2;
  float precisionLost;

  int index_ij = EtaNIndex( i+0.5, j-0.5, seqlength);
  int nNicks = etaN[ index_ij][0];
  int leftNick, rightNick;

  int isEndNicked = FALSE;
  DBL_TYPE oldValue;

  leftNick = rightNick = -1;

  if( nNicks >= 1) {
    leftNick = nicks[ etaN[ index_ij][1]];
    rightNick = nicks[ etaN[ index_ij][1]+nNicks-1];
  }

  //End nicks are important for extensibility
  if( etaN[ EtaNIndex_same(i-0.5,seqlength)][0] == 1 ||
      etaN[ EtaNIndex_same(j+0.5,seqlength)][0] == 1) {
    isEndNicked = TRUE;
  }

  //L1, L2 <= 3 (directly add to Pb)
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    for( L1 = 0; L1 <= 3; L1++) {
      d = i + L1 + 1;
      for( L2 = 0; L2 <= MIN( 3, j-d-2); L2++) {
        e = j - L2 - 1;

        smallInteriorLoop( pf_ij, seq, seqlength, i, j, d, e, leftNick, rightNick, Qb, 
            Qb_bonus, Pb, 3);
      }
    }

    /* Next consider large bulges or large asymmetric loops */
    // Case 2a  L1 = 0,1,2,3, L2 >= 4;
    for( L1 = 0; L1 <= 3; L1++) {
      d = i + L1 + 1;
      for( L2 = 4; L2 <= j - d - 2; L2++) {
        e = j - L2 - 1;

        smallInteriorLoop( pf_ij, seq, seqlength, i, j, d, e, leftNick, rightNick, Qb, 
            Qb_bonus, Pb, 4);
      }
    }

    // Case 2b L1 >= 4, L2 = 0,1,2,3;
    for( L2 = 0; L2 <= 3; L2++) {
      e = j - L2 - 1;
      for( L1 = 4;  L1 <= e - i - 2; L1++) {
        d = i + L1 + 1;

        smallInteriorLoop( pf_ij, seq, seqlength, i, j, d, e, leftNick, rightNick, Qb, 
            Qb_bonus, Pb, 5);

      }
    }
  } //end canpair(i,j)


  //Add in special "root" cases ( i == 0 || j == seqlength - 1) to Qx
  // or i-1 nicked or j nicked must be added
  if( (i == 0 || j == seqlength - 1 || isEndNicked == TRUE) &&
      j >= i+11 /* && L <= seqlength - 2*/ ) {
    for( d = i + 5; d <= j - 6; d++) {
      if( leftNick != -1 && leftNick <= d-1) break;

      for( e = d + 1; e <= j - 5; e++) {
        if( rightNick != -1 && rightNick >= e) continue;

        if( CanPair( seq[d], seq[e]) == TRUE) {

          L1 = d - i - 1;
          L2 = j - e - 1;
          size = L1 + L2;

          pf_de = pf_index( d,e,seqlength);

          energy = asymmetryEfn( L1, L2, size);
          energy += InteriorMM( seq[e], seq[d], seq[e+1], seq[d-1]);
          /*Exclude the i-j stacking energy here, just in case i-j
            don't pair */


          fbix =  fbixIndex( j-i, i, size, seqlength);

          Qx[ fbix ] +=
            EXP_FUNC(-energy/(kB*TEMP_K))*Qb[ pf_de];


        }
      }
    }
  }//end special cases


  //use Qb to calculate Px
  if( CanPair( seq[ i], seq[j]) == TRUE && Qb[ pf_ij] > 0) {
    for( size = 8; size <= L - 4; size++) {

      fbix =  fbixIndex( j-i, i, size, seqlength);
      pr = Pb[ pf_ij] * Qx[ fbix ] * Qb_bonus[pf_ij] * 
        EXP_FUNC( -InteriorMM( seq[i], seq[j], seq[i+1],
                           seq[j-1])/(kB*TEMP_K))/ Qb[ pf_ij];

      oldValue = Px[ fbix];
#ifdef NUPACK_SAMPLE
      if(!nupack_sample) {
#endif // NUPACK_SAMPLE
      Px[ fbix] += pr;
#ifdef NUPACK_SAMPLE
      } else {
        samplingGlobals.pqb += pr;
        if(!samplingGlobals.pqbset && samplingGlobals.pqb >= samplingGlobals.Zqb) {
          samplingGlobals.pqbset = 1;
          Px[fbix] = 1;
        }
      }
#endif


      if( (pr > 1.0 + NUM_PRECISION ) ) {
        printf("Numerical precision loss in pairsPr.c\n");
        printf("Error in precision due to subtractions!: ");
        printf("6 %Le %d %d %d\n", (long double) pr, i, j, seq[0]);
        exit(1);
      }


    }
  }

  //Next calculate Pb using Px and Qx for L1 == 4 || L2 == 4
  //Case 1:  L1 = 4, L2 >= 4;
  if( L >= 12 && (leftNick == -1 || leftNick >= i+5) ) {
    L1 = 4;
    d = i + L1 + 1;
    for( L2 = 4; L2 <= j - d - 2; L2++) {
      size = L1 + L2;
      e = j - L2 - 1;
      fbix =  fbixIndex( j-i, i, size, seqlength);
      if( rightNick != -1 && rightNick >= e) continue;

      if( CanPair( seq[d], seq[e]) == TRUE) {

        energy = asymmetryEfn( L1, L2, size);
        energy += InteriorMM( seq[e], seq[d], seq[e+1], seq[d-1]);
        /*Exclude the i-j stacking energy here, just in case i-j
          don't pair */

        pf_de = pf_index(d,e,seqlength);
        if( Qx[ fbix] > 0) {

          if( (pr > 1.0 + NUM_PRECISION ) ) {
            printf("Numerical precision loss in pairsPr.c\n");
            printf("6a! %Le\n", (long double) Px[ fbix]);
          }

          pr = Px[ fbix ] *
            EXP_FUNC(-energy/(kB*TEMP_K))*Qb[ pf_de] / Qx[ fbix];
#ifdef NUPACK_SAMPLE
          if(!nupack_sample) {
#endif
            Pb[ pf_de] += pr;
            Px[ fbix] -= pr;
#ifdef NUPACK_SAMPLE
          } else {
            samplingGlobals.pqx += pr;
            if( (pr > 1.0 + NUM_PRECISION ) ) {
              printf("Numerical precision loss in pairsPr.c\n");
              printf("8! Px[fbix] = %Le\n", (long double) Px[ fbix]);
              printf("pr - 1 = %Le  energy = %Le  Qb[pf_de] = %Le  Qx[fbix] = %Le \n",pr - 1.0, energy,Qb[pf_de] ,Qx[fbix]);
              exit(1);
            }
            pr = 0;

            if(!samplingGlobals.pqxset && samplingGlobals.pqx >= samplingGlobals.Zqx) {
              samplingGlobals.pqxset = 1;
              Pb[pf_de] = 1;
              Px[fbix] = 0;
            }
          }
#endif


          if( (pr > 1.0 + NUM_PRECISION ) ) {
            printf("Numerical precision loss in pairsPr.c\n");
            printf("7\n");
            exit(1);
          }
        }

        precisionLost = subtractLongDouble( &(Qx[ fbix]),
                                            //scale( j-i+d-e) *
                                            EXP_FUNC(-energy/(kB*TEMP_K))*
                                            Qb[ pf_de]);


        preX[ fbix] += precisionLost;
        if( preX[ fbix] >= MAXPRECERR) {
          recalculateQx( i, j, size, fbix, seq, seqlength, Qx, Qb, nicks, etaN,  1);
          preX[ fbix] = 0.0;
        }

      }
    }
  }
  if( L>=12 && (rightNick == -1 || rightNick <= j - 6) ) {
    //Case 2  L1 > 4, L2 = 4
    L2 = 4;
    e = j - L2 -1;
    for( L1 = 5; L1 <= e-i-2; L1++) {
      size = L1 + L2;
      d = i + L1 + 1;
      if( leftNick != -1 && leftNick <= d-1) break;

      fbix =  fbixIndex( j-i, i, size, seqlength);

      if( CanPair( seq[d], seq[e]) == TRUE) {

        energy = asymmetryEfn( L1, L2, size);
        energy += InteriorMM( seq[e], seq[d], seq[e+1], seq[d-1]);

        pf_de = pf_index(d,e,seqlength);
        if( Qx[ fbix] > 0) {


          if( (pr > 1.0 + NUM_PRECISION ) ) {
            printf("Numerical precision loss in pairsPr.c\n");
            printf("7a! %Le\n", (long double) Px[ fbix]);
          }


          pr = ((Px[ fbix ] *
            (EXP_FUNC(-energy/(kB*TEMP_K))))*Qb[ pf_de]) / Qx[ fbix];
          #ifdef NUPACK_SAMPLE
          if(!nupack_sample) {
          #endif // NUPACK_SAMPLE
            Pb[ pf_de] += pr;
            Px[ fbix] -= pr;
          #ifdef NUPACK_SAMPLE
          } else {
            samplingGlobals.pqx += pr;
            if( (pr > 1.0 + NUM_PRECISION ) ) {
              printf("Numerical precision loss in pairsPr.c\n");
              printf("8! Px[fbix] = %Le\n", (long double) Px[ fbix]);
              printf("pr - 1 = %Le  energy = %Le  Qb[pf_de] = %Le  Qx[fbix] = %Le \n",pr - 1.0, energy,Qb[pf_de] ,Qx[fbix]);
              exit(1);
            }
            pr = 0;

            if( (pr > 1.0 + NUM_PRECISION ) ) {
              printf("Numerical precision loss in pairsPr.c\n");
              printf("8! Px[fbix] = %Le\n", (long double) Px[ fbix]);
              printf("pr - 1 = %Le  energy = %Le  Qb[pf_de] = %Le  Qx[fbix] = %Le \n",pr - 1.0, energy,Qb[pf_de] ,Qx[fbix]);
              exit(1);
            }
            pr = 0;
            if(!samplingGlobals.pqxset && samplingGlobals.pqx >= samplingGlobals.Zqx) {
              samplingGlobals.pqxset = 1;
              Pb[pf_de] = 1;
              Px[fbix] = 0;
            }
          }
          #endif //NUPACK_SAMPLE



          if( (pr > 1.0 + NUM_PRECISION ) ) {
            printf("Numerical precision loss in pairsPr.c\n");
            printf("8! Px[fbix] = %Le\n", (long double) Px[ fbix]);
            printf("pr - 1 = %Le  energy = %Le  Qb[pf_de] = %Le  Qx[fbix] = %Le \n",pr - 1.0, energy,Qb[pf_de] ,Qx[fbix]);
            exit(1);
          }
        }

        precisionLost = subtractLongDouble( &(Qx[ fbix]),
                                            //scale( j-i+d-e) *
                                            EXP_FUNC(-energy/(kB*TEMP_K))*
                                            Qb[ pf_de]);

        preX[ fbix] += precisionLost;
        if( preX[ fbix] >= MAXPRECERR) {
          recalculateQx( i, j, size, fbix, seq, seqlength, Qx, Qb, nicks, etaN, 2);
          preX[ fbix] = 0.0;
        }

      }
    }
  }

  //contraction
  for( size = 10; size <= L - 4; size++) {
    if( size <= 30) {
      oldSizeEnergy = loop37[ size - 1];
    }
    else {
      oldSizeEnergy = loop37[ 30 - 1];
      oldSizeEnergy += sizeLog (size); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);
    }

    if( size - 2 <= 30) {
      newSizeEnergy = loop37[ size-2 - 1];
    }
    else {
      newSizeEnergy = loop37[ 30 - 1];
      newSizeEnergy += sizeLog (size-2); //1.75*kB*TEMP_K*LOG_FUNC( (size-2)/30.0);
    }

    fbix  = fbixIndex( j-i, i, size, seqlength);
    fbix2 = fbixIndex( j-i-2, i+1, size-2, seqlength);

    Qx_2[ fbix2] =
      Qx[ fbix] * EXP_FUNC( -(newSizeEnergy - oldSizeEnergy)/(kB*TEMP_K));

    Px_2[ fbix2] = Px[ fbix];
    preX_2[ fbix2] = preX[ fbix];

  }
}

/* ********** */
void smallInteriorLoop( int pf_ij, int seq[], int seqlength, int i, int j,
                        int d, int e, int leftNick, int rightNick, DBL_TYPE *Qb,
                        DBL_TYPE * Qb_bonus, 
                        DBL_TYPE *Pb, int error) {

  DBL_TYPE energy;
  int pf_de;
  DBL_TYPE pr;

  if( CanPair( seq[d], seq[e]) == TRUE &&
      (leftNick == -1 || leftNick >= d) &&
      (rightNick == -1 || rightNick <= e-1) ) {

    pf_de = pf_index( d, e, seqlength);

    energy = InteriorEnergy( i, j, d, e, seq);

    if( Qb[ pf_ij] > 0) {
      pr = Pb[ pf_ij] * EXP_FUNC( -energy/(kB*TEMP_K)) *
        Qb[ pf_de] * Qb_bonus[pf_ij] / Qb[ pf_ij];

      //printf("%i %i %i %i %i %Le %Le %Le\n",i,j,d,e,samplingGlobals.pqbset,samplingGlobals.pqb,samplingGlobals.Zqb,(long double)pr);
      #ifdef NUPACK_SAMPLE
      if(!nupack_sample) {
      #endif // NUPACK_SAMPLE
        Pb[ pf_de] += pr;
      #ifdef NUPACK_SAMPLE
      } else {
        samplingGlobals.pqb += pr;
        if(!samplingGlobals.pqbset && samplingGlobals.pqb >= samplingGlobals.Zqb) {
          samplingGlobals.pqbset = 1;
          Pb[pf_de] = 1;
        }
      }
      #endif // NUPACK_SAMPLE

      if( (pr > 1.0 + NUM_PRECISION ) ) {
        printf("Numerical precision loss in pairsPr.c\n");
        printf("%d %d %d %d %d %Le %Le %Le %Le\n",
               error, i,d,e,j,
               (long double) Pb[ pf_ij],
               (long double) EXP_FUNC( -energy/(kB*TEMP_K)),
               (long double) Qb[ pf_de],
               (long double) Qb[ pf_ij]);
        exit(1);
      }

    }
  }
}

/* *********** */


void recalculateQx( int i, int j, int size, int fbix, int seq[],
                    int seqlength, DBL_TYPE *Qx, DBL_TYPE *Qb,
                    int *nicks, int **etaN,
                    int side) {

  int d,e;
  DBL_TYPE tmp = 0;

  int L1min = 5;
  int L2min = 4;

  int L1max = j-i-10;
  int leftN = -1, rightN = -1;

  int etaIndex = EtaNIndex( i+0.5, j-0.5, seqlength);
  int nNicks;

  if( side == 2) L2min = 5;

  nNicks =  etaN[ etaIndex][0];
  if( nNicks >= 1) {
    leftN = nicks[ etaN[ etaIndex][1] ];
    L1max = leftN-i-1;
    rightN = nicks[ etaN[ etaIndex][1] + nNicks - 1];
  }

  for( d = i+L1min+1; d <= i+L1max+1; d++) {
    e = d-i+j-2-size;
    if( j-e-1>=L2min && e >= rightN+1 &&
        CanPair( seq[d], seq[e]) == TRUE) {

      if( leftN == -1 && e-d<= 3) continue;

      tmp += Qb[ pf_index(d,e,seqlength)]*
        EXP_FUNC( -1*(InteriorEnergyFull( i,j,d,e,seq, FALSE))/
              (kB*TEMP_K));
    }
  }

  Qx[ fbix] = tmp;
}


/* ****************************************************** */
void prManageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1,
                 DBL_TYPE **Qx_2, DBL_TYPE **Px, DBL_TYPE **Px_1,
                 DBL_TYPE **Px_2,
                 float **preX, float **preX_1, float **preX_2,
                 int len, int seqlength) {
  // Allocate and deallocate QbIx matrices
  // since flow is reversed in Pr calulations, Qx_1 and Qx_2 now
  // now represent values for L-1 and L-2 (rather than +1 and +2)

  int i;
  int maxStorage;
  DBL_TYPE *temp;
  float *tempf;
  //  int c;

  maxStorage = (seqlength - len)*(len - 1);
  for( i = len-1; i >= len - 2; i--) {
    maxStorage = MAX( maxStorage, ( seqlength - i)*(i - 1) );
  }


  if( len == seqlength-1 && len >= 11) { //first use of these matrices

    free( *Qx);
    free( *Qx_1);
    free( *Qx_2);

    *Qx = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Qx_1 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Qx_2 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    /* old
       free( *Qx_2);
       *Qx_2 =
       (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    */

    *Px = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Px_1 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Px_2 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    *preX = (float *) calloc( maxStorage, sizeof( float));
    *preX_1 = (float *) calloc( maxStorage, sizeof( float));
    *preX_2 = (float *) calloc( maxStorage, sizeof( float));

    if( *Qx_2 == NULL || *Px == NULL || *Px_1 == NULL || *Px_2 == NULL
        || *preX_1 == NULL || *preX_1 == NULL || *preX_2 == NULL) {
      printf("Error in Qx, Px, preX allocation\n");
    }

  }
  else if( len >= 11) {
    temp = *Px;
    *Px = *Px_1;
    *Px_1 = *Px_2;
    free( temp);
    *Px_2 =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    temp = *Qx;
    *Qx = *Qx_1;
    *Qx_1 = *Qx_2;
    free( temp);
    *Qx_2 =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    tempf = *preX;
    *preX = *preX_1;
    *preX_1 = *preX_2;

    free( tempf);
    tempf = NULL;

    *preX_2 =
      (float *) calloc( maxStorage, sizeof( float));

    if( *Qx_2 == NULL || *Px_2 == NULL || *preX_2 == NULL) {
      printf("Error in Qx_2, Px_2, preX_2 allocation\n");
    }

  }
}


/* ************************ */
void MakeP_Pm_N3( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qs,
                  DBL_TYPE *Qms, DBL_TYPE *Qm,
                  DBL_TYPE *P,  DBL_TYPE *Ps,
                  DBL_TYPE *Pms, DBL_TYPE *Pm,
                  int **etaN) {

  int d;//left base of rightmost base pair.
  int pf_ij = pf_index( i, j, seqlength);
  int pf_id1, pf_dj;
  DBL_TYPE pr;

  DBL_TYPE extraTerms;

  for( d = i; d <= j - 1; d++) {
    pf_id1 = pf_index(i,d-1,seqlength);
    pf_dj = pf_index(d,j, seqlength);

    if( etaN[ EtaNIndex(d-0.5, d-0.5, seqlength)][0] == 0 || d == i ) {

      if( Q[ pf_ij] > 0) {
        pr = P[ pf_ij]*Q[ pf_id1] *
          Qs[ pf_dj]/Q[ pf_ij];

#ifdef NUPACK_SAMPLE
        if(!nupack_sample) {
#endif // NUPACK_SAMPLE
          P[pf_id1] += pr;
          Ps[pf_dj] += pr;
#ifdef NUPACK_SAMPLE
        } else {
          samplingGlobals.pq += pr;
          if(!samplingGlobals.pqset && samplingGlobals.pq >= samplingGlobals.Zq) {
            samplingGlobals.pqset = 1;
            P[pf_id1] = 1;
            Ps[pf_dj] = 1;
          }
        }
#endif // NUPACK_SAMPLE
        if( (pr > 1.0 + NUM_PRECISION ) ) {
          printf("Numerical precision loss in pairsPr.c\n");
          printf("9\n");
          exit(1);
        }
      }

      extraTerms = ExplDangle( i, d-1, seq, seqlength) *
        EXP_FUNC( -(ALPHA_3)*(d-i)/(kB*TEMP_K));


      if( etaN[ EtaNIndex( d-0.5, d-0.5, seqlength)][0] == 0) {
        //otherwise Qm not possible

        if( etaN[ EtaNIndex(i+0.5, d-0.5, seqlength)][0] == 0 ) {
          if( Qm[ pf_ij] > 0) {
            pr = Pm[ pf_ij]*Qms[ pf_dj] *
              extraTerms/Qm[pf_ij]; //Single Pair
#ifdef NUPACK_SAMPLE
            if( !nupack_sample) {
#endif // NUPACK_SAMPLE
              Pms[ pf_dj] += pr;
#ifdef NUPACK_SAMPLE
            } else {
              samplingGlobals.pqm += pr;
              if(!samplingGlobals.pqmset && samplingGlobals.pqm >= samplingGlobals.Zqm) {
                samplingGlobals.pqmset = 1;
                Pms[pf_dj] = 1;
              }
            }
#endif
            if( (pr > 1.0 + NUM_PRECISION ) ) {
              printf("Numerical precision loss in pairsPr.c\n");
              printf("10\n");
              exit(1);
            }
          }
        }

        if( d >= i+2) {
          if( Qm[ pf_ij] > 0 ) {
            pr = Pm[ pf_ij]*Qm[ pf_id1 ] *
              Qms[ pf_dj ]/Qm[ pf_ij];
           
#ifdef NUPACK_SAMPLE
            if( !nupack_sample) {
#endif // NUPACK_SAMPLE
              Pm[ pf_id1] += pr;
              Pms[ pf_dj] += pr;
#ifdef NUPACK_SAMPLE
            } else {
              samplingGlobals.pqm += pr;
              if(!samplingGlobals.pqmset && samplingGlobals.pqm >= samplingGlobals.Zqm) {
                samplingGlobals.pqmset = 1;
                Pm[pf_id1] = 1;
                Pms[pf_dj] = 1;
              }
            }
#endif // NUPACK_SAMPLE

            if( (pr > 1.0 + NUM_PRECISION ) ) {
              printf("Numerical precision loss in pairsPr.c\n");
              printf("11\n");
              exit(1);
            }

          }
        }
      }
    }
  }
}

/* ********** */
/* Ps, Pms  Recursion */
void MakePs_Pms( int i, int j, int seq[], int seqlength,
                 DBL_TYPE *Qs, DBL_TYPE *Qms, DBL_TYPE *Qb,
                 DBL_TYPE *Ps, DBL_TYPE *Pms, DBL_TYPE *Pb,
                 int *nicks, int **etaN) {

  int d; //rightmost base pair is i,d
  DBL_TYPE bp_penalty = 0.0;
  int pf_ij = pf_index( i, j, seqlength);
  int pf_id;

  DBL_TYPE extraTerms;
  DBL_TYPE pr;

  int nNicks;
  int index_ij = EtaNIndex( i+0.5, j-0.5, seqlength);
  int start;

  nNicks = etaN[ index_ij][0];
  if( nNicks >= 1) {
    start = nicks[ etaN[ index_ij][1] + nNicks - 1]+1;
  }
  else {
    start = i+4;
  }

  for( d = start; d <= j; d++) {
    pf_id = pf_index(i,d,seqlength);

    bp_penalty = 0.0;

    if( CanPair( seq[i], seq[ d]) == TRUE &&
        CanWCPair(seq[i], seq[d])) {

      if( seq[i] != BASE_C && seq[d] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }

      extraTerms = EXP_FUNC( -(NickDangle( d+1,j,nicks, etaN,
                                       FALSE, seq,seqlength) +
                           bp_penalty)/(kB*TEMP_K) );

      if( Qs[ pf_ij] > 0) {
        pr = Ps[ pf_ij] * Qb[ pf_id ] *
          extraTerms/Qs[ pf_ij];

        //printf("%i %i %i %Le %Le %Le\n", i, j, d, (long double) pr, (long double) samplingGlobals.pqs, (long double) samplingGlobals.Zqs);
        
        #ifdef NUPACK_SAMPLE
        if(!nupack_sample) {
        #endif // NUPACK_SAMPLE
          Pb[ pf_id] += pr;
        #ifdef NUPACK_SAMPLE
        } else {
          samplingGlobals.pqs += pr;
          if(!samplingGlobals.pqsset && samplingGlobals.pqs >= samplingGlobals.Zqs) {
            samplingGlobals.pqsset = 1;
            Pb[pf_id] = 1;
          }
        }
        #endif // NUPACK_SAMPLE

        if( (pr > 1.0 + NUM_PRECISION ) ) {
          printf("Numerical precision loss in pairsPr.c\n");
          printf("12 %d %.16Le\n", seq[0], (long double) pr);
          exit(1);
        }

      }
      // ********************

      extraTerms =  ExplDangle( d+1, j, seq, seqlength) *
        EXP_FUNC( -(bp_penalty + ALPHA_2 + ALPHA_3*(j-d))/(kB*TEMP_K) );

      if( Qms[ pf_ij] > 0) {
        pr = Pms[ pf_ij] *Qb[ pf_id ] *
          extraTerms /Qms[ pf_ij];

#ifdef NUPACK_SAMPLE
        if(!nupack_sample) {
#endif // NUPACK_SAMPLE
          Pb[ pf_id] += pr;
#ifdef NUPACK_SAMPLE          
        } else {
          samplingGlobals.pqms += pr;
          if(!samplingGlobals.pqmsset && samplingGlobals.pqms >= samplingGlobals.Zqms) {
            samplingGlobals.pqmsset = 1;
            Pb[pf_id] = 1;
          }
        }
#endif // NUPACK_SAMPLE

        if( (pr > 1.0 + NUM_PRECISION ) ) {
          printf("Numerical precision loss in pairsPr.c\n");
          printf("13\n");
          exit(1);
        }
      }
    }
  }
}

/* ***************************** */

void prExterior_N3( int i,int j, int seq[], int seqlength,
                    DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE * Qb_bonus,
                    DBL_TYPE *P, DBL_TYPE *Pb,
                    int *nicks, int **etaN) {

  DBL_TYPE pr;
  DBL_TYPE bp_penalty = 0.0;

  int multiNick;
  int n; //n is nick under consideration

  DBL_TYPE extraTerms;
  int pf_ij = pf_index(i,j,seqlength);
  int pf_i1m, pf_m1j1;

  int index_ij = EtaNIndex(i+0.5, j-0.5, seqlength);
  int nNicks, leftIndex;

  int iNicked, jNicked;

  if( Qb[ pf_ij] <= 0) return;

  iNicked = jNicked = FALSE;
  if( etaN[ EtaNIndex(j-0.5, j-0.5, seqlength)][0] != 0) {
    jNicked = TRUE;
  }
  if( etaN[ EtaNIndex(i+0.5, i+0.5, seqlength)][0] != 0) {
    iNicked = TRUE;
  }

  if( CanWCPair(seq[i], seq[j])) {
    bp_penalty = 0.0;

    if( seq[i] != BASE_C  && seq[j] != BASE_C) {
      bp_penalty = AT_PENALTY;
    }

    nNicks = etaN[ index_ij][0];
    leftIndex = etaN[ index_ij ][1];

    //treat each nick as rightmost nick in paired interval
    //All non-nicked cases handled outside this loop
    for( n = 0; n <= nNicks-1; n++) {
      multiNick = nicks[ leftIndex + n];

      extraTerms =
        EXP_FUNC( -1*(bp_penalty)/(kB*TEMP_K)) * Qb_bonus[pf_ij];

      if( (iNicked == FALSE && jNicked == FALSE) ||
          (i == j - 1) ||
          (multiNick == i && jNicked == FALSE) ||
          (multiNick == j-1 && iNicked == FALSE ) ) {

        pf_i1m = pf_index( i+1, multiNick, seqlength);
        pf_m1j1 = pf_index( multiNick+1, j-1, seqlength);


        if( Qb[ pf_ij] > 0) {
          pr = Pb[ pf_ij] * Q[ pf_i1m]*
            Q[ pf_m1j1] * extraTerms / Qb[ pf_ij];

        // BW K
#ifdef NUPACK_SAMPLE
          if(!nupack_sample) {
#endif // NUPACK_SAMPLE
            P[ pf_i1m] += pr;
            P[ pf_m1j1] += pr;
#ifdef NUPACK_SAMPLE
          } else {
            samplingGlobals.pqb += pr;
            if(!samplingGlobals.pqbset && samplingGlobals.pqb >= samplingGlobals.Zqb) {
              samplingGlobals.pqbset = 1;
              P[pf_i1m] = 1;
              P[pf_m1j1] = 1;
            }
          }
#endif // NUPACK_SAMPLE

          if( (pr > 1.0 + NUM_PRECISION ) ) {
            printf("Numerical precision loss in pairsPr.c\n");
            printf("2 - %d %d %d %Le %Le %Le %Le %Le\n",i, multiNick, j,
                   (long double) Pb[ pf_ij], (long double) Q[ pf_i1m],
                   (long double) Q[ pf_m1j1], (long double) extraTerms,
                   (long double) Qb[ pf_ij]) ;
            exit(1);
          }
        }

      }

    }
  }

}

/* *********************************** */

void prMultiBp_N3( int i, int j, int seq[], int seqlength,
                   DBL_TYPE *Qb, DBL_TYPE *Qms, DBL_TYPE *Qm, DBL_TYPE * Qb_bonus,
                   DBL_TYPE *Pb, DBL_TYPE *Pms, DBL_TYPE *Pm, int **etaN){

  // Decomposes the region inside pair i,j into multiloops, i.e.
  // and excludes the possibility of "top level" pseudoknots

  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE extraTerms, pr;

  int d; // d is the left base of a rightmost paired base between i, j.
  int pf_ij = pf_index(i,j,seqlength);
  int pf_i1d1, pf_dj1;

  if( Qb[ pf_ij] <= 0) return;

  if( CanWCPair(seq[i], seq[j])) {
    for( d = i+3; d <= j - 2; d++) {
      pf_i1d1 = pf_index( i+1, d-1,seqlength);
      pf_dj1 = pf_index( d, j-1, seqlength);

      //reset loop parameters
      bp_penalty = 0.0;

      if( etaN[ EtaNIndex( d-0.5, d-0.5, seqlength)][0] == 0 ) {

        if( seq[i] != BASE_C  && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }

        extraTerms = EXP_FUNC( -( ALPHA_1 + ALPHA_2 + bp_penalty) /
                           (kB*TEMP_K) ) * Qb_bonus[pf_ij];

        pr = Pb[ pf_ij] * Qm[ pf_i1d1] *
          Qms[ pf_dj1] * extraTerms / Qb[ pf_ij];

#ifdef NUPACK_SAMPLE
        if(!nupack_sample) {
#endif // NUPACK_SAMPLE
          Pm[ pf_i1d1] += pr;
          Pms[ pf_dj1] += pr;
#ifdef NUPACK_SAMPLE
        } else {
          samplingGlobals.pqb += pr;
          if(!samplingGlobals.pqbset && samplingGlobals.pqb >= samplingGlobals.Zqb) {
            samplingGlobals.pqbset = 1;
            Pm[pf_i1d1] = 1;
            Pms[pf_dj1] = 1;
          }
        }
#endif // NUPACK_SAMPLE

        if( (pr > 1.0 + NUM_PRECISION ) ) {
          printf("Numerical precision loss in pairsPr.c\n");
          printf("16\n");
          exit(1);
        }
      }
    }
  }

  return;
}

/* *************************************** */
#ifdef O_N8

void calculatePairsN8( DBL_TYPE Q[], DBL_TYPE Qb[], DBL_TYPE Qm[],
                       DBL_TYPE Qp[], DBL_TYPE Qz[], DBL_TYPE Qg[],
                       DBL_TYPE P[],
                       DBL_TYPE Pb[], DBL_TYPE Pp[], DBL_TYPE Pz[],
                       DBL_TYPE Pg[], DBL_TYPE Pbg[],
                       DBL_TYPE Pm[],
                       int seqlength,
                       char seq[]) {

  int L, i, j, d, e;
  DBL_TYPE rowsum;
  int pf_ij, gap_idej;
  //DBL_TYPE prob; //probability

#ifdef FILEOUTPUT
  FILE *fp;
#endif

  int indI;

  for( L = seqlength; L >= 1; L--) {
    for( i = 0; i <= seqlength - L; i++) {
      j = i + L - 1;
      pf_ij = pf_index(i, j, seqlength);

      MakeP_Pm_Pz( i, j, seq, seqlength,
                   Q, Qm, Qz, Qb, Qp, P, Pm, Pz, Pb, Pp);
      //requires Pp( i,j) is set.


      if( L >= 6) {
        PseudoknotLoop( i, j, pf_ij, Qg, Qz, Qp, Pp, Pz, Pg, Pbg, seq,
                        seqlength);
      }


      //require Qg
      MakePg_N8( i, j, seq, seqlength, Qg, Qm,
                 Pg, Pm);


      //require Pb
      MakePb_N5(  i, j, seq, seqlength, Qm, Qb, Qp, Pm, Pb, Pp);

      //Calculate Pbg
      if( CanPair( seq[i], seq[j]) == TRUE && L >= 7) {

        for( d = i+1; d<= j-5; d++) {
          for( e = d+4; e <= j-1; e++) {
            if(  CanPair( seq[d], seq[e]) == TRUE) {
              gap_idej = gap_index(i,d,e,j,seqlength);
              Pbg[ pf_ij] += Pg[ gap_idej];

            }
          }
        }
      }


    }
  }

#ifdef FILEOUTPUT
  fp = fopen( "Pb_N8.txt", "w");
  for( i = 0; i <= seqlength - 1; i++) {
    for( j = 0; j <= seqlength - 1; j++) {
      if ( i < j) {
        pf_ij = pf_index(i,j, seqlength);

        prob = Pb[ pf_ij] +
          Pbg[ pf_ij];
        if( prob > 1) {
          printf("Error!!! %d %d %Le %Le\n", i, j,
                 (long double) Pb[ pf_ij], (long double) Pbg[ pf_ij]);
        }
        fprintf( fp, "%d %d %.10f %.10f %.10f\n", i+1, j+1,
                 prob, Pb[ pf_ij],Pbg[ pf_ij]);

      }
      else {
        fprintf( fp, "%d %d %.10f %.10f %.10f\n", i+1, j+1,
                 Pb[ pf_index( j, i, seqlength)] +
                 Pbg[ pf_index( j, i, seqlength)],
                 Pb[ pf_index( j, i, seqlength)],
                 Pbg[ pf_index( j, i, seqlength)]);
      }
    }
  }
  fclose( fp);
#endif

  //check rowsums:
  for( i = 0; i <= seqlength - 1; i++) {
    rowsum = 0;
    for ( j = 0; j <= seqlength - 1; j++) {
      if( i <= j) {
        rowsum += Pb[  pf_index(i, j, seqlength)] +
          Pbg[ pf_index(i,j,seqlength)];
      }
      else {
        rowsum += Pb[  pf_index(j, i, seqlength)] +
          Pbg[ pf_index(j,i,seqlength)];
      }
    }
    //printf( "rs: %d %Le\n", i, rowsum);
    if( rowsum > 1) {
      printf("Error!!!!\n");
      exit(1);
    }
  }
  //  printf("Row sums < 1\n");


  //store values in pairPr
  for( i = 0; i <= seqlength-1; i++) {
    rowsum = 0;
    indI = i*(seqlength+1);
    for( j = 0; j <= seqlength; j++) {
      if( j == seqlength) {
        pairPr[ indI + j] = 1 - rowsum;
      }
      else if ( i <= j) {
        pf_ij = pf_index(i,j, seqlength);
        rowsum += Pb[ pf_ij] + Pbg[ pf_ij];
        pairPr[ indI + j] = Pb[ pf_ij] + Pbg[ pf_ij];
      }
      else {
        pf_ij = pf_index(j,i, seqlength);
        rowsum += Pb[ pf_ij] + Pbg[ pf_ij];
        pairPr[ indI + j] = Pb[ pf_ij] + Pbg[ pf_ij];
      }
    }
  }

}

/* *********************************************** */


void PseudoknotLoop( int i, int j, int pf_ij,
                     DBL_TYPE Qg[], DBL_TYPE Qz[], DBL_TYPE Qp[],
                     DBL_TYPE Pp[], DBL_TYPE Pz[], DBL_TYPE Pg[], DBL_TYPE Pbg[],
                     char seq[], int seqlength) {

  int a,b,c,d,e,f;
  DBL_TYPE bp_penalty;
  DBL_TYPE prob;

  if( j - i <= 4) return;

  //first case is exactly 1 pair per Qg
  a = i;
  f = j;
  for( b = a+1;  b <= j-4; b++) {
    if( CanPair( seq[ b], seq[ j]) == TRUE
        && CanWCPair(seq[b], seq[j])) {
      c = b;
      for( d = MAX(c+1,a+4); d <= j - 1; d++) {
        if( CanPair( seq[ a], seq[ d]) == TRUE
            && CanWCPair(seq[a], seq[d])) {
          e = d;

          bp_penalty = 0.0;
          if( seq[a] != BASE_C && seq[d] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          if( seq[c] != BASE_C && seq[f] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          if( seq[i] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          if( seq[b] != BASE_C && seq[j] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          prob = Pp[ pf_ij] *
            Qg[ gap_index( i, a, d, e, seqlength) ] *
            Qg[ gap_index( b, c, f, j, seqlength) ] *
            EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                  (kB*TEMP_K) ) *
            Qz[ pf_index( e + 1, f - 1, seqlength)]*
            Qz[ pf_index( c + 1, d - 1, seqlength)]*
            Qz[ pf_index( a + 1, b - 1, seqlength)]/
            Qp[ pf_ij];

          Pz[ pf_index( a+1,b-1,seqlength)] += prob;
          Pz[ pf_index( c+1,d-1,seqlength)] += prob;
          Pz[ pf_index( e+1,f-1,seqlength)] += prob;
          Pg[ gap_index(i,a,d,e,seqlength)] += prob;
          Pg[ gap_index(b,c,f,j,seqlength)] += prob;
          Pbg[ pf_index(a, d, seqlength)] += prob;
          Pbg[ pf_index(c, f, seqlength)] += prob;

          if( prob > 1) {
            printf("%d %d %d %d %d %d %d %d pk prob = %Le\n",
                   i,a,b,c,d,e,f,j,prob);
            exit(1);
          }


        }
      }
    }
  }

  if( j - i <= 6) return;

  // a = i, f != j
  //left Qg has exactly 1 pair, the other has 2+
  a = i;
  for( b = a+1; b <= j - 6; b++) {
    if( CanPair( seq[ b], seq[ j]) == TRUE
        && CanWCPair(seq[b], seq[j])) {
      for( c = b+1; c <= j - 5; c++) {
        for( d = MAX(c+1,a+4); d <= j - 2; d++) {
          if( CanPair( seq[ a], seq[ d]) == TRUE
              && CanWCPair(seq[a], seq[d])) {
            e = d;
            if( CanPair( seq[ i], seq[ e]) == TRUE
                && CanWCPair(seq[i], seq[e])) {
              for( f = MAX(e+1,c+4); f <= j-1; f++) {
                if( CanPair( seq[c], seq[f]) == TRUE
                    && CanWCPair(seq[c], seq[f])) {

                  bp_penalty = 0.0;
                  if( seq[a] != BASE_C && seq[d] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }
                  if( seq[c] != BASE_C && seq[f] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }
                  if( seq[i] != BASE_C && seq[e] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }
                  if( seq[b] != BASE_C && seq[j] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }

                  prob = Pp[ pf_ij] *
                    Qg[ gap_index( i, a, d, e, seqlength) ] *
                    Qg[ gap_index( b, c, f, j, seqlength) ] *
                    EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                          (kB*TEMP_K) ) *
                    Qz[ pf_index( e + 1, f - 1, seqlength)]*
                    Qz[ pf_index( c + 1, d - 1, seqlength)]*
                    Qz[ pf_index( a + 1, b - 1, seqlength)]/
                    Qp[ pf_ij];

                  Pz[ pf_index( a+1,b-1,seqlength)] += prob;
                  Pz[ pf_index( c+1,d-1,seqlength)] += prob;
                  Pz[ pf_index( e+1,f-1,seqlength)] += prob;
                  Pg[ gap_index(i,a,d,e,seqlength)] += prob;
                  Pg[ gap_index(b,c,f,j,seqlength)] += prob;
                  Pbg[ pf_index(a, d, seqlength)] += prob;
                  Pbg[ pf_index(c, f, seqlength)] += prob;

                  if( prob > 1) {
                    printf("%d %d %d %d %d %d %d %d pk prob = %Le\n",
                           i,a,b,c,d,e,f,j,prob);
                    exit(1);
                  }



                }
              }
            }
          }
        }
      }
    }
  }


  // a != i, b == c
  //right Qg has exactly 1 pair, the other has 2+
  f = j;
  for( a = i+1; a <= j - 6; a++) {
    for( b = a+1; b <= j - 4; b++) {
      if( CanPair( seq[ b], seq[ j]) == TRUE
          && CanWCPair(seq[b], seq[j])) {
        c = b;
        if( CanPair( seq[c], seq[f]) == TRUE
            && CanWCPair(seq[c], seq[f])) {
          for( d = MAX(c+1,a+4); d <= j - 2; d++) {
            if( CanPair( seq[ a], seq[ d]) == TRUE
                && CanWCPair(seq[a], seq[d])) {
              for( e = d+1; e <= j - 1; e++) {
                if( CanPair( seq[ i], seq[ e]) == TRUE
                    && CanWCPair(seq[i], seq[e])) {

                  bp_penalty = 0.0;
                  if( seq[a] != BASE_C && seq[d] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }
                  if( seq[c] != BASE_C && seq[f] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }
                  if( seq[i] != BASE_C && seq[e] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }
                  if( seq[b] != BASE_C && seq[j] != BASE_C) {
                    bp_penalty += AT_PENALTY;
                  }

                  prob = Pp[ pf_ij] *
                    Qg[ gap_index( i, a, d, e, seqlength) ] *
                    Qg[ gap_index( b, c, f, j, seqlength) ] *
                    EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                          (kB*TEMP_K) ) *
                    Qz[ pf_index( e + 1, f - 1, seqlength)]*
                    Qz[ pf_index( c + 1, d - 1, seqlength)]*
                    Qz[ pf_index( a + 1, b - 1, seqlength)]/
                    Qp[ pf_ij];

                  Pz[ pf_index( a+1,b-1,seqlength)] += prob;
                  Pz[ pf_index( c+1,d-1,seqlength)] += prob;
                  Pz[ pf_index( e+1,f-1,seqlength)] += prob;
                  Pg[ gap_index(i,a,d,e,seqlength)] += prob;
                  Pg[ gap_index(b,c,f,j,seqlength)] += prob;
                  Pbg[ pf_index(a, d, seqlength)] += prob;
                  Pbg[ pf_index(c, f, seqlength)] += prob;

                  if( prob > 1) {
                    printf("%d %d %d %d %d %d %d %d pk prob = %Le\n",
                           i,a,b,c,d,e,f,j,prob);
                    exit(1);
                  }


                }
              }
            }
          }
        }
      }
    }
  }

  if( j - i <= 7) return;

  //both Qg have 2+ pairs
  for( a = i + 1; a <= j - 7; a++) {
    for( b = a+1; b <= j - 6; b++) {
      if( CanPair( seq[ b], seq[ j]) == TRUE
          && CanWCPair(seq[b], seq[j])) {
        for( c = b + 1; c <= j - 5; c++) {
          for( d = MAX(c+1,a+4); d <= j - 3; d++) {
            if( CanPair( seq[ a], seq[ d]) == TRUE
                && CanWCPair(seq[a], seq[d])) {
              for( e = d+1; e <= j - 2; e++) {
                if( CanPair( seq[ i], seq[ e]) == TRUE
                    && CanWCPair(seq[i], seq[e])) {
                  for( f = MAX(e+1,c+4); f <= j-1; f++) {
                    if( CanPair( seq[c], seq[f]) == TRUE
                        && CanWCPair(seq[c], seq[f])) {

                      bp_penalty = 0.0;
                      if( seq[a] != BASE_C && seq[d] != BASE_C) {
                        bp_penalty += AT_PENALTY;
                      }
                      if( seq[c] != BASE_C && seq[f] != BASE_C) {
                        bp_penalty += AT_PENALTY;
                      }
                      if( seq[i] != BASE_C && seq[e] != BASE_C) {
                        bp_penalty += AT_PENALTY;
                      }
                      if( seq[b] != BASE_C && seq[j] != BASE_C) {
                        bp_penalty += AT_PENALTY;
                      }

                      prob = Pp[ pf_ij] *
                        Qg[ gap_index( i, a, d, e, seqlength) ] *
                        Qg[ gap_index( b, c, f, j, seqlength) ] *
                        EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                              (kB*TEMP_K) ) *
                        Qz[ pf_index( e + 1, f - 1, seqlength)]*
                        Qz[ pf_index( c + 1, d - 1, seqlength)]*
                        Qz[ pf_index( a + 1, b - 1, seqlength)]/
                        Qp[ pf_ij];

                      Pz[ pf_index( a+1,b-1,seqlength)] += prob;
                      Pz[ pf_index( c+1,d-1,seqlength)] += prob;
                      Pz[ pf_index( e+1,f-1,seqlength)] += prob;
                      Pg[ gap_index(i,a,d,e,seqlength)] += prob;
                      Pg[ gap_index(b,c,f,j,seqlength)] += prob;
                      Pbg[ pf_index(a, d, seqlength)] += prob;
                      Pbg[ pf_index(c, f, seqlength)] += prob;

                      if( prob > 1) {
                        printf("%d %d %d %d %d %d %d %d pk prob = %Le\n",
                               i,a,b,c,d,e,f,j,prob);
                        exit(1);
                      }

                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
#endif //O_N8
/* ******************************************** */

void calculatePairsN5( DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE *Qm,
                       DBL_TYPE *Qp, DBL_TYPE *Qz, DBL_TYPE *Qg,
                       DBL_TYPE *Qgl, DBL_TYPE *Qgr, DBL_TYPE *Qgls,
                       DBL_TYPE *Qgrs,
                       DBL_TYPE **QgIx, DBL_TYPE **QgIx_1, DBL_TYPE **QgIx_2,
                       DBL_TYPE *P,
                       DBL_TYPE *Pb, DBL_TYPE *Pp, DBL_TYPE *Pz,
                       DBL_TYPE *Pg, DBL_TYPE *Pbg,
                       DBL_TYPE *Pm, DBL_TYPE *Pgl, DBL_TYPE *Pgr,
                       DBL_TYPE *Pgls, DBL_TYPE *Pgrs,
                       int seqlength,
                       int seq[]) {

  int L, i, j, d, e;
  DBL_TYPE rowsum, rowsumPb, rowsumPbg;
  int pf_ij, gap_idej;

  float *preX, *preX_1, *preX_2;

  int indI;
  extern DBL_TYPE *pairPrPb;
  extern DBL_TYPE *pairPrPbg;

  DBL_TYPE *PgIx, *PgIx_1, *PgIx_2;
  PgIx = PgIx_1 = PgIx_2 = NULL;
  preX = preX_1 = preX_2 = NULL;

  for( L = seqlength; L >= 1; L--) {

    prManageQgIx(  QgIx, QgIx_1, QgIx_2, &PgIx, &PgIx_1, &PgIx_2,
                   &preX, &preX_1, &preX_2, L-1,
                   seqlength);

    for( i = 0; i <= seqlength - L; i++) {

      j = i + L - 1;
      pf_ij = pf_index(i, j, seqlength);

      MakeP_Pm_Pz( i, j, seq, seqlength,
                   Q, Qm, Qz, Qb, Qp, P, Pm, Pz, Pb, Pp);

      //requires Pp( i,j) is set.

      if( L >= 6) {
        PseudoknotLoopN5( i, j, pf_ij, Qp, Qgl, Qgr, Qg, Qz,
                          Pp, Pgl, Pgr, Pg, Pz, Pbg, seq, seqlength);
      }

      //requires Pgr(i,j) is set
      MakePgr(i, j, seq, seqlength, Qgr, Qgl, Qz, Pgr, Pgl, Pz);
      //require Pgl set

      MakePgl(i, j, seq, seqlength, Qg, Qgl, Qz, Pg, Pgl, Pz, Pbg);

      //require Pgrs set
      MakePgrs( i, j, seq, seqlength, Qg, Qm, Qgrs, Pg, Pm, Pgrs);

      //require Pgls
      MakePgls( i, j, seq, seqlength, Qg, Qm, Qgls, Pg, Pm, Pgls);


      //require Qg
      MakePg_N5( i, j, seq, seqlength, Qg, Qm, Qgls, Qgrs, *QgIx,
                 *QgIx_2, Pg, Pm, Pgls, Pgrs, PgIx, PgIx_2, preX, preX_2);

      //require Pb
      MakePb_N5(  i, j, seq, seqlength, Qm, Qb, Qp, Pm, Pb, Pp);

      //Calculate Pbg
      if( CanPair( seq[i], seq[j]) == TRUE &&  L >= 7) {

        for( d = i+1; d<= j-5; d++) {
          for( e = d+4; e <= j-1; e++) {
            if(  CanPair( seq[d], seq[e]) == TRUE) {
              gap_idej = gap_index(i,d,e,j,seqlength);
              Pbg[ pf_ij] += Pg[ gap_idej];

            }
          }
        }
      }

    }
  }

  //check rowsums:
  for( i = 0; i <= seqlength - 1; i++) {
    rowsum = 0;
    for ( j = 0; j <= seqlength - 1; j++) {
      if( i <= j) {
        rowsum += Pb[  pf_index(i, j, seqlength)] +
          Pbg[ pf_index(i,j,seqlength)];
      }
      else {
        rowsum += Pb[  pf_index(j, i, seqlength)] +
          Pbg[ pf_index(j,i,seqlength)];
      }
    }
    //printf( "rs: %d %Le\n", i, rowsum);
    if( rowsum > 1) {
      printf("Error!!!!\n");
      exit(1);
    }
  }
  //  printf("Row sums < 1\n");

  free( PgIx); PgIx = NULL;
  free( PgIx_1); PgIx_1 = NULL;
  free( PgIx_2); PgIx_2 = NULL;

  free( preX);
  free( preX_1);
  free( preX_2);
  preX = preX_1 = preX_2 = NULL;

  //store values in pairPr (pairPrPbg, PairPrPb)
  for( i = 0; i <= seqlength-1; i++) {
    rowsum = rowsumPb = rowsumPbg = 0;

    indI = i*(seqlength+1);
    for( j = 0; j <= seqlength; j++) {
      if( j == seqlength) {
        pairPr[ indI + j] = 1 - rowsum;
        pairPrPb[ indI + j] = 1 - rowsumPb;
        pairPrPbg[ indI + j] = 1 - rowsumPbg;
      }
      else if ( i <= j) {
        pf_ij = pf_index(i,j, seqlength);
        rowsum += Pb[ pf_ij] + Pbg[ pf_ij];
        rowsumPb += Pb[ pf_ij];
        rowsumPbg += Pbg[ pf_ij];
        pairPr[ indI + j] = Pb[ pf_ij] + Pbg[ pf_ij];
        pairPrPb[ indI + j] = Pb[ pf_ij];
        pairPrPbg[ indI + j] = Pbg[ pf_ij];
      }
      else {
        pf_ij = pf_index(j,i, seqlength);
        rowsum += Pb[ pf_ij] + Pbg[ pf_ij];
        rowsumPb += Pb[ pf_ij];
        rowsumPbg += Pbg[ pf_ij];
        pairPr[ indI + j] = Pb[ pf_ij] + Pbg[ pf_ij];
        pairPrPb[ indI + j] = Pb[ pf_ij];
        pairPrPbg[ indI + j] = Pbg[ pf_ij];
      }
    }
  }



}

/* ******************************************** */

void PseudoknotLoopN5( int i, int j, int pf_ij,
                       DBL_TYPE *Qp, DBL_TYPE *Qgl, DBL_TYPE *Qgr, DBL_TYPE *Qg,
                       DBL_TYPE *Qz,
                       DBL_TYPE *Pp, DBL_TYPE *Pgl, DBL_TYPE *Pgr, DBL_TYPE *Pg,
                       DBL_TYPE *Pz, DBL_TYPE *Pbg,
                       int seq[], int seqlength) {

  int d,e,f;
  DBL_TYPE bp_penalty, new_bp_penalty;
  DBL_TYPE prob;

  int a, b, c;

  if( j - i <= 4) return;
  //case 1 both Qg are exactly 1 pair
  //first case is exactly 1 pair per Qg

  a = i;
  f = j;
  for( b = a+1;  b <= j-4; b++) {
    if( CanPair( seq[ b], seq[ j]) == TRUE
        && CanWCPair(seq[b], seq[j])) {
      c = b;
      for( d = MAX(c+1,a+4); d <= j - 1; d++) {
        if( CanPair( seq[ a], seq[ d]) == TRUE
            && CanWCPair(seq[a], seq[d])) {
          e = d;

          bp_penalty = 0.0;
          if( seq[a] != BASE_C && seq[d] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          if( seq[c] != BASE_C && seq[f] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          if( seq[i] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          if( seq[b] != BASE_C && seq[j] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          prob = Pp[ pf_ij] *
            Qg[ gap_index( i, a, d, e, seqlength) ] *
            Qg[ gap_index( b, c, f, j, seqlength) ] *
            EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                  (kB*TEMP_K) ) *
            Qz[ pf_index( e + 1, f - 1, seqlength)]*
            Qz[ pf_index( c + 1, d - 1, seqlength)]*
            Qz[ pf_index( a + 1, b - 1, seqlength)] / Qp[ pf_ij];


          if( prob > 1) {
            printf("pk5 1 -  %d %d %Le\n", i, j, (long double) prob);
          }

          Pz[ pf_index( a+1,b-1,seqlength)] += prob;
          Pz[ pf_index( c+1,d-1,seqlength)] += prob;
          Pz[ pf_index( e+1,f-1,seqlength)] += prob;
          Pg[ gap_index(i,a,d,e,seqlength)] += prob;
          Pg[ gap_index(b,c,f,j,seqlength)] += prob;
          Pbg[ pf_index(a, d, seqlength)] += prob;
          Pbg[ pf_index(c, f, seqlength)] += prob;


        }
      }
    }
  }


  if( j - i <= 6) return;

  //case 2 left Qg is exactly 1 pair, right is 2+
  for( d = i+1; d <= j-6; d++) {
    if( CanPair( seq[d], seq[j]) == TRUE &&
        CanWCPair(seq[d], seq[j])
        ) {
      bp_penalty = 0;
      if( seq[d] != BASE_C && seq[j] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }

      for( e = MAX(d+2, i+4); e <= j-2; e++) {
        f = e;
        if( CanPair( seq[i], seq[f]) == TRUE &&
            CanWCPair(seq[i], seq[f])
            ) {
          new_bp_penalty = bp_penalty;
          if( seq[i] != BASE_C && seq[f] != BASE_C) {
            new_bp_penalty += 2*AT_PENALTY; //count twice
          }

          prob = Pp[ pf_ij] *Qg[ gap_index(i, i, e, f, seqlength)] *
            Qz[ pf_index( i+1, d-1, seqlength)] *
            Qgr[ gap_index(d, e-1, f+1, j, seqlength)] *
            EXP_FUNC( -(new_bp_penalty + BETA_2)/(kB*TEMP_K) ) / Qp[ pf_ij];


          if( prob > 1) {
            printf("pk5 2 -  %d %d %Le\n", i, j, (long double) prob);
          }

          Pg[ gap_index(i,i,e,f,seqlength)] += prob;
          Pz[ pf_index( i+1,d-1,seqlength)] += prob;
          Pgr[ gap_index( d,e-1,f+1,j,seqlength)] += prob;
          Pbg[ pf_index( i,e,seqlength)] += prob;

        }
      }
    }
  }


  //case 2 left Qg is 2+ pairs, right is 1
  for( d = i + 2; d <= j-4; d++) {
    if( CanPair( seq[d], seq[j]) == TRUE &&
        CanWCPair(seq[d], seq[j])
        ) {
      bp_penalty = 0;
      if( seq[d] != BASE_C && seq[j] != BASE_C) {
        bp_penalty = 2*AT_PENALTY; //count twice
      }

      for( e = MAX(d+1, i+4); e <= j-2; e++) {
        for( f = e+1; f <= j-1; f++) {

          if( CanPair( seq[i], seq[f]) == TRUE &&
              CanWCPair(seq[i], seq[f])
              ) {
            new_bp_penalty = bp_penalty;
            if( seq[i] != BASE_C && seq[f] != BASE_C) {
              new_bp_penalty += AT_PENALTY;
            }

            prob = Pp[ pf_ij]*Qgl[ gap_index(i, d-1, e, f, seqlength)] *
              Qg[ gap_index(d, d, j, j, seqlength)] *
              Qz[ pf_index( d+1, e-1, seqlength)] *
              Qz[ pf_index( f+1,j-1, seqlength)] *
              EXP_FUNC( -(new_bp_penalty+BETA_2)/(kB*TEMP_K) )/ Qp[ pf_ij];


            if( prob > 1) {
              printf("pk5 3 -  %d %d %Le\n", i, j,(long double) prob);
            }

            Pg[ gap_index(d, d, j, j, seqlength)] += prob;
            Pgl[ gap_index(i,d-1,e,f,seqlength)] += prob;
            Pz[ pf_index( d+1, e-1, seqlength)] += prob;
            Pz[ pf_index( f+1,j-1, seqlength)] += prob;
            Pbg[ pf_index(d,j,seqlength)] += prob;

          }
        }
      }
    }
  }



  if( j - i <= 7) return;

  for( d = i + 2; d <= j-4; d++) {
    if( CanPair( seq[d], seq[j]) == TRUE &&
        CanWCPair(seq[d], seq[j])
        ) {
      bp_penalty = 0;
      if( seq[d] != BASE_C && seq[j] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }

      for( e = MAX(d+2, i+5); e <= j-3; e++) {
        for( f = e+1; f <= j-2; f++) {

          if( CanPair( seq[i], seq[f]) == TRUE &&
              CanWCPair(seq[i], seq[f])
              ) {
            new_bp_penalty = bp_penalty;
            if( seq[i] != BASE_C && seq[f] != BASE_C) {
              new_bp_penalty += AT_PENALTY;
            }

            if( Pp[ pf_ij] > 0) {
              prob = Pp[ pf_ij]*Qgl[ gap_index(i, d-1, e, f, seqlength)] *
                Qgr[ gap_index(d, e-1, f+1, j, seqlength)] *
                EXP_FUNC( -(new_bp_penalty)/(kB*TEMP_K) ) /
                Qp[ pf_ij];

              if( prob > 1) {
                printf("a %d %d %Le\n", i, j, (long double) prob);
              }

              Pgl[ gap_index(i,d-1,e,f,seqlength)] += prob;
              Pgr[ gap_index(d,e-1,f+1,j,seqlength)] += prob;
            }
          }
        }
      }
    }
  }

}

/* ******************************************** */

void MakeP_Pm_Pz( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qm, DBL_TYPE *Qz,
                  DBL_TYPE *Qb, DBL_TYPE *Qp, DBL_TYPE *P,
                  DBL_TYPE *Pm, DBL_TYPE *Pz, DBL_TYPE *Pb,
                  DBL_TYPE *Pp) {

  int d, e; // d - e is internal basepair or pk boundary
  DBL_TYPE bp_penalty;
  int pf_ij = pf_index(i, j, seqlength);
  int pf_de, pf_id1;
  DBL_TYPE prob;

  for( d = i; d <= j - 4; d++) {
    pf_id1 = pf_index( i, d-1, seqlength);
    for( e = d + 4; e <= j; e++) {
      pf_de = pf_index( d, e, seqlength);
      if( CanPair( seq[d], seq[e]) == TRUE &&
          CanWCPair(seq[d], seq[e]) ) {
        bp_penalty = 0;
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }

        if( P[ pf_ij] > 0) {
          prob = P[ pf_ij]*Q[ pf_id1] *
            Qb[ pf_de ] *
            EXP_FUNC( -bp_penalty/(kB*TEMP_K) ) *
            ExplDangle( e+1, j, seq, seqlength) / Q[ pf_ij];
          if( prob > 1) {
            printf("b %d %d %Le\n", i, j, (long double) prob);
          }

          P[ pf_id1] += prob;
          Pb[ pf_de] += prob;

        }

        if( i != 0 && j != seqlength - 1) {
          if( Pm[ pf_ij] > 0) {
            prob = Pm[ pf_ij] *
              EXP_FUNC( -(ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty)/
                    (kB*TEMP_K) )*
              Qb[ pf_de] *
              ExplDangle( e+1, j, seq, seqlength) *
              ExplDangle( i, d-1, seq, seqlength) / Qm[ pf_ij];
            if( prob > 1) {
              printf("c %d %d %Le\n", i, j, (long double) prob);
            }

            Pb[ pf_de] += prob;
          }

          if( d >= i+5) {
            if( Pm[ pf_ij] > 0) {
              prob = Pm[ pf_ij]*Qm[ pf_id1] *
                Qb[ pf_de] *
                EXP_FUNC( -(ALPHA_2 + ALPHA_3*(j-e) + bp_penalty)/
                      (kB*TEMP_K) )*
                ExplDangle( e+1, j, seq, seqlength) / Qm[ pf_ij];
              Pm[ pf_id1] += prob;
              Pb[ pf_de] += prob;
              if( prob > 1) {
                printf("d %d %d %Le\n", i, j, (long double) prob);
              }
            }
          }

          if( Pz[ pf_ij] > 0) {
            prob = Pz[ pf_ij] *Qz[ pf_id1]*
              Qb[ pf_de] * EXP_FUNC( -(BETA_2 + BETA_3*(j-e) + bp_penalty)/
                                 (kB*TEMP_K) )*
              ExplDangle( e+1, j, seq, seqlength)/ Qz[ pf_ij];
            Pz[ pf_id1] += prob;
            Pb[ pf_de] += prob;
            if( prob > 1) {
              printf("e %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }

      }

    }
  }

  for( d = i; d <= j - 5; d++) {
    pf_id1 = pf_index( i, d-1, seqlength);
    for( e = d + 5; e <= j; e++) {
      pf_de = pf_index( d, e, seqlength);

      if( P[ pf_ij] > 0) {
        prob = P[ pf_ij] * Q[ pf_id1] *
          Qp[ pf_de ] *
          EXP_FUNC( -BETA_1/(kB*TEMP_K) ) *
          ExplDangle( e+1, j, seq, seqlength) / Q[ pf_ij];
        P[ pf_id1] += prob;
        Pp[ pf_de] += prob;

        if( prob > 1) {
          printf("f %d %d %Le\n", i, j, (long double) prob);
        }

      }

      if( i != 0 && j != seqlength - 1) {
        if( Pm[ pf_ij] > 0) {
          prob = Pm[ pf_ij] *
            EXP_FUNC( -(BETA_1M + 2*ALPHA_2 + ALPHA_3*(d-i + j-e))/
                  (kB*TEMP_K) )*
            Qp[ pf_de] *
            ExplDangle( e+1, j, seq, seqlength) *
            ExplDangle( i, d-1, seq, seqlength) / Qm[ pf_ij];
          Pp[ pf_de] += prob;
          if( prob > 1) {
            printf("g %d %d %Le\n", i, j,(long double)  prob);
          }

        }

        if( d >= i+5) {
          if( Pm[ pf_ij] > 0) {
            prob = Pm[ pf_ij] * Qm[ pf_id1] *
              Qp[ pf_de] *
              EXP_FUNC( -(BETA_1M + 2*ALPHA_2 + ALPHA_3*(j-e))/
                    (kB*TEMP_K) )*
              ExplDangle( e+1, j, seq, seqlength) / Qm[ pf_ij];
            Pm[ pf_id1] += prob;
            Pp[ pf_de] += prob;

            if( prob > 1) {
              printf("h %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }

        if( Pz[ pf_ij] > 0) {
          prob = Pz[ pf_ij] *Qz[ pf_id1]*
            Qp[ pf_de] * EXP_FUNC( -(BETA_1P + 2*BETA_2 + BETA_3*(j-e) )/
                               (kB*TEMP_K) )*
            ExplDangle( e+1, j, seq, seqlength)/ Qz[ pf_ij];
          Pz[ pf_id1] += prob;
          Pp[ pf_de] += prob;

          if( prob > 1) {
            printf("i %d %d %Le\n", i, j, (long double) prob);
          }

        }

      }
    }
  }
}

/* ****************************** */

void MakePgr( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Qgr, DBL_TYPE *Qgl, DBL_TYPE *Qz,
              DBL_TYPE *Pgr, DBL_TYPE *Pgl, DBL_TYPE *Pz) {

  //make sure to call this BEFORE MakeQg and MakeQgl

  int d, e, f;
  int gap_idej, gap_idfj, pf_ef1;
  DBL_TYPE prob;

  for( d = i+1; d <= j-3; d++) {
    for( e = d+2; e <= j-1; e++) {
      gap_idej =  gap_index(i, d, e, j, seqlength);
      for( f = e; f <= j-1; f++) {
        gap_idfj = gap_index( i,d,f,j,seqlength);
        pf_ef1 = pf_index( e, f-1,seqlength);

        if( Pgr[ gap_idej] > 0) {
          prob = Pgr[ gap_idej] *
            Qgl[ gap_idfj]*
            Qz[ pf_ef1]/ Qgr[ gap_idej];
          Pgl[ gap_idfj] += prob;
          Pz[ pf_ef1] += prob;

          if( prob > 1) {
            printf("j %d %d %Le\n", i, j, (long double) prob);
          }

        }
      }
    }
  }
}
/* ********************************** */

void MakePgl( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Qg, DBL_TYPE *Qgl, DBL_TYPE *Qz,
              DBL_TYPE *Pg, DBL_TYPE *Pgl, DBL_TYPE *Pz,
              DBL_TYPE *Pbg) {
  //make sure to call this BEFORE MakeQg and AFTER MakeQgr

  int d, e, f;
  int gap_idfj, gap_iefj, pf_d1e, pf_df;
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE prob;

  for( d = i+1; d <= j-5; d++) {
    for( f = d+4; f <= j-1; f++) {
      if( CanPair( seq[d], seq[f]) == TRUE &&
          CanWCPair(seq[d], seq[f])
          ) {
        bp_penalty = 0.0;
        if( seq[d] != BASE_C && seq[f] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }

        pf_df = pf_index(d, f, seqlength);
        gap_idfj =  gap_index(i, d, f, j, seqlength);
        for( e = d; e <= f-2; e++) {
          gap_iefj = gap_index( i,e,f,j,seqlength);
          pf_d1e = pf_index( d+1,e,seqlength);

          if( Pgl[ gap_iefj] > 0) {
            prob = Pgl[ gap_iefj] *
              Qg[ gap_idfj]*Qz[ pf_d1e] *
              EXP_FUNC( -(BETA_2 + bp_penalty)/(kB*TEMP_K) ) /
              Qgl[ gap_iefj];
            Pg[ gap_idfj] += prob;
            Pz[ pf_d1e] += prob;
            Pbg[ pf_df] += prob;
            if( prob > 1) {
              printf("k %d %d %Le\n", i, j, (long double) prob);
            }
          }

        }
      }
    }
  }
}

/* ******************************** */

void MakePgrs( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm, DBL_TYPE *Qgrs, DBL_TYPE *Pg, DBL_TYPE *Pm,
               DBL_TYPE *Pgrs) {

  int d, e, f;
  int gap_idej, gap_idef, pf_f1j;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE prob;

  for( d = i + 1; d <= j-10; d++) {
    for( e = d+4; e <= j-6; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE) {
        gap_idej = gap_index(i,d,e,j,seqlength);
        for( f = e+1; f <= j-5; f++) {
          if( CanPair( seq[i], seq[f]) == TRUE &&
              CanWCPair(seq[i], seq[f])
              ) {

            bp_penalty = 0;
            if( seq[i] != BASE_C && seq[f] != BASE_C) {
              bp_penalty = AT_PENALTY;
            }
            pf_f1j = pf_index( f+1, j, seqlength);
            gap_idef = gap_index( i,d,e,f,seqlength);
            if( Pgrs[ gap_idej] > 0) {
              prob = Pgrs[ gap_idej] *
                Qm[ pf_f1j] *
                Qg[ gap_idef] *
                EXP_FUNC( -(ALPHA_2 + bp_penalty)/(kB*TEMP_K) )/
                Qgrs[ gap_idej];
              Pm[ pf_f1j] += prob;
              Pg[ gap_idef] += prob;
              if( prob > 1) {
                printf("l %d %d %Le\n", i, j, (long double) prob);
              }

            }
          }
        }
      }
    }
  }
}
/* *************************************** */
void MakePgls( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm,  DBL_TYPE *Qgls, DBL_TYPE *Pg, DBL_TYPE *Pm,
               DBL_TYPE *Pgls) {
  int c, d, e;
  int pf_ic1, gap_idej, gap_cdej;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE prob;

  for( c = i + 5; c <= j-6; c++) {
    if( CanPair( seq[c], seq[j]) == TRUE
        && CanWCPair(seq[c], seq[j])
        ) {

      bp_penalty = 0.0;
      if( seq[c] != BASE_C && seq[j] != BASE_C) {
        bp_penalty = AT_PENALTY;
      }

      pf_ic1 = pf_index( i, c-1, seqlength);
      for( d = c+1; d <= j-5; d++) {

        for( e = d+4; e <= j-1; e++) {
          if( CanPair( seq[d], seq[e]) == TRUE) {

            gap_idej = gap_index( i,d,e,j,seqlength);
            gap_cdej = gap_index( c,d,e,j, seqlength);

            if( Pgls[ gap_idej] > 0) {
              prob =
                Pgls[ gap_idej] *
                Qm[ pf_ic1] *
                Qg[ gap_cdej] *
                EXP_FUNC( -(ALPHA_2 + bp_penalty)/(kB*TEMP_K) )/
                Qgls[ gap_idej];
              Pm[ pf_ic1] += prob;
              Pg[ gap_cdej] += prob;

              if( prob > 1) {
                printf("m %d %d %Le\n", i, j, (long double) prob);
              }

            }
          }
        }
      }
    }
  }
}
/* ************************** */

void MakePg_N5( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
                DBL_TYPE *Qm, DBL_TYPE *Qgls, DBL_TYPE *Qgrs,
                DBL_TYPE *QgIx, DBL_TYPE *QgIx_2,
                DBL_TYPE *Pg, DBL_TYPE *Pm, DBL_TYPE *Pgls, DBL_TYPE *Pgrs,
                DBL_TYPE *PgIx, DBL_TYPE *PgIx_2,
                float *preX, float *preX_2) {

  //  Make the gap matrix for Qg
  int c,d,e,f;
  int gap_idej,  pf_i1d1, pf_e1j1, gap_i1def, gap_cdej1, pf_i1c1;
  DBL_TYPE IJ_bp_penalty = 0.0;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE prob;

  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }


  prFastILoopsN5(i, j, seq, seqlength, Qg, QgIx, QgIx_2, Pg, PgIx, PgIx_2, preX, preX_2);

  if( CanPair( seq[i], seq[j]) == TRUE &&
      CanWCPair(seq[i], seq[j])) {
    //Case 2:  Multiloop Left
    for( d = i+6; d <= j-5; d++) {
      pf_i1d1 = pf_index(i+1,d-1,seqlength);
      for( e = d+4; e <= j-1; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE
            && CanWCPair(seq[d], seq[e])) {
          gap_idej = gap_index(i,d,e,j,seqlength);

          bp_penalty = IJ_bp_penalty;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          if( Pg[ gap_idej] > 0) {
            prob =
              Pg[ gap_idej] *
              Qm[ pf_i1d1] *
              EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 +
                      (j-e-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
              ExplDangle( e+1, j-1, seq, seqlength) /
              Qg[ gap_idej];
            Pm[ pf_i1d1] += prob;
            if( prob > 1) {
              printf("n %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }
      }
    }

    //Case 3:  Multiloop Right
    for( d = i+1; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE
            && CanWCPair(seq[d], seq[e])
            ) {
          gap_idej = gap_index(i,d,e,j,seqlength);

          bp_penalty = IJ_bp_penalty;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          pf_e1j1 = pf_index( e+1,j-1,seqlength);

          if( Pg[ gap_idej] > 0) {
            prob =
              Pg[ gap_idej] *
              Qm[ pf_e1j1] *
              EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 +
                      (d-i-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
              ExplDangle( i+1, d-1, seq, seqlength) /
              Qg[ gap_idej];
            Pm[ pf_e1j1] += prob;

            if( prob > 1) {
              printf("o %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }
      }
    }

    //Case 4: Multiloop Both Sides
    for( d = i+6; d <= j-10; d++) {
      pf_i1d1 = pf_index( i+1,d-1,seqlength);
      for( e = d+4; e <= j-6; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE
            && CanWCPair(seq[d], seq[e])
            ) {
          gap_idej = gap_index(i,d,e,j,seqlength);

          bp_penalty = IJ_bp_penalty;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          pf_e1j1 = pf_index( e+1,j-1,seqlength);

          if( Pg[ gap_idej] > 0) {
            prob =
              Pg[ gap_idej] *
              Qm[ pf_i1d1] *
              Qm[ pf_e1j1] *
              EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + bp_penalty)/(kB*TEMP_K)) /
              Qg[ gap_idej];
            Pm[ pf_i1d1] += prob;
            Pm[ pf_e1j1] += prob;

            if( prob > 1) {
              printf("p %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }
      }
    }

    //Case 5: Interior loop + multi left
    for( d = i+7; d <= j-6; d++) {
      for( e = d+4; e <= j-2; e++) {
        gap_idej = gap_index(i,d,e,j,seqlength);
        if( CanPair( seq[d], seq[e] ) == TRUE) {

          for( f = e+1; f <= j-1; f++) {
            bp_penalty = IJ_bp_penalty;

            gap_i1def = gap_index( i+1,d,e,f,seqlength);
            if( Pg[ gap_idej] > 0) {
              prob =
                Pg[ gap_idej] *
                EXP_FUNC( -(ALPHA_1+ALPHA_2 + (j-f-1)*ALPHA_3+ bp_penalty)/
                      (kB*TEMP_K)) *
                Qgls[ gap_i1def] *
                ExplDangle( f+1, j-1, seq, seqlength)/
                Qg[ gap_idej];
              Pgls[ gap_i1def] += prob;
              if( prob > 1) {
                printf("q %d %d %Le\n", i, j, (long double) prob);
              }

            }
          }
        }
      }
    }

    //Case6: Interior loop + multi right
    for( d = i+2; d <= j-11; d++) {
      for( e = d+4; e <= j-7; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE) {

          gap_idej = gap_index(i,d,e,j,seqlength);
          for( c = i+1; c <= d-1; c++) {
            bp_penalty = IJ_bp_penalty;

            gap_cdej1 = gap_index( c,d,e,j-1,seqlength);
            if( Pg[ gap_idej] > 0) {
              prob =
                Pg[ gap_idej] *
                EXP_FUNC( -(ALPHA_1+ALPHA_2 + (c-i-1)*ALPHA_3 +
                        bp_penalty)/(kB*TEMP_K)) *
                Qgrs[ gap_cdej1] *
                ExplDangle( i+1, c-1, seq, seqlength)/
                Qg[ gap_idej];
              Pgrs[ gap_cdej1] += prob;
              if( prob > 1) {
                printf("r %d %d %Le\n", i, j, (long double) prob);
              }


            }
          }
        }
      }
    }

    //Case 7: Interior loop + multi both sides
    for( d = i+7; d <= j-11; d++) {
      for( e = d+4; e <= j-7; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE) {

          gap_idej = gap_index(i,d,e,j,seqlength);
          for( c = i+6; c <= d-1; c++) {
            bp_penalty = IJ_bp_penalty;

            pf_i1c1 = pf_index( i+1,c-1,seqlength);
            gap_cdej1 = gap_index( c,d,e,j-1,seqlength);
            if( Pg[ gap_idej] > 0) {
              prob =
                Pg[ gap_idej] *
                Qm[ pf_i1c1] *
                EXP_FUNC( -(ALPHA_1+ALPHA_2+bp_penalty)/(kB*TEMP_K)) *
                Qgrs[ gap_cdej1]/
                Qg[ gap_idej];
              Pm[ pf_i1c1] += prob;
              Pgrs[ gap_cdej1] += prob;
              if( prob > 1) {
                printf("s %d %d %Le\n", i, j, (long double) prob);
              }

            }
          }
        }
      }
    }

  } // Can Pair

}

/* ************************************** */
void MakePb_N5( int i, int j, int seq[], int seqlength,
                DBL_TYPE *Qm, DBL_TYPE *Qb, DBL_TYPE *Qp,
                DBL_TYPE *Pm, DBL_TYPE *Pb, DBL_TYPE *Pp){
  // This finds all possible internal loops (no pseudoknots)
  // closed on the "outside" by bases i and j, as well as all
  // multiloops


  DBL_TYPE prob;
  int d, e; // d - e is internal basepair
  DBL_TYPE bp_penalty = 0;
  int pf_ij = pf_index( i,j,seqlength);
  int pf_de, pf_i1d1;

  if( CanPair( seq[i], seq[j]) == TRUE && j-i+1 >= 5) {
    for( d = i+1; d <= j - 5; d++) {
      pf_i1d1 = pf_index( i+1,d-1,seqlength);
      for( e = d + 4; e <= j - 1; e++) {
        pf_de = pf_index( d,e,seqlength);

        if( CanPair( seq[d], seq[e]) == TRUE) {
          bp_penalty = 0.0;
          if( Pb[ pf_ij] > 0) {
            prob = Pb[ pf_ij]*ExplInternal( i, j, d, e, seq)*
              Qb[ pf_de ] / Qb[ pf_ij];
            Pb[ pf_de] += prob;
            if( prob > 1) {
              printf("t %d %d %Le\n", i, j,  (long double) prob);
            }

          }

          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty = AT_PENALTY;
          }
          if( seq[i] != BASE_C && seq[j] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          if( d>= i+6 && CanWCPair(seq[d], seq[e]) &&
              CanWCPair(seq[i], seq[j])) {
            if( Pb[ pf_ij] > 0) {
              prob = Pb[ pf_ij] * Qm[ pf_i1d1] *
                Qb[ pf_de] *
                EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + ALPHA_3*(j-e-1) + bp_penalty)/
                      (kB*TEMP_K) )*
                ExplDangle( e+1, j-1, seq, seqlength)/ Qb[ pf_ij];
              Pm[ pf_i1d1] += prob;
              Pb[ pf_de] += prob;
              if( prob > 1) {
                printf("u %d %d %Le\n", i, j, (long double) prob);
              }

            }
          }
        }
        if( CanWCPair(seq[i], seq[j]) && d<=j-6 &&
            e >= d+5) {

          bp_penalty = 0;
          if( seq[i] != BASE_C && seq[j] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          if( Pb[ pf_ij] > 0) {

            prob =
              Pb[ pf_ij] *
              EXP_FUNC( -( BETA_1M + ALPHA_1 + 3*ALPHA_2 +
                       (j-e-1 + d-i-1)*ALPHA_3 + bp_penalty)
                    /(kB*TEMP_K) )*
              ExplDangle( i+1, d-1, seq, seqlength) *
              ExplDangle( e+1, j-1, seq, seqlength) *
              Qp[ pf_de] / Qb[ pf_ij];
            Pp[ pf_de] += prob;

            if( prob > 1) {
              printf("v %d %d %Le\n", i, j, (long double) prob);
            }

            if( d>= i+6) {
              prob = Pb[ pf_ij] *
                Qm[ pf_i1d1] *
                Qp[ pf_de] *
                EXP_FUNC( -( BETA_1M + ALPHA_1 + 3*ALPHA_2 + bp_penalty +
                         (j - e - 1)*ALPHA_3) / (kB*TEMP_K) ) *
                ExplDangle( e+1, j-1, seq, seqlength)/ Qb[ pf_ij];

              Pm[ pf_i1d1] += prob;
              Pp[ pf_de] += prob;

              if( prob > 1) {
                printf("w %d %d %Le\n", i, j, (long double) prob);
              }
            }

          }

        }


      }

    }
  }

}
/* *********************************************** */

#ifdef O_N8

void MakePg_N8( int i, int j, char seq[], int seqlength, DBL_TYPE Qg[],
                DBL_TYPE Qm[], DBL_TYPE Pg[], DBL_TYPE Pm[]) {
  //  Make the gap matrix for Qg
  int c,d,e,f;
  int gap_idej, gap_cdef, pf_i1d1, pf_e1j1;
  DBL_TYPE IJ_bp_penalty = 0.0;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE prob;


  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }

  if( CanPair( seq[i], seq[j]) == TRUE) {

    //Case 1:  Simple Interior Loop (O_N6)
    for( d = i+2; d <= j-6; d++) {
      for( e = d+4; e <= j-2; e++) {
        gap_idej = gap_index(i,d,e,j,seqlength);
        if( CanPair( seq[d], seq[e] ) == TRUE ) {

          for( c = i+1; c <= d-1; c++) {
            for( f = e+1; f <= j-1; f++) {
              if( CanPair( seq[c], seq[f]) == TRUE) {

                gap_cdef = gap_index( c,d,e,f,seqlength);
                if( Pg[ gap_idej] > 0) {
                  prob = Pg[ gap_idej] *
                    EXP_FUNC( -InteriorEnergy( i, j, c, f, seq)/(kB*TEMP_K)) *
                    Qg[ gap_cdef] /
                    Qg[ gap_idej];

                  Pg[ gap_cdef] += prob;


                  if( prob > 1) {
                    printf("n %d %d %Le\n", i, j, (long double) prob);
                  }
                }

              }
            }
          }
        }
      }
    }
  }

  if( CanPair( seq[i], seq[j]) == TRUE &&
      CanWCPair(seq[i], seq[j])) {
    //Case 2:  Multiloop Left
    for( d = i+6; d <= j-5; d++) {
      pf_i1d1 = pf_index(i+1,d-1,seqlength);

      for( e = d+4; e <= j-1; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE
            && CanWCPair(seq[d], seq[e])) {
          gap_idej = gap_index(i,d,e,j,seqlength);

          bp_penalty = IJ_bp_penalty;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          if( Pg[ gap_idej] > 0) {
            prob =
              Pg[ gap_idej] *
              Qm[ pf_i1d1] *
              EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 +
                      (j-e-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
              ExplDangle( e+1, j-1, seq, seqlength) /
              Qg[ gap_idej];
            Pm[ pf_i1d1] += prob;
            if( prob > 1) {
              printf("n %d %d %Le\n", i, j,(long double) prob);
            }

          }
        }
      }
    }

    //Case 3:  Multiloop Right
    for( d = i+1; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE
            && CanWCPair(seq[d], seq[e])
            ) {
          gap_idej = gap_index(i,d,e,j,seqlength);

          bp_penalty = IJ_bp_penalty;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }
          pf_e1j1 = pf_index( e+1,j-1,seqlength);

          if( Pg[ gap_idej] > 0) {
            prob =
              Pg[ gap_idej] *
              Qm[ pf_e1j1] *
              EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 +
                      (d-i-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
              ExplDangle( i+1, d-1, seq, seqlength) /
              Qg[ gap_idej];
            Pm[ pf_e1j1] += prob;

            if( prob > 1) {
              printf("o %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }
      }
    }

    //Case 4: Multiloop Both Sides
    for( d = i+6; d <= j-10; d++) {
      pf_i1d1 = pf_index( i+1,d-1,seqlength);
      for( e = d+4; e <= j-6; e++) {
        if( CanPair( seq[d], seq[e] ) == TRUE
            && CanWCPair(seq[d], seq[e])
            ) {
          gap_idej = gap_index(i,d,e,j,seqlength);

          bp_penalty = IJ_bp_penalty;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty += AT_PENALTY;
          }

          pf_e1j1 = pf_index( e+1,j-1,seqlength);

          if( Pg[ gap_idej] > 0) {
            prob =
              Pg[ gap_idej] *
              Qm[ pf_i1d1] *
              Qm[ pf_e1j1] *
              EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + bp_penalty)/(kB*TEMP_K)) /
              Qg[ gap_idej];
            Pm[ pf_i1d1] += prob;
            Pm[ pf_e1j1] += prob;

            if( prob > 1) {
              printf("p %d %d %Le\n", i, j, (long double) prob);
            }

          }
        }
      }
    }

    //Case 5: Multiloop Left + More Qg
    for( d = i+7; d <= j-6; d++) {
      for( e = d+4; e <= j-2; e++) {
        gap_idej = gap_index(i,d,e,j,seqlength);
        if( CanPair( seq[d], seq[e] ) == TRUE) {

          for( c = i+6; c <= d-1; c++) {
            for( f = e+1; f <= j-1; f++) {
              if( CanPair( seq[c], seq[f]) == TRUE &&
                  CanWCPair(seq[c], seq[f])) {

                bp_penalty = IJ_bp_penalty;
                if( seq[c] != BASE_C && seq[f] != BASE_C) {
                  bp_penalty += AT_PENALTY;
                }

                if( Pg[ gap_idej] > 0) {
                  prob = Pg[ gap_idej]*
                    Qm[ pf_index( i+1, c-1, seqlength)]*
                    EXP_FUNC( -(ALPHA_1+2*ALPHA_2+(j-f-1)*ALPHA_3 + bp_penalty)
                          /(kB*TEMP_K)) *
                    Qg[ gap_index( c, d, e, f, seqlength)] *
                    ExplDangle( f+1, j-1, seq, seqlength)/ Qg[gap_idej];
                  Pm[ pf_index( i+1, c-1, seqlength)] += prob;
                  Pg[ gap_index( c, d, e, f, seqlength)] += prob;
                }

              }
            }
          }
        }

      }
    }
    //Case 6: Multiloop Right + More Qg
    for( d = i+2; d <= j-11; d++) {
      for( e = d+4; e <= j-7; e++) {
        gap_idej = gap_index(i,d,e,j,seqlength);
        if( CanPair( seq[d], seq[e] ) == TRUE) {

          for( c = i+1; c <= d-1; c++) {
            for( f = e+1; f <= j-6; f++) {
              if( CanPair( seq[c], seq[f]) == TRUE &&
                  CanWCPair(seq[c], seq[f])) {

                bp_penalty = IJ_bp_penalty;
                if( seq[c] != BASE_C && seq[f] != BASE_C) {
                  bp_penalty += AT_PENALTY;
                }

                if( Pg[ gap_idej] > 0) {
                  prob = Pg[ gap_idej] *
                    Qm[ pf_index( f+1, j-1, seqlength)]*
                    EXP_FUNC( -(ALPHA_1+2*ALPHA_2+(c-i-1)*ALPHA_3+bp_penalty)
                          /(kB*TEMP_K)) *
                    Qg[ gap_index( c, d, e, f, seqlength)] *
                    ExplDangle( i+1, c-1, seq, seqlength) /Qg[ gap_idej];
                  Pm[  pf_index( f+1, j-1, seqlength)] += prob;
                  Pg[ gap_index( c, d, e, f, seqlength)] += prob;
                }

              }
            }
          }
        }
      }
    }


    //Case 7: Multiloop Both sides + More Qg
    for( d = i+7; d <= j-11; d++) {
      for( e = d+4; e <= j-7; e++) {
        gap_idej = gap_index(i,d,e,j,seqlength);
        if( CanPair( seq[d], seq[e] ) == TRUE) {

          for( c = i+6; c <= d-1; c++) {
            for( f = e+1; f <= j-6; f++) {
              if( CanPair( seq[c], seq[f]) == TRUE &&
                  CanWCPair(seq[c], seq[f])) {

                bp_penalty = IJ_bp_penalty;
                if( seq[c] != BASE_C && seq[f] != BASE_C) {
                  bp_penalty += AT_PENALTY;
                }

                if( Pg[ gap_idej] > 0) {
                  prob = Pg[ gap_idej] *
                    Qm[ pf_index( i+1, c-1, seqlength)]*
                    Qm[ pf_index( f+1, j-1, seqlength)]*
                    EXP_FUNC( -(ALPHA_1+2*ALPHA_2+bp_penalty)
                          /(kB*TEMP_K)) *
                    Qg[ gap_index( c, d, e, f, seqlength)]/Qg[ gap_idej];
                  Pm[ pf_index( i+1, c-1, seqlength)] += prob;
                  Pm[ pf_index( f+1, j-1, seqlength)] += prob;
                  Pg[ gap_index( c, d, e, f, seqlength)] += prob;
                }

              }
            }
          }
        }
      }
    }

  } // Can Pair
}
#endif //O_N8


/* ***************************************** */

void recalculateQgIx( int i, int j, int d, int e, int size, int qgix, int seq[],
                      int seqlength, DBL_TYPE *QgIx, DBL_TYPE *Qg,
                      int side) {


  int c,f;
  DBL_TYPE tmp = 0;

  int L1min = 5;
  int L2min = 4;

  if( side == 2) L2min = 5;

  for( c = i+L1min+1; c <= d-1; c++) {
    f = c-i+j-2-size;
    if( j-f-1>=L2min && f >= e+1 &&
        CanPair( seq[c], seq[f]) == TRUE) {
      tmp += Qg[ gap_index(c,d,e,f,seqlength)]*
        EXP_FUNC( -1*(InteriorEnergyFull( i,j,c,f,seq, FALSE))/
              (kB*TEMP_K));

    }
  }
  QgIx[ qgix] = tmp;

}

/* *********** */

void prManageQgIx( DBL_TYPE **QgIx, DBL_TYPE **QgIx_1,
                   DBL_TYPE **QgIx_2, DBL_TYPE **PgIx, DBL_TYPE **PgIx_1,
                   DBL_TYPE **PgIx_2,
                   float **preX, float **preX_1, float **preX_2,
                   int d, int seqlength) {
  // Allocate and deallocate QgIx matrices

  int maxStorage;
  int dTest;
  DBL_TYPE *temp;

  float *tempf;

  maxStorage = 0;
  for( dTest = d; dTest >= d-2; dTest--) {
    maxStorage = MAX( maxStorage, (seqlength-dTest)*(dTest-5)*
                      ((dTest-1)*(dTest-2)/2) );
  }

  if( d == seqlength - 1 && d>=16) {

    free( *QgIx);
    free( *QgIx_1);
    free( *QgIx_2);

    *QgIx =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *QgIx_1 =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *QgIx_2 =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    *PgIx = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *PgIx_1 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *PgIx_2 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    *preX = (float *) calloc( maxStorage, sizeof( float));
    *preX_1 = (float *) calloc( maxStorage, sizeof( float));
    *preX_2 = (float *) calloc( maxStorage, sizeof( float));

    if( *QgIx_2 == NULL || *PgIx == NULL || *PgIx_1 == NULL ||
        *PgIx_2 == NULL || *preX == NULL || *preX_1 == NULL ||
        *preX_2 == NULL) {
      printf("!! Error in QgIx_2, PgIx, preX allocation\n");

      exit(1);
    }
  }
  else if( d >= 16) { // every case beyond the first

    temp = *QgIx;
    *QgIx = *QgIx_1;
    *QgIx_1 = *QgIx_2;
    free(temp);
    temp = NULL;

    *QgIx_2 = (DBL_TYPE *) calloc( maxStorage, sizeof(DBL_TYPE) );

    temp = *PgIx;
    *PgIx = *PgIx_1;
    *PgIx_1 = *PgIx_2;
    free( temp);
    temp = NULL;

    *PgIx_2 =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    tempf = *preX;
    *preX = *preX_1;
    *preX_1 = *preX_2;

    free( tempf);
    tempf = NULL;

    *preX_2 =
      (float *) calloc( maxStorage, sizeof( float));

    if( *QgIx_2 == NULL || *PgIx_2 == NULL || *preX_2 == NULL) {
      printf("Error in QgIx_2, PgIx_2, preX_2 allocation\n");
      exit(1);
    }

  }
}

/* ************************************************ */


void prFastILoopsN5( int i, int j, int seq[], int seqlength,
                     DBL_TYPE *Qg, DBL_TYPE *QgIx, DBL_TYPE *QgIx_2,
                     DBL_TYPE *Pg,
                     DBL_TYPE *PgIx, DBL_TYPE *PgIx_2,
                     float *preX, float *preX_2) {

  int d, e, size;
  int L = j - i + 1;
  int qgix, qgix2;
  DBL_TYPE ExplInteriorMM = NAD_INFINITY;
  int gap_idej;
  DBL_TYPE energy;
  int c, f; //Internal pair.(c, f will be restricted to special cases)
  int L1, L2; //size parameters: L1 + L2 = size, L1 = c-i-1, L2 = j-f-1
  int gap_cdef;
  //extern DBL_TYPE loop37[];
  //extern DBL_TYPE asymmetry_penalty[];
  //extern DBL_TYPE max_asymmetry;
  //int asymmetry, asymmetry_index;
  DBL_TYPE pr;
  DBL_TYPE oldSizeEnergy, newSizeEnergy;

  float precisionLost;

  if( CanPair( seq[i], seq[j]) == TRUE) {
    ExplInteriorMM = EXP_FUNC( -InteriorMM( seq[i], seq[j], seq[i+1],
                                        seq[j-1])/(kB*TEMP_K));

    for( d = i+1; d <= j-5; d++) {
      for( e = d+4; e <= j-1; e++) {
        if(  CanPair( seq[d], seq[e]) == TRUE) {

          gap_idej = gap_index( i,d,e,j, seqlength);

          MakePg_Inextensible( i,j,d,e, seq, seqlength, Qg, Pg);
        }
      }
    }
  }


  //Add in special "root" cases( i ==0 || j == selqength -1) to QgIx
  if( (i == 0 || j == seqlength - 1) && L>=17 /* && L <= seqlength - 2 */) {
    for( d = i+6; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {
        if(  CanPair( seq[d], seq[e]) == TRUE) {

          for( c = i+5; c <= d-1; c++) {
            for( f = e+1; f <= j - 5; f++) {
              if( CanPair( seq[c], seq[f]) == TRUE) {

                L1 = c-i-1;
                L2 = j-f-1;
                size = L1 + L2;


                energy = asymmetryEfn( L1, L2, size);

                energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
                /*Exclude the i-j stacking energy here, just in case i-j
                  don't pair */

                gap_cdef = gap_index( c,d,e,f, seqlength);
                qgix = QgIxIndex( j-i, i, size, d, e, seqlength);

                QgIx[ qgix ] +=
                  EXP_FUNC(-energy/(kB*TEMP_K))*Qg[ gap_cdef];


              }
            }
          }
        }
      }
    }
  }


  //Calculate PgIx using Pb
  if( CanPair(seq[i], seq[j]) == TRUE) {
    for( d = i+1; d <= j-5; d++) {
      for( e = d+4; e <= j-1; e++) {
        if(  CanPair( seq[d], seq[e]) == TRUE) {

          gap_idej = gap_index( i,d,e,j, seqlength);

          for( size = 8; size <= L - 9; size++) {
            qgix = QgIxIndex( j-i,i,size,d,e,seqlength);

            pr = Pg[ gap_idej]* QgIx[ qgix] *
              ExplInteriorMM / Qg[ gap_idej];
            PgIx[ qgix] += pr;

          }
        }
      }
    }
  }

  if( L >= 17) {
    for( d = i+6; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {
        if(  CanPair( seq[d], seq[e]) == TRUE) {

          //Next calculate Pb using PgIx and QgIx for L1 == 4 || L2 == 4
          //Case 1:  L1 = 4, L2 >= 4;
          L1 = 4;
          c = i + L1 + 1;
          for( L2 = 4; L2 <= j - e - 2; L2++) {
            size = L1 + L2;
            f = j - L2 - 1;
            qgix =  QgIxIndex( j-i, i, size, d,e,seqlength);

            if( CanPair( seq[c], seq[f]) == TRUE) {

              energy = asymmetryEfn( L1, L2, size);

              energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
              /*Exclude the i-j stacking energy here, just in case i-j
                don't pair */

              gap_cdef = gap_index(c,d,e,f,seqlength);
              if( QgIx[ qgix] > 0) {
                pr = PgIx[ qgix ] *
                  EXP_FUNC(-energy/(kB*TEMP_K))*Qg[ gap_cdef] / QgIx[ qgix];
                Pg[ gap_cdef] += pr;

                PgIx[ qgix] -= pr;
              }

              //Recalculate QgIx
              precisionLost = subtractLongDouble( &(QgIx[ qgix]),
                                                  EXP_FUNC(-energy/(kB*TEMP_K))*
                                                  Qg[ gap_cdef]);

              //printf("%d %d %d %d %d %.2f %.10Le ", i, j, d, e, size, precisionLost, (long double) QgIx[ qgix]);

              preX[ qgix] += precisionLost;
              if( preX[ qgix] >= MAXPRECERR) {
                recalculateQgIx( i, j, d, e, size, qgix, seq, seqlength, QgIx, Qg, 1);
                preX[ qgix] = 0.0;
                //printf("%.10Le\n", (long double) QgIx[ qgix]);
              }
              //else { printf("\n");}

            }
          }

          //Case 2  L1 > 4, L2 = 4
          L2 = 4;
          f = j - L2 -1;
          for( L1 = 5; L1 <= d-i-2; L1++) {
            size = L1 + L2;
            c = i + L1 + 1;

            qgix =  QgIxIndex( j-i, i, size, d,e,seqlength);

            if( CanPair( seq[c], seq[f]) == TRUE) {

              energy = asymmetryEfn( L1, L2, size);

              energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
              /*Exclude the i-j stacking energy here, just in case i-j
                don't pair */

              gap_cdef = gap_index(c,d,e,f,seqlength);
              if( QgIx[ qgix] > 0) {
                pr = PgIx[ qgix ] *
                  EXP_FUNC(-energy/(kB*TEMP_K))*Qg[ gap_cdef] / QgIx[ qgix];
                Pg[ gap_cdef] += pr;

                PgIx[ qgix] -= pr;
              }

              precisionLost = subtractLongDouble( &(QgIx[ qgix]),
                                                  EXP_FUNC(-energy/(kB*TEMP_K))*
                                                  Qg[ gap_cdef]);
              //printf("-- %d %d %d %d %d %.2f %.10Le ", i, j, d, e, size, precisionLost, (long double) QgIx[ qgix]);

              preX[ qgix] += precisionLost;
              if( preX[ qgix] >= MAXPRECERR) {
                recalculateQgIx( i, j, d, e, size, qgix, seq, seqlength, QgIx, Qg, 2);
                preX[ qgix] = 0.0;

                //printf("%.10Le\n", (long double) QgIx[ qgix]);
              }
              //else { printf("\n");}

            }
          }
        }

        for( size = 10; size <= L - 9; size++) {
          if( size <= 30) {
            oldSizeEnergy = loop37[ size - 1];
          }
          else {
            oldSizeEnergy = loop37[ 30 - 1];
            oldSizeEnergy += sizeLog (size); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);
          }
          if( size - 2 <= 30) {
            newSizeEnergy = loop37[ size-2 - 1];
          }
          else {
            newSizeEnergy = loop37[ 30 - 1];
            newSizeEnergy += sizeLog (size-2); //1.75*kB*TEMP_K*LOG_FUNC( (size-2)/30.0);
          }

          qgix  = QgIxIndex( j-i, i, size, d, e, seqlength);
          qgix2 = QgIxIndex( j-i-2, i+1, size-2, d, e, seqlength);

          QgIx_2[ qgix2] =
            QgIx[ qgix] *
            EXP_FUNC( -(newSizeEnergy - oldSizeEnergy)/(kB*TEMP_K));

          PgIx_2[ qgix2] = PgIx[ qgix];
          preX_2[ qgix2] = preX[ qgix];

        }
      }
    }
  }

}


/* *********************************** */

void MakePg_Inextensible( int i, int j,int d, int e, int seq[], int seqlength,
                          DBL_TYPE *Qg, DBL_TYPE *Pg) {

  DBL_TYPE energy;
  int c, f; //Internal pair.(c, f will be restricted to special cases)
  int L1, L2; //size parameters: L1 + L2 = size, L1 = c-i-1, L2 = j-f-1
  DBL_TYPE pr;
  int gap_idej = gap_index(i,d,e,j,seqlength);
  int gap_cdef;

  if( CanPair( seq[i], seq[j]) == FALSE) {
    return;
  }

  /* First consider small loops */
  // L1 <= 3, L2 <= 3
  for( L1 = 0; L1 <= MIN(3, d-i-2); L1++) {
    c = i + L1 + 1;
    for( L2 = 0; L2 <= MIN( 3, j-e-2); L2++) {
      f = j - L2 - 1;
      if( CanPair( seq[c], seq[f]) == TRUE) { //c < d && f > e

        energy = InteriorEnergy( i, j, c, f, seq);
        gap_cdef = gap_index(c,d,e,f,seqlength);

        pr = Pg[ gap_idej]*EXP_FUNC( -energy/(kB*TEMP_K)) *
          Qg[ gap_cdef]/ Qg[ gap_idej];
        Pg[ gap_cdef] += pr;

      }
    }
  }

  /* Next consider large bulges or large asymmetric loops */
  // Case 2a  L1 = 0,1,2,3, L2 >= 4;
  for( L1 = 0; L1 <= MIN(3, d-i-2); L1++) {
    c = i + L1 + 1;
    for( L2 = 4; L2 <= j - e - 2; L2++) {
      f = j - L2 - 1;
      if( CanPair( seq[c], seq[f]) == TRUE) { //c < d && f > e

        gap_cdef = gap_index(c,d,e,f,seqlength);

        energy = InteriorEnergy( i, j, c, f, seq);
        pr = Pg[ gap_idej] * EXP_FUNC( -energy/(kB*TEMP_K)) *
          Qg[ gap_cdef]/Qg[ gap_idej];

        Pg[ gap_cdef] += pr;

      }
    }
  }

  // Case 2b L1 >= 4, L2 = 0,1,2,3;
  for( L2 = 0; L2 <= MIN(3,j-e-2) ; L2++) {
    f = j - L2 - 1;
    for( L1 = 4; L1 <= d-i-2; L1++) {
      c = i + L1 + 1;
      if( CanPair( seq[c], seq[f]) == TRUE) { //c < d && f > e

        gap_cdef = gap_index(c,d,e,f,seqlength);

        energy = InteriorEnergy( i, j, c, f, seq);
        pr = Pg[ gap_idej]*EXP_FUNC( -energy/(kB*TEMP_K)) *
          Qg[ gap_cdef]/Qg[ gap_idej];
        Pg[ gap_cdef] += pr;

      }
    }
  }

}

/* ****************************************** */

void prInteriorLoopsN4MS(int i, int j, int seq[], int seqlength,
                         DBL_TYPE *Qb, DBL_TYPE *Pb, int *nicks) {


  int d, e, pf_de;
  int pf_ij = pf_index(i,j,seqlength);
  double pr;
  double energy;

  int leftNick = -1;
  int rightNick = -1;
  int leftIndex, rightIndex;

  findNicks( nicks, &leftIndex, &rightIndex, NULL, i+1, j-1);

  if( leftIndex != -1) leftNick = nicks[ leftIndex];
  if( rightIndex != -1) rightNick = nicks[ rightIndex];

  if( CanPair( seq[ i], seq[j]) == TRUE) {
    for( d = i+1; d <= j-2; d++) {
      if( leftNick != -1 && leftNick <= d-1) break;
      for( e = d+1; e <= j-1; e++) {
        if( rightNick != -1 && rightNick >= e) continue;
        if( CanPair( seq[d], seq[e]) == TRUE) {

          pf_de = pf_index( d, e, seqlength);

          energy = InteriorEnergy( i, j, d, e, seq);

          if( Qb[ pf_ij] > 0) {
            pr = Pb[ pf_ij] * EXP_FUNC( -energy/(kB*TEMP_K)) *
              Qb[ pf_de] / Qb[ pf_ij];
            Pb[ pf_de] += pr;
          }


        }
      }
    }
  }
}

void findNicks( int *nicks, int *leftNickIndex, int *rightNickIndex,
                int *nNicks, int leftEdge, int rightEdge) {

  int k;

  *leftNickIndex = *rightNickIndex = -1;

  for( k = 0; k < MAXSTRANDS; k++) {
    if( nicks[k] == -1) return;
    if( nicks[k] >= leftEdge-1 && nicks[k] <= rightEdge) {
      if( nicks[k] != leftEdge-1 && nNicks != NULL) {
        (*nNicks)++;
      }
      if( *leftNickIndex == -1) {
        *leftNickIndex = k;
      }
      *rightNickIndex = k;
    }
  }

  return;
}

/* ***** */
