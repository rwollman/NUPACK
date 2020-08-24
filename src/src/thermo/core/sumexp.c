/*
  sumexp.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2006, Justin Bois 1/2007

  This file contains the functions that calculate the sum of exponentials
  in order to calculate a given partition function matrix. 

  See pfuncUtilsHeader.h for more specific descriptions of each function.
*/

#include "sumexp.h"

extern unsigned int seqHash;
/* ******************************************* */
DBL_TYPE ExplHairpin( int i, int j, int seq[], int seqlength, int **etaN) {
  //this version disallows nicks here

  DBL_TYPE energy = 0;
  int index;
  int nNicks = 0;

  index = EtaNIndex(i+0.5, j-0.5, seqlength);
  nNicks = etaN[ index][0];

  if( nNicks >= 1) return 0;

  if( nNicks == 0 && j-i <= 3) {
    return 0;
  }

  if( nNicks == 0) {
    energy = HairpinEnergy( i, j, seq);
  }

  if( energy == NAD_INFINITY) {
    return 0.0;
  }

  return EXP_FUNC( -energy/( kB*TEMP_K) );
}

/* ********************* */
DBL_TYPE SumExpMultiloops( int i, int j, int seq[], 
                          DBL_TYPE *Qms, DBL_TYPE *Qm, int seqlength,
                          int **etaN){
  // Decomposes the region inside pair i,j into multiloops, i.e.
  // and excludes the possibility of "top level" nicks

  DBL_TYPE sum_exp = 0.0;
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE extraTerms;

  int d; // d is the left base of a rightmost paired base between i, j.
 
  if( CanWCPair(seq[i], seq[j])) {
    for( d = i+3; d <= j - 2; d++) {
      //reset loop parameters
      bp_penalty = 0.0;
      if( etaN[EtaNIndex_same( d-0.5, seqlength)][0] == 0 ) {
        
        if( seq[i] != BASE_C  && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }

        extraTerms = EXP_FUNC( -( ALPHA_1 + ALPHA_2 + bp_penalty) 
                          / (kB*TEMP_K) );
        if( DNARNACOUNT == COUNT) 
          extraTerms = 1;

        sum_exp += Qm[ pf_index( i+1, d-1, seqlength)] *
          Qms[ pf_index(d, j-1, seqlength)] * extraTerms;
      }
      
    }
  }

  return sum_exp;  
}

/* *********************************************** */

DBL_TYPE SumExpExteriorLoop( int i,int j, int seq[], int seqlength, 
                            DBL_TYPE *Q, int *nicks, int **etaN) {

  DBL_TYPE sumExp = 0.0;
  DBL_TYPE bp_penalty = 0.0;
  int multiNick = -1;
  int index_ij;
  int leftIndex;
  int nNicks;
  int n;
  int iNicked, jNicked;
  DBL_TYPE extraTerms;

  index_ij = EtaNIndex(i+0.5, j-0.5, seqlength);
  iNicked = jNicked = FALSE;

  if( etaN[ EtaNIndex_same(j-0.5, seqlength)][0] != 0) {
    jNicked = TRUE;
  }

  if( etaN[ EtaNIndex_same(i+0.5, seqlength)][0] != 0) {
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
        EXP_FUNC( -1*(bp_penalty)/(kB*TEMP_K));

      if( DNARNACOUNT == COUNT) 
        if( extraTerms != 0) extraTerms = 1;

      if( (iNicked == FALSE && jNicked == FALSE) ||
         (i == j - 1) || 
         (multiNick == i && jNicked == FALSE) ||
         (multiNick == j-1 && iNicked == FALSE ) ) {

           sumExp +=
             Q[ pf_index(i+1, multiNick, seqlength)]*
             Q[ pf_index( multiNick+1, j-1, seqlength)] *
             extraTerms;
      }
    }
  }
  return sumExp;  
}

/* *********************************************** */

void fastILoops( int i, int j, int L, int seqlength, int seq[],
                 int **etaN, DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2,
                 DBL_TYPE * Qb_bonus) {

  int size;
  int pf_ij = pf_index( i, j, seqlength);
  DBL_TYPE extraTerms;

  int isEndNicked = FALSE;
  if( etaN[ EtaNIndex( i-0.5,i-0.5, seqlength)][0] == 1 || 
     etaN[ EtaNIndex( j+0.5,j+0.5, seqlength)][0] == 1) 
    isEndNicked = TRUE;
  if( L >= 12) {
    makeNewQx( i, j, seq, seqlength, etaN, Qb, Qx);
  }

  //Use extensible cases              
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    for( size = 8; size <= L - 4; size++) {
      
      extraTerms = EXP_FUNC( -InteriorMM( seq[i], seq[j], seq[i+1], 
                                     seq[j-1])/(kB*TEMP_K));
      if( DNARNACOUNT == COUNT) 
        extraTerms = 1;
      
      Qb[ pf_ij] += Qb_bonus[pf_ij] * 
        Qx[ fbixIndex( j-i, i, size, seqlength)] * 
        extraTerms;
    }
  }

  if( L >= 12 && i != 0 && j != seqlength -1 && isEndNicked == FALSE) {
    extendOldQx( i, j, seqlength, Qx,Qx_2);
  }
  
  /* Add in inextensible cases */  
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    //first check inextensible cases

    Qb[ pf_ij] += Qb_bonus[pf_ij] * 
      SumExpInextensibleIL( i,j, seq, seqlength, Qb,  etaN);

  } 
}

/* *************** */

/* Qs, Qms  Recursion */

void MakeQs_Qms( int i, int j, int seq[], int seqlength, 
                DBL_TYPE *Qs, DBL_TYPE *Qms, DBL_TYPE *Qb,
                int *nicks, int **etaN) {

  int d; //base pair is i,d
  DBL_TYPE bp_penalty = 0.0;
  int pf_ij = pf_index( i, j, seqlength);

  DBL_TYPE extraTerms;
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
    bp_penalty = 0.0;
    
    if( CanPair( seq[i], seq[ d]) == TRUE &&
       CanWCPair(seq[i], seq[d])) {
         
         if( seq[i] != BASE_C && seq[d] != BASE_C) {
           bp_penalty = AT_PENALTY;
         }

         extraTerms = EXP_FUNC( -(NickDangle( d+1,j,nicks, etaN,
                              FALSE, seq,seqlength) + 
                              bp_penalty)/(kB*TEMP_K) );

         if( DNARNACOUNT == COUNT) 
           extraTerms = 1;

         Qs[ pf_ij] += Qb[ pf_index( i, d, seqlength) ] * 
           extraTerms;

         // ******************** 

         extraTerms =  ExplDangle( d+1, j, seq, seqlength) * 
           EXP_FUNC( -(bp_penalty + ALPHA_2 + ALPHA_3*(j-d))/(kB*TEMP_K) );

         if( DNARNACOUNT == COUNT) 
           extraTerms = 1;
         Qms[ pf_ij] += Qb[ pf_index( i, d, seqlength) ] * 
           extraTerms;
    }
  }
}


/* ******************************* */
/* Q, Qm Recursions */
void MakeQ_Qm_N3( int i, int j, int seq[], int seqlength, 
                 DBL_TYPE *Q, DBL_TYPE *Qs, 
                 DBL_TYPE *Qms, DBL_TYPE *Qm,
                 int *nicks, int **etaN) {
  // static DBL_TYPE *ExplDanglePre;
  // static int ExplInited=0;
  int d; // ,e;//left base of rightmost base pair.
  int pf_ij = pf_index( i, j, seqlength);

  DBL_TYPE extraTerms;
  Q[ pf_ij] = NickedEmptyQ( i, j, nicks, seq, seqlength, etaN);

  for( d = i; d <= j - 1; d++) {
    if( etaN[ EtaNIndex_same(d-0.5, seqlength)][0] == 0 || d == i ) {

      Q[ pf_ij] += Q[ pf_index(i, d-1, seqlength)] *
        Qs[ pf_index( d, j, seqlength)];

      if( DNARNACOUNT == COUNT) 
        extraTerms = 1;
      else 
        extraTerms = ExplDangle( i, d-1, seq, seqlength) *
          EXP_FUNC( -(ALPHA_3)*(d-i)/(kB*TEMP_K));

      if( etaN[ EtaNIndex_same( d-0.5, seqlength)][0] == 0) { 
        //otherwise Qm not possible
       
        if( etaN[ EtaNIndex(i+0.5, d-0.5, seqlength)][0] == 0 ) {
          Qm[ pf_ij] += Qms[ pf_index( d, j, seqlength)] *
           extraTerms; //Single Pair
        }

        if( d >= i+2) {
         Qm[ pf_ij]+= Qm[ pf_index( i, d - 1, seqlength) ] *
           Qms[ pf_index( d, j, seqlength) ];
        }
      }
    }
  }
}

/* ******************************************* */

/* Functions in Q recursion */
// must be calculated after Qb, Qpk of same length

void makeNewQx( int i, int j, int seq[], int seqlength, 
               int **etaN, DBL_TYPE Qb[], DBL_TYPE Qx[]) {
                 
  /*Determine the new entries of Qx(i,j,size) that are not extended 
  versions of Qx(i+1, j-1, size-2) */

  DBL_TYPE energy;
  int d, e; //Internal pair.(d, e will be restricted to special cases)

  int size, L1, L2; //size parameters: L1 + L2 = size, L1 = h-i-1, L2 = j-m-1

  //Add in all the cases that are not an extended version of an
  //extensible case.

  //Case 1:  L1 = 4, L2 >= 4;
  L1 = 4;
  d = i + L1 + 1;

  for( L2 = 4; L2 <= j - d - 2; L2++) {
    size = L1 + L2;
    e = j - L2 - 1;

    if( CanPair( seq[d], seq[e]) == TRUE &&
      (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
      (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

        energy = asymmetryEfn( L1, L2, size) + InteriorMM( seq[e], seq[d], seq[e+1], seq[d-1]);
        /*Exclude the i-j stacking energy here, just in case i-j 
        don't pair */

        if( DNARNACOUNT == COUNT) 
          energy = 0;

        Qx[ fbixIndex( j-i, i, size, seqlength) ] += 
          EXP_FUNC(-energy/(kB*TEMP_K))*Qb[ pf_index(d, e, seqlength)];
    }
  }

  //Case 2  L1 > 4, L2 = 4
  L2 = 4;
  e = j - L2 -1;
  for( L1 = 5; L1 <= e-i-2; L1++) {
    size = L1 + L2;
    d = i + L1 + 1;

    if( CanPair( seq[d], seq[e]) == TRUE &&
      (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
      (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

        energy = asymmetryEfn( L1, L2, size) + InteriorMM( seq[e], seq[d], seq[e+1], seq[d-1]);
        /*Exclude the i-j stacking energy here, just in case i-j 
        don't pair */

        if( DNARNACOUNT == COUNT) 
          energy = 0.0;

        Qx[ fbixIndex( j-i, i, size, seqlength)] +=
          EXP_FUNC(-energy/(kB*TEMP_K))*Qb[ pf_index(d, e, seqlength)];
    }
  }
}

/* ************************** */
void extendOldQx( int i, int j, int seqlength, DBL_TYPE Qx[], DBL_TYPE Qx_2[]) {
  /* Extends all entries of Qx */
  
  int size;
  DBL_TYPE oldSizeEnergy;
  DBL_TYPE newSizeEnergy;
  
  for( size = 8; size <= (j - i + 1) - 4; size++) {
    if( size <= 30) {
      oldSizeEnergy = loop37[ size - 1];
    }
    else {
      oldSizeEnergy = loop37[ 30 - 1];
      oldSizeEnergy += sizeLog(size);  //1.75*kB*TEMP_K*log( size/30.0); 
    }

    if( size + 2 <= 30) {
      newSizeEnergy = loop37[ size+2 - 1];
    }
    else {
      newSizeEnergy = loop37[ 30 - 1];
      newSizeEnergy += sizeLog (size+2); //1.75*kB*TEMP_K*log( (size+2)/30.0); 
    }

    if( DNARNACOUNT == COUNT) 
      newSizeEnergy = oldSizeEnergy;
    
    Qx_2[ fbixIndex( j-i+2, i-1, size+2, seqlength)] = 
      Qx[ fbixIndex( j-i, i, size, seqlength)] * 
      EXP_FUNC( -(newSizeEnergy - oldSizeEnergy)/(kB*TEMP_K));
  }
}


/* ************************ */
DBL_TYPE SumExpInextensibleIL( int i, int j, int seq[], int seqlength, 
                              DBL_TYPE Qb[], int **etaN) {
  /* This finds the minimum energy IL that has a special energy 
  calculation, i.e. small loops, bulge loops or GAIL case.  None of 
  these are allowed to be nicked
  */

  DBL_TYPE energy;

  // int nse=0;
  int d, e; //Internal pair.(h, m will be restricted to special cases)  
  int L1, L2; //size parameters: L1 + L2 = size, L1 = h-i-1, L2 = j-m-1

  DBL_TYPE sumexp = 0.0;

  /* Consider "small" loops with special energy functions */

  for( L1 = 0; L1 <= 3; L1++) {
    d = i + L1 + 1;
    for( L2 = 0; L2 <= MIN( 3, j-d-2); L2++) {
      e = j - L2 - 1;

      if( CanPair( seq[d], seq[e]) == TRUE &&
         (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
         (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

           energy = InteriorEnergy( i, j, d, e, seq);

           sumexp += EXP_FUNC( -energy/(kB*TEMP_K)) *
             Qb[ pf_index( d, e, seqlength)];
      }
    }
  }

  /* Next consider large bulges or large asymmetric loops */
  // Case 2a  L1 = 0,1,2,3, L2 >= 4;
  for( L1 = 0; L1 <= 3; L1++) {
    d = i + L1 + 1;
    for( L2 = 4; L2 <= j - d - 2; L2++) {
      e = j - L2 - 1;

      if( CanPair( seq[d], seq[e]) == TRUE &&
         (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
         (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

           energy = InteriorEnergy( i, j, d, e, seq);

           sumexp += EXP_FUNC( -energy/(kB*TEMP_K)) *
             Qb[ pf_index( d, e, seqlength)]; 
      }
    }
  }

  // Case 2b L1 >= 4, L2 = 0,1,2,3;
  for( L2 = 0; L2 <= 3; L2++) {
    e = j - L2 - 1;
    for( L1 = 4;  L1 <= e - i - 2; L1++) {
      d = i + L1 + 1;

      if( CanPair( seq[d], seq[e]) == TRUE &&
         (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
         (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {

           energy = InteriorEnergy( i, j, d, e, seq);

           sumexp += EXP_FUNC( -energy/(kB*TEMP_K)) *
             Qb[ pf_index( d, e, seqlength)];
      }
    }
  }

  return sumexp;
}


/* ******************** */
DBL_TYPE SumExpInterior_Multi( int i, int j, int seq[], int seqlength, 
                              DBL_TYPE Qm[], DBL_TYPE Qb[] ){
  // This finds all possible internal loops (no pseudoknots)
  // closed on the "outside" by bases i and j, as well as all 
  // multiloops.  Ignores nicks

  DBL_TYPE sum_expl = 0.0;
  int d, e; // d - e is internal basepair 
  DBL_TYPE bp_penalty = 0;
  // int S1 = j-i+1;
  // int S2;
  // int S3;

  for( d = i+1; d <= j - 5; d++) {
    for( e = d + 4; e <= j - 1; e++) {
      //  S2 = e-d+1;
      //  S3 = d-i-1;
      
      if( CanPair( seq[d], seq[e]) == TRUE) {
        bp_penalty = 0.0;

        sum_expl += 
          ExplInternal( i, j, d, e, seq) *
          Qb[ pf_index( d, e, seqlength) ];

        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }
        if( seq[i] != BASE_C && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }

        if( d>= i+6 && CanWCPair(seq[d], seq[e]) &&
           CanWCPair(seq[i], seq[j])) {

             sum_expl += 
               Qm[ pf_index(i+1, d-1, seqlength)] *
               Qb[ pf_index( d, e, seqlength)] *
               EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + ALPHA_3*(j-e-1) + bp_penalty)/
                    (kB*TEMP_K) )*
               ExplDangle( e+1, j-1, seq, seqlength);
        }
      }
    }
  }
  return sum_expl;
}

/* *********************************************** */
#ifdef N4
void MakeQ_Qm_N4( int i, int j, int seq[], int seqlength, 
                 DBL_TYPE *Q, DBL_TYPE *Qm, DBL_TYPE *Qb ){

  int d, e; // d - e is internal basepair 
  DBL_TYPE bp_penalty = 0;
  int pf_ij = pf_index(i, j, seqlength);
  //int S1 = j-i+1;
  // int S2;
  // int S3;

  Q[ pf_ij] = /*scale(S1)* */ ExplDangle(i, j, seq, seqlength);  //Empty Graph

  for( d = i; d <= j - 4; d++) {
   // S3 = d-i;
   for( e = d + 4; e <= j; e++) {
     //S2 = e-d+1;

     if( CanPair( seq[d], seq[e]) == TRUE &&
          CanWCPair(seq[d], seq[e]) ) {
        bp_penalty = 0;
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }

        Q[ pf_ij] += //scale(S1-S2-S3)*
          Q[ pf_index(i, d-1, seqlength)] *
          Qb[ pf_index( d, e, seqlength) ] *
          EXP_FUNC( -bp_penalty/(kB*TEMP_K) ) *
          ExplDangle( e+1, j, seq, seqlength); 

        Qm[ pf_ij] += //scale( S1-S2) *
          EXP_FUNC( -(ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty)/
               (kB*TEMP_K) )*
          Qb[ pf_index( d, e, seqlength)] *
          ExplDangle( e+1, j, seq, seqlength) *
          ExplDangle( i, d-1, seq, seqlength);

        if( d >= i+5) {
          Qm[ pf_ij] += 
            //scale( S1 - S2 - S3)*
            Qm[ pf_index(i, d-1, seqlength)] *
            Qb[ pf_index( d, e, seqlength)] *
            EXP_FUNC( -(ALPHA_2 + ALPHA_3*(j-e) + bp_penalty)/
                 (kB*TEMP_K) )*
            ExplDangle( e+1, j, seq, seqlength);
        }
      }
    }
  }
}
#endif
