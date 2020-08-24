/** \file backtrack.c
 *  backtrack.c is part of the NUPACK software suite
 *  Copyright (c) 2007 Caltech. All rights reserved.
 *  Coded by: Robert Dirks 3/2006, Justin Bois 1/2007, Brian Wolfe 6/2011
 *
 *  The purpose of this collection of subroutines is to take the
 *  calculated values for the minimum energies on each contiguous
 *  subsequence and to find a structure that has the minimum possible
 *  free energy (mfe).
 *
 *  This code is definitely not optimal, but still runs fast.  Future
 *  versions of this software will hopefully clean it up.
 */

#include "backtrack.h"

//DBL_TYPE mfeEpsilon = 1.0; //allowed deviation from mfe

/* ********************************************* */

void bktrF_Fm_N3( int i, int j, int seq[], int seqlength, 
                 const DBL_TYPE *F, const DBL_TYPE *Fb, 
                 const DBL_TYPE *Fm, 
                 const DBL_TYPE *Fs, const DBL_TYPE *Fms, 
                 const int *nicks, 
                 int **etaN, dnaStructures *dnaStr, 
                 const char *type,
                 const int *maxILoopSize, const DBL_TYPE mfeEpsilon,
                 const int onlyOne){
  int d;//left base of rightmost base pair.
  int pf_ij = pf_index( i, j, seqlength);
  int pf_id1;

  DBL_TYPE extraTerms;
  DBL_TYPE tempMin;
  DBL_TYPE energy;

  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
#ifdef DEBUG
  printf("bktr %s %d %d ",  type, i, j);
#endif

  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);

  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;

#ifdef DEBUG
  printf("allowable = %.2Lf\n", (DBL_TYPE) allowableError);
  PrintS( &rootStr);
#endif

  energy = NickedEmptyF( i, j, nicks, seq, seqlength, etaN);
  if( strcmp( type, "F") == 0 && WithinEps( F[pf_ij], energy, ENERGY_TOLERANCE,
                                          allowableError) ) {
    //copyDnaStructures( &newStr, &rootStr);
    DBL_TYPE newErr = energy - F[pf_ij];
    addDnaStructures( dnaStr, &rootStr, newErr, mfeEpsilon, onlyOne); 
    //add new structures that are within mfeEpsilon of mfe

#ifdef DEBUG
    if( i == 0 && j == seqlength - 1) {
      printf("** %s add\n", type);
      PrintS( &rootStr);
      printf("** end\n");
    }
#endif

    matches++;
  }

  for( d = i; d <= j - 1; d++) {
   if( etaN[ EtaNIndex(d-0.5, d-0.5, seqlength)][0] == 0 || d == i ) {
     
     pf_id1 = pf_index(i,d-1,seqlength);
     
     tempMin = F[ pf_id1 ] + Fs[ pf_index( d, j, seqlength)];
     if( strcmp( type, "F") == 0 && 
        WithinEps( F[ pf_ij], tempMin, ENERGY_TOLERANCE, 
                  allowableError)) {
        
        copyDnaStructures( &newStr, &rootStr);
        bktrF_Fm_N3( i, d-1, seq, seqlength, F, Fb, Fm, Fs, 
                    Fms, nicks, etaN, &newStr, "F", maxILoopSize, mfeEpsilon, onlyOne);
        bktrFs_Fms( d, j, seq, seqlength, F, Fb, Fm, Fs, Fms, nicks, etaN, &newStr, "Fs",
                   maxILoopSize, mfeEpsilon, onlyOne);
        addDnaStructures( dnaStr, &newStr, tempMin - F[pf_ij], mfeEpsilon, onlyOne); 
        if( onlyOne) {
          clearDnaStructures( &rootStr);
          clearDnaStructures( &newStr);
          return;
        }
        
#ifdef DEBUG
        if( i == 0 && j == seqlength - 1) {
          printf("** %s add\n", type);
          PrintS( &newStr);
          printf("** end\n");
        }
#endif
        matches++;
      }

     extraTerms = DangleEnergy( i, d-1, seq, seqlength) +
       (ALPHA_3)*(d-i);

     if( etaN[ EtaNIndex( d-0.5, d-0.5, seqlength)][0] == 0) { 
       //otherwise Qm not possible

       if( etaN[ EtaNIndex(i+0.5, d-0.5, seqlength)][0] == 0 ) {
         tempMin = Fms[ pf_index( d, j, seqlength)] + extraTerms;
         if( strcmp( type, "Fm") == 0 && WithinEps( Fm[ pf_ij], tempMin, 
             ENERGY_TOLERANCE, allowableError)) {
           copyDnaStructures( &newStr, &rootStr);
           bktrFs_Fms( d, j, seq, seqlength, F, Fb, Fm, Fs, Fms, nicks, etaN, 
                      &newStr, "Fms", maxILoopSize, mfeEpsilon, onlyOne);
           addDnaStructures( dnaStr, &newStr, tempMin - Fm[pf_ij], mfeEpsilon, onlyOne);
           if( onlyOne) {
             clearDnaStructures( &rootStr);
             clearDnaStructures( &newStr);
             return;
           }

#ifdef DEBUG
           if( i == 0 && j == seqlength - 1) {
             printf("** %s add\n", type);
             PrintS( &newStr);
             printf("** end\n");
           }
#endif
           matches++;
         }
       }

       if( d >= i+2) {
         tempMin = Fm[ pf_index( i, d - 1, seqlength) ] +
           Fms[ pf_index( d, j, seqlength) ];

         if( strcmp( type, "Fm") == 0 && WithinEps( Fm[ pf_ij], tempMin, 
             ENERGY_TOLERANCE, allowableError)) {
           copyDnaStructures( &newStr, &rootStr);
           bktrF_Fm_N3( i, d-1, seq, seqlength, F, Fb, Fm, Fs, 
                       Fms, nicks, etaN, &newStr, "Fm",
                       maxILoopSize, mfeEpsilon, onlyOne);
           bktrFs_Fms( d, j, seq, seqlength, F, Fb, Fm, Fs, Fms, nicks, etaN, 
                       &newStr, "Fms", maxILoopSize, mfeEpsilon, onlyOne);
           addDnaStructures( dnaStr, &newStr, tempMin - Fm[pf_ij], mfeEpsilon, onlyOne);

           if( onlyOne) {
             clearDnaStructures( &rootStr);
             clearDnaStructures( &newStr);
             return;
           }

#ifdef DEBUG
           if( i == 0 && j == seqlength - 1) {
             printf("** %s add\n", type);
             PrintS( &newStr);
             printf("** end\n");
           }
#endif
           matches++;
         }
       }
     }
   }
  }

  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);

  if( matches == 0) {
   printf("Error in backtrack F_Fm_N3, %d %d\n", i, j);
   exit(1);
  }

  #ifdef DEBUG
  printf( "r%s %d %d\n", type, i, j );
  PrintS( dnaStr);
  #endif
 
}

/* **************** */

void bktrFs_Fms( int i, int j, int seq[], int seqlength,
                const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
                const DBL_TYPE *Fs, const DBL_TYPE *Fms, 
                const int *nicks, int **etaN, dnaStructures *dnaStr, const char *type,
                const int *maxILoopSize, const DBL_TYPE mfeEpsilon, 
                const int onlyOne) {
  
  int d; //base pair is i,d
  DBL_TYPE bp_penalty = 0.0;
  int pf_ij = pf_index( i, j, seqlength);
  
  DBL_TYPE extraTerms;
  int nNicks;
  int index_ij = EtaNIndex( i+0.5, j-0.5, seqlength);
  int start;
  DBL_TYPE tempMin;
  
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
  //printf("ae = %.2f\n", (DBL_TYPE) allowableError); 
#ifdef DEBUG
    printf("%s %d %d ", type, i, j);
#endif
  
  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  
#ifdef DEBUG
  printf("%.2lf\n", (double) mfeEpsilon);
  printf("%.2lf\n", (double) rootStr.minError);
  printf("%.2lf\n", (double) allowableError);
  PrintS( &rootStr);
#endif
  
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
         
         extraTerms = NickDangle( d+1,j,nicks, etaN,
                                 FALSE, seq,seqlength) + bp_penalty;    
         
         tempMin = Fb[ pf_index( i, d, seqlength) ] +
           extraTerms;
         
         if( strcmp( type, "Fs") == 0 && WithinEps( Fs[ pf_ij], tempMin, ENERGY_TOLERANCE,
                                                   allowableError) ) {
           matches++;
           copyDnaStructures( &newStr, &rootStr);
           bktrFb_N3(i, d, seq, seqlength, F, Fb, Fm, Fs, Fms, nicks, etaN, &newStr,
                     maxILoopSize, mfeEpsilon, onlyOne);
           addDnaStructures( dnaStr, &newStr, tempMin - Fs[pf_ij], mfeEpsilon, onlyOne);
        }

         // ******************** 
         
         extraTerms =  DangleEnergy( d+1, j, seq, seqlength)+
           bp_penalty + ALPHA_2 + ALPHA_3*(j-d);

         tempMin = Fb[ pf_index( i, d, seqlength) ] + 
           extraTerms;

         if( strcmp( type, "Fms") == 0 && WithinEps( Fms[ pf_ij], tempMin, ENERGY_TOLERANCE,
                                                    allowableError) ) {
            matches++;
            copyDnaStructures( &newStr, &rootStr);
            bktrFb_N3(i, d, seq, seqlength, F, Fb, Fm, Fs, Fms, nicks, etaN, &newStr,
                      maxILoopSize, mfeEpsilon, onlyOne);
            addDnaStructures( dnaStr, &newStr, tempMin - Fms[pf_ij], mfeEpsilon, onlyOne);
          }
       }
  }

  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);

  if( matches == 0) {
    printf("Error in backtrack Fs_Fms_N3, %d %d\n", i, j);
    exit(1);
  }

#ifdef DEBUG
  printf( "r%s %d %d\n", type, i, j );
  PrintS( dnaStr);
#endif
}

/* ************** */
void bktrFb_N3( int i, int j, int seq[], int seqlength, const DBL_TYPE *F, const DBL_TYPE *Fb, 
               const DBL_TYPE *Fm, const DBL_TYPE *Fs, const DBL_TYPE *Fms,
               const int *nicks, int **etaN, dnaStructures *dnaStr,
               const int *maxILoopSize, const DBL_TYPE mfeEpsilon,
               const int onlyOne) {
 
  int pf_ij = pf_index( i, j, seqlength);

  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  DBL_TYPE hairpinEnergy;
  int tmpMatch;

#ifdef DEBUG
  printf("%s %d %d %i %i ", "Fb", i, j, seq[i], seq[j]);	
#endif

  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  SetPairs( i, j, &rootStr);

  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  //printf("ae = %.2f\n", (DBL_TYPE) allowableError); 
  #ifdef DEBUG
  printf("%.2lf\n", (double) allowableError);
  PrintS( &rootStr);
  #endif

  hairpinEnergy = MinHairpin( i, j, seq, seqlength, etaN);
  if( WithinEps( Fb[ pf_ij], hairpinEnergy, ENERGY_TOLERANCE, allowableError) ) {
    //copyDnaStructures( &newStr, &rootStr);
    addDnaStructures( dnaStr, &rootStr, hairpinEnergy - Fb[pf_ij], mfeEpsilon, onlyOne);
    matches++;
  }

  tmpMatch = 0;
  copyDnaStructures( &newStr, &rootStr);
  if( bktrMinMultiloops(i, j, seq, seqlength,
                       F, Fb, Fm, Fs, Fms, nicks, 
                       etaN, &newStr, maxILoopSize, mfeEpsilon, onlyOne) ) {
    matches++;
    addDnaStructures( dnaStr, &newStr, 0, mfeEpsilon, onlyOne);
    tmpMatch = 1;

#ifdef DEBUG
    printf( "rml2%s %d %d\n", "Fb", i, j );
    PrintS( dnaStr);
#endif
  }

  if( tmpMatch == 1) {
   copyDnaStructures( &newStr, &rootStr);
   tmpMatch = 0;
  }

  if( bktrMinExteriorLoop( i, j, seq, seqlength,
                         F, Fb, Fm, Fs, Fms, nicks, 
                         etaN, &newStr, maxILoopSize, mfeEpsilon, onlyOne)) {
    matches++;
    addDnaStructures( dnaStr, &newStr, 0, mfeEpsilon, onlyOne);
    tmpMatch = 1;

#ifdef DEBUG
    printf( "rel2%s %d %d\n", "Fb", i, j );
    PrintS( dnaStr);
#endif
  }

  if( tmpMatch == 1) {
   copyDnaStructures( &newStr, &rootStr);
   tmpMatch = 0;
  }

  if( bktrMinInteriorLoop( i, j, seq, seqlength,
                         F, Fb, Fm, Fs, Fms, nicks, 
                         etaN, &newStr, maxILoopSize, mfeEpsilon, onlyOne)) {
    matches++;

    addDnaStructures( dnaStr, &newStr, 0, mfeEpsilon, onlyOne);
    tmpMatch = 1;
#ifdef DEBUG
    printf( "ril2%s %d %d %lf %lf\n", "Fb", i, j , (double)dnaStr->validStructs[0].error, (double)dnaStr->validStructs[0].correctedEnergy);
    PrintS( dnaStr);
#endif
  }

  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);

#ifdef DEBUG
  printf("Leaving check Fb %d %d\n", i, j);
#endif
   if( matches == 0) {
     printf("Error in backtrack Fb_N3, %d %d\n", i, j);
     exit(1);
   }

#ifdef DEBUG
  printf( "r Fb %d %d\n", i, j);
  PrintS( dnaStr);
#endif
 }

/* ************** */
int bktrMinMultiloops( int i, int j, int seq[], int seqlength,
                      const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
                      const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks, 
                      int **etaN, dnaStructures *dnaStr,
                      const int *maxILoopSize, const DBL_TYPE mfeEpsilon,
                      const int onlyOne){
  // Decomposes the region inside pair i,j into multiloops, i.e.
  // and excludes the possibility of "top level" nicks

  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE extraTerms;
  DBL_TYPE tempMin;
  int pf_ij = pf_index(i, j, seqlength);

  int d; // d is the left base of a rightmost paired base between i, j.

  dnaStructures rootStr = {NULL, 0, 0, 0, dnaStr->minError};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;

#ifdef DEBUG
  printf("%s %d %d\n", "Multiloop check", i, j);
#endif

  //copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  //clearDnaStructures( dnaStr);

  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  //printf("ae = %.2f\n", (DBL_TYPE) allowableError); 
  //printf("%.2f\n", (DBL_TYPE) allowableError);
  //PrintS( &rootStr);

  if( CanWCPair(seq[i], seq[j])) {
    for( d = i+3; d <= j - 2; d++) {
      //reset loop parameters
      bp_penalty = 0.0;

      if( etaN[ EtaNIndex( d-0.5, d-0.5, seqlength)][0] == 0 ) {
        
        if( seq[i] != BASE_C  && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }

        extraTerms = ( ALPHA_1 + ALPHA_2 + bp_penalty);
        
        tempMin = Fm[ pf_index( i+1, d-1, seqlength)] +
          Fms[ pf_index(d, j-1, seqlength)] + extraTerms;
        if( WithinEps( Fb[ pf_ij], tempMin, ENERGY_TOLERANCE, allowableError)) {
          if( matches++ == 0) {
            copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
            clearDnaStructures( dnaStr);
          }

          copyDnaStructures( &newStr, &rootStr);
          bktrF_Fm_N3( i+1, d-1, seq, seqlength,
                      F, Fb, Fm, Fs, 
                      Fms, nicks, etaN, &newStr, "Fm", maxILoopSize, mfeEpsilon,
                      onlyOne);
          bktrFs_Fms( d, j-1, seq, seqlength,
                     F, Fb, Fm, Fs, 
                     Fms, nicks, etaN, &newStr, "Fms", maxILoopSize, mfeEpsilon, onlyOne);
          addDnaStructures( dnaStr, &newStr, tempMin - Fb[pf_ij], mfeEpsilon, onlyOne);
        }
      }
    }
  }

  if( matches > 0) {
    clearDnaStructures( &rootStr);
    clearDnaStructures( &newStr);
  }

#ifdef DEBUG
  printf( "r Multiloop check %d %d\n", i, j);
  PrintS( dnaStr);
#endif

  return matches;  
}

/* ******** */
int bktrMinExteriorLoop( int i, int j, int seq[], int seqlength,
                        const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
                        const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks, 
                        int **etaN, dnaStructures *dnaStr,
                        const int *maxILoopSize, const DBL_TYPE mfeEpsilon,
                        const int onlyOne) {

  DBL_TYPE tempMin;
  DBL_TYPE bp_penalty = 0.0;
  int multiNick = -1;
  int index_ij;
  int leftIndex;
  int nNicks;
  int n;
  int iNicked, jNicked;
  int pf_ij = pf_index(i,j,seqlength);
  
  DBL_TYPE extraTerms;
  
  dnaStructures rootStr = {NULL, 0, 0, 0, dnaStr->minError};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
#ifdef DEBUG
  printf("%s %d %d\n", "Exterior Check", i, j);
#endif

  //Don't copy or clear until needed
  //copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  //clearDnaStructures( dnaStr);

  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  //printf("%.2f\n", (DBL_TYPE) allowableError);
  //PrintS( &rootStr);
  
  
  index_ij = EtaNIndex(i+0.5, j-0.5, seqlength);
  
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

      extraTerms = (bp_penalty);

      if( (iNicked == FALSE && jNicked == FALSE) ||
         (i == j - 1) || 
         (multiNick == i && jNicked == FALSE) ||
         (multiNick == j-1 && iNicked == FALSE ) ) {

        tempMin = F[ pf_index(i+1, multiNick, seqlength)] + 
         F[ pf_index( multiNick+1, j-1, seqlength)] + extraTerms;
        if( WithinEps( Fb[ pf_ij], tempMin, ENERGY_TOLERANCE, allowableError)) {

          if( matches++ == 0) {
           copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
           clearDnaStructures( dnaStr);
          }

          copyDnaStructures( &newStr, &rootStr);

          bktrF_Fm_N3( i+1, multiNick, seq, seqlength,
                      F, Fb, Fm, Fs, Fms, nicks, etaN, &newStr, "F", 
                      maxILoopSize, mfeEpsilon, onlyOne);
          bktrF_Fm_N3( multiNick+1, j-1, seq, seqlength,
                      F, Fb, Fm, Fs,  Fms, nicks, etaN, &newStr, "F", 
                      maxILoopSize, mfeEpsilon, onlyOne);
          addDnaStructures( dnaStr, &newStr, tempMin - Fb[pf_ij], mfeEpsilon, onlyOne);
        }
      }
    }
  }
  
  if( matches > 0) { 
    clearDnaStructures( &rootStr);
    clearDnaStructures( &newStr);
  }
  
#ifdef DEBUG
  printf( "r Exterior Check %d %d\n", i, j);
  PrintS( dnaStr);
#endif
  
  return matches;  
}


/* ********** */
int bktrMinInteriorLoop( int i, int j, int seq[], int seqlength,
                        const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
                        const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks, 
                        int **etaN, dnaStructures *dnaStr,
                        const int *maxILoopSize, const DBL_TYPE mfeEpsilon,
                        const int onlyOne) {
  int L, d, e;
  DBL_TYPE energy;
  int pf_ij = pf_index(i,j,seqlength);
  
  dnaStructures rootStr = {NULL, 0, 0, 0, dnaStr->minError};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;

#ifdef DEBUG
  printf("%s %d %d\n", "Interior loop check", i, j);
  PrintS( dnaStr);
#endif

  /*copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr); */

  if( j - i <= 1) {
#ifdef DEBUG
    printf( "r Interior Loop %d %d\n", i, j);
    PrintS( dnaStr);
#endif	
    return matches;
  }

  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;

  //Currently could be as bad as N^4?
  //for( L = j-i-2; L >= 1; L--) { //L = e-d}
  for( L = j-i-2; L >= j-i-2-maxILoopSize[pf_ij]; L--) { 
    for( d = i+1; d <= j-L-1; d++) {
      e = d+L;
      if( CanPair( seq[d], seq[e]) && etaN[ EtaNIndex( i+0.5, d-0.5, seqlength)][0] == 0 &&
         etaN[ EtaNIndex( e+0.5, j-0.5, seqlength)][0] == 0) {
        energy = InteriorEnergy( i,j,d,e,seq) + Fb[ pf_index(d,e,seqlength)];
        if( WithinEps( Fb[ pf_ij], energy, ENERGY_TOLERANCE, allowableError)) {

          if( matches++ == 0) {
           copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
           clearDnaStructures( dnaStr);
          };

          copyDnaStructures( &newStr, &rootStr);

          bktrFb_N3(d, e, seq, seqlength, F, Fb, Fm, Fs, Fms, nicks, etaN, &newStr, maxILoopSize, 
                   mfeEpsilon, onlyOne);

          addDnaStructures( dnaStr, &newStr, energy - Fb[pf_ij], mfeEpsilon, onlyOne);
        }
      }
    }
  }

  if( matches > 0) {
    clearDnaStructures( &rootStr); 
    clearDnaStructures( &newStr);
  }

#ifdef DEBUG
  printf( "r Interior Loop %d %d\n",i,j);
  PrintS( dnaStr);
#endif

  return matches;
}


/** ******************* */
void bktrF_Fm_FzN5( int i, int j, int seq[], int seqlength, 
                   const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                   const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
                   const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
                   const DBL_TYPE *Fgr, dnaStructures *dnaStr, const int *nicks, 
                   int **etaN, DBL_TYPE mfeEpsilon, const char *type) {

  int d, e; // d - e is internal basepair or pk boundary
  DBL_TYPE bp_penalty;
  int pf_de;
  int pf_ij = pf_index(i,j,seqlength);
  DBL_TYPE energy;

  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;

  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;

  energy = DangleEnergy(i, j, seq, seqlength);
  if( strcmp( "F", type) == 0 &&
    WithinEps( F[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) == TRUE) {
      matches++;
      addDnaStructures( dnaStr, &rootStr, energy - F[pf_ij], mfeEpsilon, 0);
  }
  else if( strcmp( "Fz", type) == 0 && i != 0 && j != seqlength - 1) { 
    energy = DangleEnergy(i, j, seq, seqlength) + BETA_3*(j-i+1);

    if( WithinEps( Fz[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) == TRUE) {
       matches++;
       addDnaStructures( dnaStr, &rootStr, energy - Fz[pf_ij], mfeEpsilon, 0);
    }
  }

  for( d = i; d <= j - 4; d++) {
    for( e = d + 4; e <= j; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE &&
        CanWCPair(seq[d], seq[e]) ) {
        bp_penalty = 0;
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }

        pf_de = pf_index( d, e, seqlength);

        energy = F[ pf_index(i, d-1, seqlength)] +
          Fb[ pf_de ] + bp_penalty +
          DangleEnergy( e+1, j, seq, seqlength);

        if( strcmp( "F", type) == 0 &&
           WithinEps( F[ pf_ij], energy, ENERGY_TOLERANCE,
                     allowableError) ) { 

          copyDnaStructures( &newStr, &rootStr);
          bktrF_Fm_FzN5( i, d-1, seq, seqlength,
                       F, Fb, Fm, Fp, Fz, Fg, 
                       Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                       etaN, mfeEpsilon, "F");

          bktrFbN5( d, e, seq, seqlength, 
                  F, Fb, Fm, Fp, Fz, Fg,
                  Fgls, Fgrs, Fgl, Fgr, 
                  &newStr, nicks,
                  etaN, mfeEpsilon);

          addDnaStructures( dnaStr, &newStr, energy - F[pf_ij], mfeEpsilon, 0);
          matches++;
        }

        if( i != 0 && j != seqlength - 1) {
          energy = (ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty) +
            Fb[ pf_de] + DangleEnergy( e+1, j, seq, seqlength) +
            DangleEnergy( i, d-1, seq, seqlength);

          if( strcmp( "Fm", type) == 0 &&
             WithinEps( Fm[ pf_ij], energy, ENERGY_TOLERANCE,
                       allowableError)== TRUE) {
            copyDnaStructures( &newStr, &rootStr);

            bktrFbN5(  d, e, seq, seqlength, F, Fb, Fm, Fp, Fz, Fg,
                    Fgls, Fgrs, Fgl, Fgr, 
                    &newStr, nicks,
                    etaN, mfeEpsilon);
            addDnaStructures( dnaStr, &newStr,  energy - Fm[pf_ij], mfeEpsilon, 0);
            matches++;
          }

          if( d >= i+5) {
            energy = Fm[ pf_index(i, d-1, seqlength)] +
              Fb[ pf_de] + (ALPHA_2 + ALPHA_3*(j-e) + bp_penalty) +
              DangleEnergy( e+1, j, seq, seqlength);

            if( strcmp( "Fm", type) == 0 &&
               WithinEps( Fm[ pf_ij], energy, ENERGY_TOLERANCE,
                         allowableError)== TRUE) {

              copyDnaStructures( &newStr, &rootStr);
              bktrFbN5(  d, e, seq, seqlength, F, Fb, Fm, 
                      Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                      &newStr, nicks,
                      etaN, mfeEpsilon);
              bktrF_Fm_FzN5( i,d-1,seq,seqlength, 
                           F, Fb, Fm, Fp, Fz, Fg, 
                           Fgls, Fgrs, Fgl, Fgr, 
                           &newStr, nicks,
                           etaN, mfeEpsilon, "Fm");
              addDnaStructures( dnaStr, &newStr,  energy - Fm[pf_ij], mfeEpsilon, 0);
              matches++;
            }
          }

          energy = Fz[ pf_index(i, d-1, seqlength)] + 
            Fb[ pf_de] + (BETA_2 + BETA_3*(j-e) + bp_penalty) +
            DangleEnergy( e+1, j, seq, seqlength);

          if( strcmp( "Fz", type) == 0 &&
             WithinEps( Fz[ pf_ij], energy, ENERGY_TOLERANCE,
                       allowableError) == TRUE) {

            copyDnaStructures( &newStr, &rootStr);
            bktrFbN5(  d, e, seq, seqlength, F, Fb, Fm, Fp, Fz, Fg,
                    Fgls, Fgrs, Fgl, Fgr, 
                    &newStr, nicks,
                    etaN, mfeEpsilon);

            bktrF_Fm_FzN5( i,d-1,seq,seqlength,
                         F, Fb, Fm, Fp, Fz, Fg, 
                         Fgls, Fgrs, Fgl, Fgr, 
                         &newStr, nicks,
                         etaN, mfeEpsilon, "Fz");

            addDnaStructures( dnaStr, &newStr,  energy - Fz[pf_ij], mfeEpsilon, 0);
            matches++;
          }
        }
      }
    }
  }

  for( d = i; d <= j - 5; d++) {
    for( e = d + 5; e <= j; e++) {

      pf_de = pf_index( d, e, seqlength);

      energy = F[ pf_index(i, d-1, seqlength)] +
       Fp[ pf_de ] + BETA_1 +
       DangleEnergy( e+1, j, seq, seqlength); 

      if( strcmp( "F", type) == 0 &&
        WithinEps( F[ pf_ij], energy, ENERGY_TOLERANCE,
                  allowableError)) {

        copyDnaStructures( &newStr, &rootStr);
        bktrF_Fm_FzN5( i, d-1, seq, seqlength,
                      F, Fb, Fm, Fp, Fz, Fg, 
                      Fgls, Fgrs, Fgl, Fgr,
                      &newStr, nicks,
                      etaN, mfeEpsilon, "F");

        bktrFpN5( d, e, seq, seqlength, F, Fb, Fm, Fp,
                 Fz, Fg, Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                 etaN, mfeEpsilon);
        addDnaStructures( dnaStr, &newStr, energy - F[pf_ij], mfeEpsilon, 0);
        matches++;
      }

      if( i != 0 && j != seqlength - 1) {
        energy =
         (BETA_1M + 2*ALPHA_2 + ALPHA_3*(d-i + j-e)) +
         Fp[ pf_de] +
         DangleEnergy( e+1, j, seq, seqlength) +
         DangleEnergy( i, d-1, seq, seqlength);

        if( strcmp( "Fm", type) == 0 &&
          WithinEps( Fm[ pf_ij], energy, ENERGY_TOLERANCE,
                    allowableError)) {

          copyDnaStructures( &newStr, &rootStr);
          bktrFpN5( d, e, seq, seqlength, F, Fb, Fm,
                   Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                   etaN, mfeEpsilon);
          addDnaStructures( dnaStr, &newStr, energy - Fm[pf_ij], mfeEpsilon, 0);
          matches++;
        }

        if( d >= i+5) {
          energy = Fm[ pf_index(i, d-1, seqlength)] +
           Fp[ pf_de] + (BETA_1M + 2*ALPHA_2 + ALPHA_3*(j-e)) +
           DangleEnergy( e+1, j, seq, seqlength);

          if( strcmp( "Fm", type) == 0 &&
            WithinEps( Fm[ pf_ij], energy, ENERGY_TOLERANCE,
                      allowableError) ) {

            copyDnaStructures( &newStr, &rootStr);
            bktrF_Fm_FzN5( i,d-1,seq,seqlength, 
                          F, Fb, Fm, Fp, Fz, Fg, 
                          Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                          etaN, mfeEpsilon, "Fm");
            bktrFpN5( d, e, seq, seqlength, F, Fb, Fm,
                     Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                     etaN, mfeEpsilon);
            
            addDnaStructures( dnaStr, &newStr, energy - Fm[pf_ij], mfeEpsilon, 0);
            matches++;
          }
        }

       energy = Fz[ pf_index(i, d-1, seqlength)] +
         Fp[ pf_de] + (BETA_1P + 2*BETA_2 + BETA_3*(j-e)) +
         DangleEnergy( e+1, j, seq, seqlength) ;
       
       if( strcmp( "Fz", type) == 0 &&
          WithinEps( Fz[ pf_ij], energy, ENERGY_TOLERANCE,
        allowableError)) {
          
          copyDnaStructures( &newStr, &rootStr);
          bktrF_Fm_FzN5( i,d-1,seq,seqlength,
                        F, Fb, Fm, Fp, Fz, Fg, 
                        Fgls, Fgrs, Fgl, Fgr,  &newStr, nicks,
                        etaN, mfeEpsilon, "Fz");
          bktrFpN5( d, e, seq, seqlength, F, Fb, Fm,
                   Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                   etaN, mfeEpsilon);
          
          addDnaStructures( dnaStr, &newStr, energy - Fz[pf_ij], mfeEpsilon, 0);
          matches++;
        }
      }
    }
  }

  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);

  if( matches == 0) {
   printf("Error in backtrack F_Fm_FzN5, %s: %d %d\n", type, i, j);
   exit(1);
  }
}


/* ***************************** */

void bktrFbN5( int i, int j, int seq[], int seqlength, 
              const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
              const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
              const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl, 
              const DBL_TYPE *Fgr, dnaStructures *dnaStr,
              const int *nicks, int **etaN, 
              const DBL_TYPE mfeEpsilon) {
  
  int d, e;
  int pf_de;
  DBL_TYPE energy;
  DBL_TYPE bp_penalty = 0;
  int pf_ij = pf_index( i, j, seqlength);
  
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
  if( CanPair( seq[i], seq[j]) == FALSE) {
    printf( "Error in pairing for bktrFbN5( %d %d)!!!\n", i, j);
    exit(1);
  }
  
  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  SetPairs( i, j, &rootStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  
  energy = MinHairpin( i, j, seq, seqlength, etaN);
  if( WithinEps( Fb[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
    //copyDnaStructures( &newStr, &rootStr);
    addDnaStructures( dnaStr, &rootStr, energy - Fb[pf_ij], mfeEpsilon, 0);
    matches++;		
    
  }
  
  //copyDnaStructures( &newStr, &rootStr);
  //next check multiloops and interior loops
  for( d = i+1; d <= j - 5; d++) {
    for( e = d + 4; e <= j - 1; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE) {
        bp_penalty = 0.0;
        
        pf_de = pf_index( d,e,seqlength);
        
        energy = InteriorEnergy( i, j, d, e, seq) +
          Fb[ pf_de ];
        
        if( WithinEps( Fb[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
          
          copyDnaStructures( &newStr, &rootStr);
          bktrFbN5( d, e, seq, seqlength,
                   F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                   Fgl, Fgr, &newStr, nicks,
                   etaN, mfeEpsilon);
          addDnaStructures( dnaStr, &newStr, energy - Fb[pf_ij], mfeEpsilon, 0);
          matches++;
        }
        
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }
        if( seq[i] != BASE_C && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }
        
        if( d>= i+6 && CanWCPair(seq[d], seq[e]) &&
           CanWCPair(seq[i], seq[j])) {
             
             energy = Fm[ pf_index(i+1, d-1, seqlength)] +
               Fb[ pf_de] +
               (ALPHA_1 + 2*ALPHA_2 + ALPHA_3*(j-e-1) + bp_penalty) +
               DangleEnergy( e+1, j-1, seq, seqlength);
             
             if( WithinEps( Fb[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
               
               copyDnaStructures( &newStr, &rootStr);
               
               bktrF_Fm_FzN5( i+1,d-1,seq,seqlength, 
                             F, Fb, Fm, Fp, Fz, Fg, 
                             Fgls, Fgrs, Fgl, Fgr, 
                             &newStr, nicks,
                             etaN, mfeEpsilon, 
                             "Fm");
               bktrFbN5( d, e, seq, seqlength,
                        F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                        Fgl, Fgr, &newStr, nicks,
                        etaN, mfeEpsilon);
               addDnaStructures( dnaStr, &newStr, energy - Fb[pf_ij], mfeEpsilon, 0);
               matches++;
             }
             
           }
      }
    }
  }
  
  //next check pseudoknots
  if( seq[i] != BASE_C  && seq[j] != BASE_C) {
    bp_penalty = AT_PENALTY;
  }
  
  for( d = i+1; d <= j - 6; d++) {
    for( e = d + 5; e <= j - 1; e++) {
      pf_de = pf_index( d, e, seqlength);
      
      energy = 
        (BETA_1M + ALPHA_1 + 3*ALPHA_2 + 
         (j-e-1 + d-i-1)*ALPHA_3 + bp_penalty) +
        DangleEnergy( i+1, d-1, seq, seqlength) +
        DangleEnergy( e+1, j-1, seq, seqlength) +
        Fp[ pf_de];
      
      
      if( WithinEps( Fb[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
        
        copyDnaStructures( &newStr, &rootStr);
        
        bktrFpN5( d, e, seq,seqlength,
                 F, Fb, Fm, Fp, Fz, Fg, 
                 Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                 etaN, mfeEpsilon);
        
        addDnaStructures( dnaStr, &newStr, energy - Fb[pf_ij], mfeEpsilon, 0);
        matches++;
      }
      
      energy = Fm[ pf_index( i+1, d-1, seqlength)] + 
        Fp[ pf_de] + ( BETA_1M + ALPHA_1 + 3*ALPHA_2 + bp_penalty + 
                      (j - e - 1)*ALPHA_3 ) +
        DangleEnergy( e+1, j-1, seq, seqlength);
      
      if( WithinEps( Fb[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
        copyDnaStructures( &newStr, &rootStr);
        
        bktrFpN5( d, e, seq,seqlength, 
                 F, Fb, Fm, Fp, Fz, Fg, 
                 Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                 etaN, mfeEpsilon);
        bktrF_Fm_FzN5( i+1,d-1,seq,seqlength, 
                      F, Fb, Fm, Fp, Fz, Fg, 
                      Fgls, Fgrs, Fgl, Fgr, 
                      &newStr, nicks,
                      etaN, mfeEpsilon, 
                      "Fm");
        
        addDnaStructures( dnaStr, &newStr, energy - Fb[pf_ij], mfeEpsilon, 0);
        matches++;
        
      }
    }
  }
  
  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);
  
  if( matches == 0) {
    printf("Error in backtrack Fb_N5, %d %d\n", i, j);
    exit(1);
  }
  
}
/* ***************** */

void bktrFgN5( int i, int d, int e, int j, int seq[], int seqlength, 
              const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
              const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
              const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
              const DBL_TYPE *Fgr,  
              dnaStructures *dnaStr,
              const int *nicks, int **etaN, 
              const DBL_TYPE mfeEpsilon) {
                
  int c, f;
  DBL_TYPE energy;
  DBL_TYPE bp_penalty, IJ_bp_penalty;
  int gap_idej = gap_index( i,d,e,j,seqlength);
  
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  
  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }
  
  if( i == d && e == j) {
    if( WithinEps( Fg[ gap_idej], 0, ENERGY_TOLERANCE, allowableError) ) {
      addDnaStructures( dnaStr, &rootStr, 0 - Fg[gap_idej], mfeEpsilon, 0);
      matches++;		
    }
    
    clearDnaStructures( &rootStr);
    clearDnaStructures( &newStr);
    
    return; //only possibility
  }
  
  //put the following after the single pair case 
  //(de already counted as interior pair)
  SetPairs( i, j, &rootStr);
  
  //Interior loop
  energy = InteriorEnergy(i,j,d,e,seq);
  if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
    addDnaStructures( dnaStr, &rootStr, energy - Fg[gap_idej], mfeEpsilon, 0);
    matches++;		
  }
  
  //new interior loop (O_N6)
  for( c = i+1; c <= d-1; c++) {
    for( f = e+1; f <= j-1; f++) {
      if( CanPair( seq[c], seq[f]) == TRUE) {
        
        energy = 
          InteriorEnergy( i, j, c, f, seq) + 
          Fg[ gap_index( c, d, e, f, seqlength)];
        
        if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
          copyDnaStructures( &newStr, &rootStr);
          bktrFgN5( c,d,e,f, seq, seqlength,
                   F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                   Fgl, Fgr, &newStr, nicks,
                   etaN, mfeEpsilon);
          
          
          addDnaStructures( dnaStr, &newStr, energy - Fg[ gap_idej], mfeEpsilon, 0);
          matches++;
          
        }
        
      }
    }
  }
  
  
  //Multiloops
  if( CanWCPair(seq[d], seq[e])) {
    
    bp_penalty = IJ_bp_penalty;
    if( seq[d] != BASE_C && seq[e] != BASE_C) {
      bp_penalty += AT_PENALTY;
    }
    
    //multiloop left
    energy = 
      Fm[ pf_index( i+1, d-1, seqlength)] +
      (ALPHA_1 + 2*ALPHA_2 + (j-e-1)*ALPHA_3 + bp_penalty) +
      DangleEnergy( e+1, j-1, seq, seqlength);
    
    if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
      copyDnaStructures( &newStr, &rootStr);
      
      bktrF_Fm_FzN5( i+1,d-1, seq, seqlength,
                    F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                    Fgl, Fgr, &newStr, nicks,
                    etaN, mfeEpsilon,  "Fm");
      
      addDnaStructures( dnaStr, &newStr, energy -  Fg[ gap_idej], mfeEpsilon, 0);
      matches++;
      
    }
    
    //multiloop right
    energy = 
      Fm[ pf_index( e+1, j-1, seqlength)] +
      (ALPHA_1 + 2*ALPHA_2 + (d-i-1)*ALPHA_3 + bp_penalty) +
      DangleEnergy( i+1, d-1, seq, seqlength) ;
    
    if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
      copyDnaStructures( &newStr, &rootStr);
      
      bktrF_Fm_FzN5( e+1,j-1, seq, seqlength,
                    F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                    Fgl, Fgr, &newStr, nicks,
                    etaN, mfeEpsilon,  "Fm");
      
      addDnaStructures( dnaStr, &newStr, energy -  Fg[ gap_idej], mfeEpsilon, 0);
      matches++;
      
    }
    
    
    //multiloop both
    energy = 
      Fm[ pf_index( i+1, d-1, seqlength)] +
      Fm[ pf_index( e+1, j-1, seqlength)] +
      (ALPHA_1 + 2*ALPHA_2 + bp_penalty);
    
    if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
      copyDnaStructures( &newStr, &rootStr);
      
      bktrF_Fm_FzN5( i+1,d-1, seq, seqlength,
                    F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                    Fgl, Fgr, &newStr, nicks,
                    etaN, mfeEpsilon,  "Fm");
      
      bktrF_Fm_FzN5( e+1,j-1, seq, seqlength,
                    F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                    Fgl, Fgr, &newStr, nicks,
                    etaN, mfeEpsilon,  "Fm");
      
      addDnaStructures( dnaStr, &newStr, energy -  Fg[ gap_idej], mfeEpsilon, 0);
      matches++;
      
    }
  }
  
  bp_penalty = IJ_bp_penalty;
  //Interior + multi left
  if( d > i+1) {
    for( f = e+1; f <= j-1; f++) {
      energy = 
        (ALPHA_1+ALPHA_2 + (j-f-1)*ALPHA_3+ bp_penalty) + 
        Fgls[ gap_index( i+1, d, e, f, seqlength)] +
        DangleEnergy( f+1, j-1, seq, seqlength);
      
      if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
        copyDnaStructures( &newStr, &rootStr);
        
        bktrFgls( i+1,d,e,f, seq, seqlength,
                 F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                 Fgl, Fgr, &newStr, nicks,
                 etaN, mfeEpsilon);
        
        addDnaStructures( dnaStr, &newStr, energy -  Fg[ gap_idej], mfeEpsilon, 0);
        matches++;
        
      }
    }
  }
  
  
  //Interior + multi right
  if( j-1 > e) {
    for( c = i+1; c <= d-1; c++) {
      
      energy = 
        (ALPHA_1+ALPHA_2 + (c-i-1)*ALPHA_3 + 
         bp_penalty) + 
        Fgrs[ gap_index( c, d, e, j-1, seqlength)] +
        DangleEnergy( i+1, c-1, seq, seqlength);
      if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
        copyDnaStructures( &newStr, &rootStr);
        
        bktrFgrs( c,d,e,j-1, seq, seqlength,
                 F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                 Fgl, Fgr, &newStr, nicks,
                 etaN, mfeEpsilon);
        
        addDnaStructures( dnaStr, &newStr, energy -  Fg[ gap_idej], mfeEpsilon, 0);
        matches++;
        
      }
    }	
    
    //Interior + multi both
    for( c = i+6; c <= d-1; c++) {	
      
      energy = 
        Fm[ pf_index( i+1, c-1, seqlength)] +
        (ALPHA_1+ALPHA_2+bp_penalty) +
        Fgrs[ gap_index( c, d, e, j-1, seqlength)];
      
      if( WithinEps( Fg[ gap_idej], energy, ENERGY_TOLERANCE, allowableError) ) {
        copyDnaStructures( &newStr, &rootStr);
        
        bktrFgrs( c,d,e,j-1, seq, seqlength,
                 F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                 Fgl, Fgr, &newStr, nicks,
                 etaN, mfeEpsilon);
        bktrF_Fm_FzN5( i+1, c-1, seq, seqlength,
                      F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                      Fgl, Fgr, &newStr, nicks,
                      etaN, mfeEpsilon, "Fm");
        
        
        addDnaStructures( dnaStr, &newStr, energy -  Fg[ gap_idej], mfeEpsilon, 0);
        matches++;
        
      }
      
    }
  }
  
  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);
  
  if( matches == 0) {
    printf("Error in backtrack Fg_N5, %d %d %d %d\n", i,d,e, j);
    exit(1);
  }
  
}
/* *************** */

void bktrFpN5( int i, int j, int seq[], int seqlength, 
              const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
              const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
              const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl, 
              const DBL_TYPE *Fgr, dnaStructures *dnaStr, const int *nicks, 
              int **etaN, DBL_TYPE mfeEpsilon) {
                
  int d, e, f;
  DBL_TYPE bp_penalty, new_bp_penalty;
  DBL_TYPE energy;
  
  int a, b, c;
  
  int pf_ij = pf_index(i,j,seqlength);
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  
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
                
                energy = 
                  Fg[ gap_index( i, a, d, e, seqlength) ]+
                  Fg[ gap_index( b, c, f, j, seqlength) ]+
                  (bp_penalty + 2*BETA_2) +
                  Fz[ pf_index( e + 1, f - 1, seqlength)]+
                  Fz[ pf_index( c + 1, d - 1, seqlength)]+
                  Fz[ pf_index( a + 1, b - 1, seqlength)];
                
                if( WithinEps( Fp[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
                  copyDnaStructures( &newStr, &rootStr);
                  
                  SetPairs( a, d, &newStr);
                  SetPairs( c, f, &newStr);
                  
                  bktrFgN5( i,a,d,e, seq, seqlength, 
                           F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                           &newStr, nicks,
                           etaN, mfeEpsilon);
                  bktrFgN5( b,c,f,j, seq, seqlength, 
                           F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                           &newStr, nicks,
                           etaN, mfeEpsilon);
                  bktrF_Fm_FzN5( c+1, d-1, seq, seqlength, 
                                F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgr,
                                Fgl, &newStr, nicks,
                                etaN, mfeEpsilon, "Fz");
                  bktrF_Fm_FzN5( e+1, f-1, seq, seqlength, 
                                F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgr,
                                Fgl, &newStr, nicks,
                                etaN, mfeEpsilon, "Fz");
                  bktrF_Fm_FzN5( a+1, b-1, seq, seqlength, 
                                F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgr,
                                Fgl, &newStr, nicks,
                                etaN, mfeEpsilon, "Fz");
                  
                  
                  addDnaStructures( dnaStr, &newStr, energy -  Fp[ pf_ij], mfeEpsilon, 0);
                  matches++;
                  
                }
              }
         }
       }
  }
  
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
                
                energy = Fg[ gap_index(i, i, e, f, seqlength)] +
                  Fz[ pf_index( i+1, d-1, seqlength)] +
                  Fgr[ gap_index(d, e-1, f+1, j, seqlength)] +
                  new_bp_penalty + BETA_2;
                
                
                
                if( WithinEps( Fp[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
                  copyDnaStructures( &newStr, &rootStr);
                  
                  SetPairs( i, e, &newStr);
                  
                  bktrFgN5( i,i,e,f, seq, seqlength, 
                           F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                           &newStr, nicks,
                           etaN, mfeEpsilon);
                  bktrF_Fm_FzN5( i+1,d-1,seq,seqlength, 
                                F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                                &newStr, nicks,
                                etaN, mfeEpsilon, "Fz");
                  bktrFgr( d, e-1, f+1, j, seq, seqlength,
                          F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                          &newStr, nicks,
                          etaN, mfeEpsilon);
                  
                  addDnaStructures( dnaStr, &newStr, energy -  Fp[ pf_ij], mfeEpsilon, 0);
                  matches++;
                  
                }
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
                  
                  energy = Fgl[ gap_index(i, d-1, e, f, seqlength)] +
                    Fg[ gap_index(d, d, j, j, seqlength)] +
                    Fz[ pf_index( d+1, e-1, seqlength)] +
                    Fz[ pf_index( f+1,j-1, seqlength)] +
                    new_bp_penalty+BETA_2;
                  
                  if( WithinEps( Fp[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
                    copyDnaStructures( &newStr, &rootStr);
                    
                    SetPairs( d, j, &newStr);
                    
                    bktrFgl( i,d-1,e,f,seq,seqlength, 
                            F, Fb, Fm, Fp, Fz, Fg, 
                            Fgls, Fgrs, Fgl, Fgr, &newStr, nicks,
                            etaN, mfeEpsilon);
                    bktrFgN5( d,d,j,j,seq,seqlength,
                             F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                             Fgl, Fgr, &newStr, nicks,
                             etaN, mfeEpsilon);
                    bktrF_Fm_FzN5( d+1,e-1,seq,seqlength,
                                  F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, 
                                  Fgl, Fgr, &newStr, nicks,
                                  etaN, mfeEpsilon, "Fz");
                    bktrF_Fm_FzN5( f+1,j-1,seq,seqlength,
                                  F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, 
                                  Fgl, Fgr, &newStr, nicks,
                                  etaN, mfeEpsilon, "Fz");
                    
                    addDnaStructures( dnaStr, &newStr, energy -  Fp[ pf_ij], mfeEpsilon, 0);
                    matches++;
                    
                  }
                  
                }
           }
         }
       }
  }
  
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
                  
                  
                  energy = Fgl[ gap_index(i, d-1, e, f, seqlength)] +
                    Fgr[ gap_index(d, e-1, f+1, j, seqlength)] +
                    (new_bp_penalty);
                  
                  
                  if( WithinEps( Fp[ pf_ij], energy, ENERGY_TOLERANCE, allowableError) ) {
                    copyDnaStructures( &newStr, &rootStr);
                    
                    bktrFgl( i,d-1,e,f,seq,seqlength,
                            F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                            &newStr, nicks,
                            etaN, mfeEpsilon);
                    bktrFgr( d,e-1,f+1,j,seq,seqlength,
                            F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                            &newStr, nicks,
                            etaN, mfeEpsilon);
                    
                    addDnaStructures( dnaStr, &newStr, energy -  Fp[ pf_ij], mfeEpsilon, 0);
                    matches++;
                    
                  }
                  
                }
           }
         }
       }
  }
  
  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);
  
  if( matches == 0) {
    printf("Error in bktrFpN5(%d %d).  No energy found!\n", i,j);
    exit(1);
  }
}

/* ********************** */
void bktrFgls( int i, int d, int e, int j, int seq[], int seqlength,
              const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
              const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
              const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl, 
              const DBL_TYPE *Fgr, 
              dnaStructures *dnaStr,
              const int *nicks, int **etaN, const DBL_TYPE mfeEpsilon) {
                
  int c;
  int pf_ic1;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE energy;
  
  DBL_TYPE Renergy = Fgls[ gap_index( i,d,e,j,seqlength)];
  
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  
  for( c = i + 5; c <= d-1; c++) {
    if( CanPair( seq[c], seq[j]) == TRUE 
       && CanWCPair(seq[c], seq[j]) 
       ) {
         
         bp_penalty = 0.0;
         if( seq[c] != BASE_C && seq[j] != BASE_C) {
           bp_penalty = AT_PENALTY;
         }
         
         pf_ic1 = pf_index( i, c-1, seqlength);
         
         energy =
           Fm[ pf_ic1] +
           Fg[ gap_index(c, d, e, j, seqlength)] +
           (ALPHA_2 + bp_penalty);
         
         if( WithinEps( Renergy, energy, ENERGY_TOLERANCE, allowableError) ) {
           copyDnaStructures( &newStr, &rootStr);
           
           bktrF_Fm_FzN5( i, c-1, seq,seqlength, 
                         F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                         Fgl, Fgr, &newStr, nicks,
                         etaN, mfeEpsilon, "Fm");
           bktrFgN5( c,d,e,j,seq,seqlength, 
                    F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs,
                    Fgl, Fgr, &newStr, nicks,
                    etaN, mfeEpsilon);
           
           addDnaStructures( dnaStr, &newStr, energy -  Renergy, mfeEpsilon, 0);
           matches++;
           
         }
         
       }
  }
  
  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);
  
  if( matches == 0) {
    printf("Error in bktrFgls(%d %d %d %d).  No energy found!\n", i,d,e,j);
    exit(1);
  }
  
}

/* ********************* */
void bktrFgrs( int i, int d, int e, int j, int seq[], int seqlength,
              const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
              const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
              const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl, const DBL_TYPE *Fgr, 
              dnaStructures *dnaStr,
              const int *nicks, int **etaN, const DBL_TYPE mfeEpsilon) {
                
  int f;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE energy;
  
  DBL_TYPE Renergy = Fgrs[ gap_index( i,d,e,j,seqlength)];
  
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;
  
  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;
  
  for( f = e+1; f <= j-5; f++) {
    if( CanPair( seq[i], seq[f]) == TRUE &&
       CanWCPair(seq[i], seq[f]) 
       ) {
         
         bp_penalty = 0;
         if( seq[i] != BASE_C && seq[f] != BASE_C) {
           bp_penalty = AT_PENALTY;
         }
         
         energy = 
           Fm[ pf_index( f+1, j, seqlength)] + 
           Fg[ gap_index( i, d, e, f, seqlength)] +
           (ALPHA_2 + bp_penalty);
         
         if( WithinEps( Renergy, energy, ENERGY_TOLERANCE, allowableError) ) {
           copyDnaStructures( &newStr, &rootStr);
           
           bktrF_Fm_FzN5( f+1, j, seq, seqlength, 
                         F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, 
                         Fgl, Fgr, &newStr, nicks,
                         etaN, mfeEpsilon, "Fm");
           bktrFgN5( i,d,e,f,seq,seqlength,
                    F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, 
                    Fgl, Fgr, &newStr, nicks,
                    etaN, mfeEpsilon);
           
           addDnaStructures( dnaStr, &newStr, energy -  Renergy, mfeEpsilon, 0);
           matches++;
           
         }
         
       }
  }
  
  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);
  
  if( matches == 0) {
    printf("Error in bktrFgrs(%d %d %d %d).  No energy found!\n", i,d,e,j);
    exit(1);
  }
}

/* *** */

void bktrFgl(  int i, int e, int f, int j, int seq[], int seqlength,
             const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
             const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
             const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
             const DBL_TYPE *Fgr, dnaStructures *dnaStr,
             const int *nicks, int **etaN, 
             const DBL_TYPE mfeEpsilon) {
               
  int d;
  int gap_idfj;
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE energy;

  DBL_TYPE Renergy = Fgl[ gap_index( i,e,f,j,seqlength)];
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;

  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;

  for( d = i+1; d <= MIN(f-4,e); d++) {
   if( CanPair( seq[d], seq[f]) == TRUE &&
      CanWCPair(seq[d], seq[f])
      ) {
        bp_penalty = 0.0;
        if( seq[d] != BASE_C && seq[f] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }
        
        gap_idfj = gap_index( i, d, f, j, seqlength);
        
        energy = 
          Fg[ gap_idfj]+Fz[ pf_index( d+1, e, seqlength)] +
          (BETA_2 + bp_penalty);
        
        if( WithinEps( Renergy, energy, ENERGY_TOLERANCE, allowableError) ) {
          copyDnaStructures( &newStr, &rootStr);
          
          SetPairs( d, f, &newStr);
          
          bktrFgN5( i, d, f, j, seq, seqlength, 
                   F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                   &newStr, nicks,
                   etaN, mfeEpsilon);
          bktrF_Fm_FzN5( d+1,e,seq,seqlength, 
                        F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                        &newStr, nicks,
                        etaN, mfeEpsilon, "Fz");
          
          addDnaStructures( dnaStr, &newStr, energy -  Renergy, mfeEpsilon, 0);
          matches++;
          
        }
      }
  }

  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);

  if( matches == 0) {
   printf("Error in bktrFgl(%d %d %d %d).  No energy found!\n", i,e,f,j);
   exit(1);
  }
 
}
/* ******************* */


void bktrFgr(  int i, int d, int e, int j, int seq[], int seqlength,
             const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm, 
             const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg, 
             const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
             const DBL_TYPE *Fgr, dnaStructures *dnaStr,
             const int *nicks, int **etaN, 
             const DBL_TYPE mfeEpsilon) {
 
  int f;
  DBL_TYPE energy;

  DBL_TYPE Renergy = Fgr[ gap_index( i,d,e,j,seqlength)];
  dnaStructures rootStr = {NULL, 0, 0, 0, NAD_INFINITY};
  dnaStructures newStr = {NULL, 0, 0, 0, NAD_INFINITY};
  DBL_TYPE allowableError;
  int matches = 0;

  copyDnaStructures( &rootStr, dnaStr);  //copy dnaStr to rootStr
  clearDnaStructures( dnaStr);
  allowableError = ENERGY_TOLERANCE+mfeEpsilon-rootStr.minError;

  for( f = e; f <= j-1; f++) {
   energy =
     Fgl[ gap_index(i,d,f,j,seqlength)]+
     Fz[ pf_index( e, f-1, seqlength)];
   
   if( WithinEps( Renergy, energy, ENERGY_TOLERANCE, allowableError) ) {
     copyDnaStructures( &newStr, &rootStr);
     
     bktrFgl( i,d,f,j,seq,seqlength,
             F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
             &newStr, nicks,
             etaN, mfeEpsilon);
     bktrF_Fm_FzN5( e,f-1,seq,seqlength, 
                   F, Fb, Fm, Fp, Fz, Fg, Fgls, Fgrs, Fgl, Fgr,
                   &newStr, nicks,
                   etaN, mfeEpsilon, "Fz");
     
     addDnaStructures( dnaStr, &newStr, energy -  Renergy, mfeEpsilon, 0);
     matches++;
     
   }
  }

  clearDnaStructures( &rootStr);
  clearDnaStructures( &newStr);

  if( matches == 0) {
   printf("Error in bktrFgr(%d %d %d %d).  No energy found!\n", i,d,e,j);
   exit(1);
  }
 
}

/* ***************************************************** */


void SetPair( int i, int j, int thepairs[]) {
  
  if( thepairs[i] != -1 || thepairs[j] != -1 ) {
    printf("Error in assigning pairs in backtrack %d %d!\n", i, j);
    exit(1);
  }
  thepairs[ i] = j;
  thepairs[ j] = i;
  
}

/* ******************** */

void SetPairs( int i, int j, dnaStructures *ds) {
  int n;
  int *thepairs;
  
  for( n = 0; n < ds->nStructs; n++) {
    thepairs = (ds->validStructs)[n].theStruct;
    SetPair( i, j, thepairs);
  }
}




/* ********************************** */
#ifdef O_N8
//EVERYTHING BELOW HERE IS NOT USED (NOR UPDATED!)

void BacktrackN8( int seqlength, char thefold[], int thepairs[], 
                 DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[], 
                 DBL_TYPE Fp[], DBL_TYPE Fz[], DBL_TYPE Fg[],
                 DBL_TYPE mfe, char seq[] ) {
                   
  // This program invokes the backtrack scheme for the O_N8

  bktrF_Fm_FzN8( 0, seqlength - 1, seq, seqlength, mfe, F, Fb, Fm, Fp,
               Fz, Fg, thepairs, "F");


  PrintStructure( thefold, thepairs, seqlength);

  printf( "mfe = %0.2f kcal/mol\n", (DBL_TYPE) mfe);
 
}

/** ******************* */
void bktrF_Fm_FzN8( int i, int j, char seq[], int seqlength, DBL_TYPE Renergy,
                   DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[],
                   DBL_TYPE Fp[], DBL_TYPE Fz[], DBL_TYPE Fg[], int thepairs[],
                   char *type) {
 
  int d, e; // d - e is internal basepair or pk boundary
  DBL_TYPE bp_penalty;
  int pf_de;
  DBL_TYPE energy;

  energy = DangleEnergy(i, j, seq, seqlength);
  if( strcmp( "F", type) == 0 &&
    Equal( energy, Renergy) == TRUE) 
   return;

  if( i != 0 && j != seqlength - 1) { 
   energy = DangleEnergy(i, j, seq, seqlength) +
     BETA_3*(j-i+1);
   
   if( strcmp( "Fz", type) == 0 &&
      Equal( energy, Renergy) == TRUE) 
     return;
  }  

  for( d = i; d <= j - 4; d++) {
   for( e = d + 4; e <= j; e++) {
     if( CanPair( seq[d], seq[e]) == TRUE &&
        CanWCPair(seq[d], seq[e]) ) {
          bp_penalty = 0;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty = AT_PENALTY;
          }
          
          pf_de = pf_index( d, e, seqlength);
          
          energy = F[ pf_index(i, d-1, seqlength)] +
            Fb[ pf_de ] + bp_penalty +
            DangleEnergy( e+1, j, seq, seqlength);
          
          if( strcmp( "F", type) == 0 &&
             Equal( energy, Renergy) == TRUE) { 
               
               bktrF_Fm_FzN8( i, d-1, seq, seqlength, F[ pf_index(i,d-1,seqlength)],
                             F, Fb, Fm, Fp, Fz, Fg, thepairs, "F");
               bktrFbN8( d, e, seq, seqlength, Fb[ pf_de], F, Fb, Fm, Fp, Fz, Fg,
                        thepairs);
               return;
             }
          
          if( i != 0 && j != seqlength - 1) {
            
            energy =
              (ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty) +
              Fb[ pf_de] +
              DangleEnergy( e+1, j, seq, seqlength) +
              DangleEnergy( i, d-1, seq, seqlength);
            if( strcmp( "Fm", type) == 0 &&
               Equal( energy, Renergy) == TRUE) {
                 bktrFbN8(  d, e, seq, seqlength, Fb[ pf_de], F, Fb, Fm, Fp, Fz, Fg,
                          thepairs);
                 return;
               }
            
            if( d >= i+5) {
              energy = Fm[ pf_index(i, d-1, seqlength)] +
                Fb[ pf_de] + (ALPHA_2 + ALPHA_3*(j-e) + bp_penalty) +
                DangleEnergy( e+1, j, seq, seqlength);
              
              if( strcmp( "Fm", type) == 0 &&
                 Equal( energy, Renergy) == TRUE) {
                   bktrFbN8(  d, e, seq, seqlength, Fb[ pf_de], F, Fb, Fm, 
                            Fp, Fz, Fg,
                            thepairs);
                   bktrF_Fm_FzN8( i,d-1,seq,seqlength, 
                                 Fm[ pf_index(i,d-1,seqlength)],
                                 F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
                   return;
                 }
            }
            
            
            energy = Fz[ pf_index(i, d-1, seqlength)] + 
              Fb[ pf_de] + (BETA_2 + BETA_3*(j-e) + bp_penalty) +
              DangleEnergy( e+1, j, seq, seqlength);
            
            if( strcmp( "Fz", type) == 0 &&
               Equal( energy, Renergy) == TRUE) {
                 bktrFbN8(  d, e, seq, seqlength, Fb[ pf_de], F, Fb, Fm, Fp, Fz, Fg,
                          thepairs);
                 bktrF_Fm_FzN8( i,d-1,seq,seqlength, Fz[ pf_index(i,d-1,seqlength)],
                               F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                 return;
               }
          }
        }
   }
  }

  for( d = i; d <= j - 5; d++) {
   for( e = d + 5; e <= j; e++) {
     
     pf_de = pf_index( d, e, seqlength);
     
     energy = F[ pf_index(i, d-1, seqlength)] +
       Fp[ pf_de ] + BETA_1 +
       DangleEnergy( e+1, j, seq, seqlength); 
     
     if( strcmp( "F", type) == 0 &&
        Equal( energy, Renergy) == TRUE) {
          bktrF_Fm_FzN8( i, d-1, seq, seqlength, F[ pf_index(i,d-1,seqlength)],
                        F, Fb, Fm, Fp, Fz, Fg, thepairs, "F");
          bktrFpN8( d, e, seq, seqlength, Fp[ pf_de], F, Fb, Fm, Fp,
                   Fz, Fg, thepairs);
          return;
        }
     
     if( i != 0 && j != seqlength - 1) {
       energy =
         (BETA_1M + 2*ALPHA_2 + ALPHA_3*(d-i + j-e)) +
         Fp[ pf_de] +
         DangleEnergy( e+1, j, seq, seqlength) +
         DangleEnergy( i, d-1, seq, seqlength);
       
       if( strcmp( "Fm", type) == 0 &&
          Equal( energy, Renergy) == TRUE) {
            bktrFpN8( d, e, seq, seqlength, Fp[ pf_de], F, Fb, Fm,
                     Fp, Fz, Fg, thepairs);
            return;
          }
       
       if( d >= i+5) {
         energy = Fm[ pf_index(i, d-1, seqlength)] +
           Fp[ pf_de] + (BETA_1M + 2*ALPHA_2 + ALPHA_3*(j-e)) +
           DangleEnergy( e+1, j, seq, seqlength);
         
         if( strcmp( "Fm", type) == 0 &&
            Equal( energy, Renergy) == TRUE) {
              bktrF_Fm_FzN8( i,d-1,seq,seqlength, Fm[ pf_index(i,d-1,seqlength)],
                            F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
              bktrFpN8( d, e, seq, seqlength, Fp[ pf_de], F, Fb, Fm,
                       Fp, Fz, Fg, thepairs);
              return;
            }
         
       }
       
       
       energy = Fz[ pf_index(i, d-1, seqlength)] +
         Fp[ pf_de] + (BETA_1P + 2*BETA_2 + BETA_3*(j-e)) +
         DangleEnergy( e+1, j, seq, seqlength) ;
       
       if( strcmp( "Fz", type) == 0 &&
          Equal( energy, Renergy) == TRUE) {
            bktrF_Fm_FzN8( i,d-1,seq,seqlength, Fz[ pf_index(i,d-1,seqlength)],
                          F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
            bktrFpN8( d, e, seq, seqlength, Fp[ pf_de], F, Fb, Fm,
                     Fp, Fz, Fg, thepairs);
            return;
          }
     }
   }    
  }
  printf("Error in bktrF_Fm_FzN8(%d %d).  No energy found!\n", i,j);
  exit(1);
  return;
}

/* ***************************** */
void bktrFbN8( int i, int j, char seq[], int seqlength, DBL_TYPE Renergy,
              DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[],
              DBL_TYPE Fp[], DBL_TYPE Fz[], DBL_TYPE Fg[], int thepairs[]) {
                
  int d, e;
  int pf_de;
  DBL_TYPE energy;
  DBL_TYPE bp_penalty;
  
  if( CanPair( seq[i], seq[j]) == FALSE) {
    printf( "Error in pairing for bktrFbN8( %d %d)!!!\n", i, j);
    exit(1);
  }
  
  SetPair( i, j, thepairs);
  
  if( Equal( HairpinEnergy( i, j, seq), Renergy) == TRUE) return;
  
  //next check multiloops and interior loops
  for( d = i+1; d <= j - 5; d++) {
    for( e = d + 4; e <= j - 1; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE) {
        bp_penalty = 0.0;
        
        pf_de = pf_index( d,e,seqlength);
        
        energy = InteriorEnergy( i, j, d, e, seq) +
          Fb[ pf_de ];
        
        if( Equal( energy, Renergy) == TRUE) {
          
          bktrFbN8( d, e, seq, seqlength, Fb[ pf_de],
                   F, Fb, Fm, Fp, Fz, Fg, thepairs);
          return;
        }
        
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }
        if( seq[i] != BASE_C && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }
        
        if( d>= i+6 && CanWCPair(seq[d], seq[e]) &&
           CanWCPair(seq[i], seq[j])) {
             
             energy = Fm[ pf_index(i+1, d-1, seqlength)] +
               Fb[ pf_de] +
               (ALPHA_1 + 2*ALPHA_2 + ALPHA_3*(j-e-1) + bp_penalty) +
               DangleEnergy( e+1, j-1, seq, seqlength);
             
             if( Equal( energy, Renergy) == TRUE) {
               bktrF_Fm_FzN8( i+1,d-1,seq,seqlength, 
                             Fm[pf_index(i+1,d-1,seqlength)],
                             F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
               
               bktrFbN8( d, e, seq, seqlength, Fb[ pf_de],
                        F, Fb, Fm, Fp, Fz, Fg, thepairs);
               return;
             }
           }
      }
    }
  }
  
  //next check pseudoknots
  if( seq[i] != BASE_C  && seq[j] != BASE_C) {
    bp_penalty = AT_PENALTY;
  }
  
  for( d = i+1; d <= j - 6; d++) {
    for( e = d + 5; e <= j - 1; e++) {
      pf_de = pf_index( d, e, seqlength);
      
      energy = 
        (BETA_1M + ALPHA_1 + 3*ALPHA_2 + 
         (j-e-1 + d-i-1)*ALPHA_3 + bp_penalty) +
        DangleEnergy( i+1, d-1, seq, seqlength) +
        DangleEnergy( e+1, j-1, seq, seqlength) +
        Fp[ pf_de];
      
      if( Equal( energy, Renergy) == TRUE) {
        bktrFpN8( d, e, seq,seqlength, Fp[ pf_de],
                 F, Fb, Fm, Fp, Fz, Fg, thepairs);
        return;
      }
      
      energy = Fm[ pf_index( i+1, d-1, seqlength)] + 
        Fp[ pf_de] + ( BETA_1M + ALPHA_1 + 3*ALPHA_2 + bp_penalty + 
                      (j - e - 1)*ALPHA_3 ) +
        DangleEnergy( e+1, j-1, seq, seqlength);
      
      if( Equal( energy, Renergy) == TRUE) {
        bktrF_Fm_FzN8( i+1,d-1,seq,seqlength, Fm[ pf_index(i+1,d-1,seqlength)],
                      F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
        bktrFpN8( d, e, seq,seqlength, Fp[ pf_de],
                 F, Fb, Fm, Fp, Fz, Fg, thepairs);
        return;
      }
    }
  }
  
  printf("Error in bktrFbN8(%d %d).  No energy found!\n", i,j);
  exit(1);
  return;
}

/* ***************** */
void bktrFpN8( int i, int j, char seq[], int seqlength, DBL_TYPE Renergy,
              DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[],
              DBL_TYPE Fp[], DBL_TYPE Fz[], DBL_TYPE Fg[], int thepairs[]) {
                
  int a, d, c, f, e, b; // 4 pairs: i-e, a-d, b-j, c-f
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE energy;
  
  FILE *pairfile;
  
  extern int containsPk;
  
#ifdef DEBUG
  printf("bktr_pk! %d %d\n", i, j);
#endif
    containsPk = TRUE;
  
  pairfile = fopen( "out.ene", "a");
  fprintf( pairfile, "p %d %d\n", i+1, j+1);
  fclose( pairfile); 
  
  
  if( j - i <= 4) {// not enough room for pk
    printf("Error in bktrFpN8(%d %d).  No energy found!\n", i,j);
    exit(1);
  }
  
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
                
                energy = 
                  Fg[ gap_index( i, a, d, e, seqlength) ]+
                  Fg[ gap_index( b, c, f, j, seqlength) ]+
                  bp_penalty + 2*BETA_2  +
                  Fz[ pf_index( e + 1, f - 1, seqlength)]+
                  Fz[ pf_index( c + 1, d - 1, seqlength)]+
                  Fz[ pf_index( a + 1, b - 1, seqlength)];
                
                if( Equal( energy, Renergy) == TRUE) {
                  SetPair( a, d, thepairs);
                  SetPair( c, f, thepairs);
                  
                  bktrFgN8( i, a, d, e, seq, seqlength,
                           Fg[ gap_index(i,a,d,e,seqlength)],
                           F, Fb, Fm, Fp, Fz, Fg,
                           thepairs);
                  bktrFgN8( b, c, f, j, seq, seqlength,
                           Fg[ gap_index(b,c,f,j,seqlength)],
                           F, Fb, Fm, Fp, Fz, Fg,
                           thepairs);
                  bktrF_Fm_FzN8( e+1,f-1,seq,seqlength,
                                Fz[ pf_index(e+1,f-1,seqlength)],
                                F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                  bktrF_Fm_FzN8( c+1,d-1,seq,seqlength,
                                Fz[ pf_index(c+1,d-1,seqlength)],
                                F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                  bktrF_Fm_FzN8( a+1,b-1,seq,seqlength,
                                Fz[ pf_index(a+1,b-1,seqlength)],
                                F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                  return;
                }
              }
         }
       }
  }
  
  if( j - i <= 6) {// not enough room for pk
    printf("Error in bktrFpN8(%d %d).  No energy found!\n", i,j);
    exit(1);
  }
  
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
                              
                              energy = 
                                Fg[ gap_index( i, a, d, e, seqlength) ]+
                                Fg[ gap_index( b, c, f, j, seqlength) ]+
                                bp_penalty + 2*BETA_2 +
                                Fz[ pf_index( e + 1, f - 1, seqlength)]+
                                Fz[ pf_index( c + 1, d - 1, seqlength)]+
                                Fz[ pf_index( a + 1, b - 1, seqlength)];
                              
                              if( Equal( energy, Renergy) == TRUE) {
                                SetPair( a, d, thepairs);
                                SetPair( c, f, thepairs);
                                
                                bktrFgN8( i, a, d, e, seq, seqlength,
                                         Fg[ gap_index(i,a,d,e,seqlength)],
                                         F, Fb, Fm, Fp, Fz, Fg,
                                         thepairs);
                                bktrFgN8( b, c, f, j, seq, seqlength,
                                         Fg[ gap_index(b,c,f,j,seqlength)],
                                         F, Fb, Fm, Fp, Fz, Fg,
                                         thepairs);
                                bktrF_Fm_FzN8( e+1,f-1,seq,seqlength,
                                              Fz[ pf_index(e+1,f-1,seqlength)],
                                              F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                bktrF_Fm_FzN8( c+1,d-1,seq,seqlength,
                                              Fz[ pf_index(c+1,d-1,seqlength)],
                                              F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                bktrF_Fm_FzN8( a+1,b-1,seq,seqlength,
                                              Fz[ pf_index(a+1,b-1,seqlength)],
                                              F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                return;
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
                              
                              energy = 
                                Fg[ gap_index( i, a, d, e, seqlength) ]+
                                Fg[ gap_index( b, c, f, j, seqlength) ]+
                                bp_penalty + 2*BETA_2 +
                                Fz[ pf_index( e + 1, f - 1, seqlength)]+
                                Fz[ pf_index( c + 1, d - 1, seqlength)]+
                                Fz[ pf_index( a + 1, b - 1, seqlength)];
                              
                              if( Equal( energy, Renergy) == TRUE) {
                                SetPair( a, d, thepairs);
                                SetPair( c, f, thepairs);
                                
                                bktrFgN8( i, a, d, e, seq, seqlength,
                                         Fg[ gap_index(i,a,d,e,seqlength)],
                                         F, Fb, Fm, Fp, Fz, Fg,
                                         thepairs);
                                bktrFgN8( b, c, f, j, seq, seqlength,
                                         Fg[ gap_index(b,c,f,j,seqlength)],
                                         F, Fb, Fm, Fp, Fz, Fg,
                                         thepairs);
                                bktrF_Fm_FzN8( e+1,f-1,seq,seqlength,
                                              Fz[ pf_index(e+1,f-1,seqlength)],
                                              F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                bktrF_Fm_FzN8( c+1,d-1,seq,seqlength,
                                              Fz[ pf_index(c+1,d-1,seqlength)],
                                              F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                bktrF_Fm_FzN8( a+1,b-1,seq,seqlength,
                                              Fz[ pf_index(a+1,b-1,seqlength)],
                                              F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                return;
                              }
                              
                            }
                       }
                     }
                }
              }
         }
    }
  }
  
  
  if( j - i <= 7) {// not enough room for pk
    printf("Error in bktrFpN8(%d %d).  No energy found!\n", i,j);
    exit(1);
  }
  
  //at least 2 pairs per Qg
  for( a = i+1; a <= j - 7; a++) {
    for( b = a+1; b <= j - 6; b++) {
      if( CanPair( seq[ b], seq[ j]) == TRUE	  
         && CanWCPair(seq[b], seq[j])) {
           for( c = b+1; c <= j - 5; c++) {
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
                                  
                                  energy = 
                                    Fg[ gap_index( i, a, d, e, seqlength) ]+
                                    Fg[ gap_index( b, c, f, j, seqlength) ]+
                                    bp_penalty + 2*BETA_2 +
                                    Fz[ pf_index( e + 1, f - 1, seqlength)]+
                                    Fz[ pf_index( c + 1, d - 1, seqlength)]+
                                    Fz[ pf_index( a + 1, b - 1, seqlength)];
                                  
                                  if( Equal( energy, Renergy) == TRUE) {
                                    SetPair( a, d, thepairs);
                                    SetPair( c, f, thepairs);
                                    
                                    bktrFgN8( i, a, d, e, seq, seqlength,
                                             Fg[ gap_index(i,a,d,e,seqlength)],
                                             F, Fb, Fm, Fp, Fz, Fg,
                                             thepairs);
                                    bktrFgN8( b, c, f, j, seq, seqlength,
                                             Fg[ gap_index(b,c,f,j,seqlength)],
                                             F, Fb, Fm, Fp, Fz, Fg,
                                             thepairs);
                                    bktrF_Fm_FzN8( e+1,f-1,seq,seqlength,
                                                  Fz[ pf_index(e+1,f-1,seqlength)],
                                                  F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                    bktrF_Fm_FzN8( c+1,d-1,seq,seqlength,
                                                  Fz[ pf_index(c+1,d-1,seqlength)],
                                                  F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                    bktrF_Fm_FzN8( a+1,b-1,seq,seqlength,
                                                  Fz[ pf_index(a+1,b-1,seqlength)],
                                                  F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fz");
                                    return;
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
  
  printf("Error in bktrFpN8(%d %d).  No energy found!\n", i,j);
  exit(1);
  return;
}

/* ******************************** */
void bktrFgN8( int i, int d, int e, int j, char seq[], int seqlength, 
              DBL_TYPE Renergy,
              DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[],
              DBL_TYPE Fp[], DBL_TYPE Fz[], DBL_TYPE Fg[], int thepairs[]) {
                
  int c, f;
  DBL_TYPE energy;
  DBL_TYPE bp_penalty, IJ_bp_penalty;
  
  
  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }
  
  if( i == d && e == j) {
    if( Equal( 0, Renergy) == TRUE) return;
    printf("Error in Qg backtrack %d %d!!\n", i, d);
  }
  
  
  //put the following after the single pair case 
  //(already counted as interior pair)
  SetPair( i, j, thepairs);
  
  //Interior loop
  if( Equal( InteriorEnergy( i,j,d,e,seq), Renergy) == TRUE ) {
    return;
  }
  
  //new interior loop
  for( c = i+1; c <= d-1; c++) {
    for( f = e+1; f <= j-1; f++) {
      if( CanPair( seq[c], seq[f]) == TRUE) {
        
        energy = 
          InteriorEnergy( i, j, c, f, seq) + 
          Fg[ gap_index( c, d, e, f, seqlength)];
        
        if( Equal( energy, Renergy) == TRUE) {
          bktrFgN8( c,d,e,f, seq, seqlength,
                   Fg[ gap_index(c,d,e,f,seqlength)],
                   F, Fb, Fm, Fp, Fz, Fg,
                   thepairs);
          return;
        }
      }
    }
  }
  
  //multiloops
  if( CanWCPair(seq[i], seq[j])) {
    //multiloop left
    if( CanWCPair(seq[d], seq[e])) {
      bp_penalty = IJ_bp_penalty;
      if( seq[d] != BASE_C && seq[e] != BASE_C) {
        bp_penalty += AT_PENALTY;
      }
      
      energy = 
        Fm[ pf_index( i+1, d-1, seqlength)] +
        (ALPHA_1 + 2*ALPHA_2 + (j-e-1)*ALPHA_3 + bp_penalty) +
        DangleEnergy( e+1, j-1, seq, seqlength);
      
      if( Equal( energy, Renergy) == TRUE) {
        bktrF_Fm_FzN8( i+1,d-1,seq,seqlength, 
                      Fm[ pf_index(i+1,d-1,seqlength)],
                      F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
        return;
      }
      
      ///multiloop right
      energy = 
        Fm[ pf_index( e+1, j-1, seqlength)] +
        (ALPHA_1 + 2*ALPHA_2 + (d-i-1)*ALPHA_3 + bp_penalty) +
        DangleEnergy( i+1, d-1, seq, seqlength);
      
      if( Equal( energy, Renergy) == TRUE) {
        bktrF_Fm_FzN8( e+1,j-1,seq,seqlength, 
                      Fm[ pf_index(e+1,j-1,seqlength)],
                      F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
        return;
      }
      
      //multiloop both
      energy = 
        Fm[ pf_index( i+1, d-1, seqlength)] +
        Fm[ pf_index( e+1, j-1, seqlength)] +
        (ALPHA_1 + 2*ALPHA_2 + bp_penalty);
      if( Equal( energy, Renergy) == TRUE) {
        bktrF_Fm_FzN8( i+1,d-1,seq,seqlength, 
                      Fm[ pf_index(i+1,d-1,seqlength)],
                      F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
        bktrF_Fm_FzN8( e+1,j-1,seq,seqlength, 
                      Fm[ pf_index(e+1,j-1,seqlength)],
                      F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
        return;
      }
    }
    
    //Multiloop Left + More Qg
    for( c = i+6; c <= d-1; c++) {
      for( f = e+1; f <= j-1; f++) {
        if( CanPair( seq[c], seq[f]) == TRUE &&
           CanWCPair(seq[c], seq[f])) {
             
             bp_penalty = IJ_bp_penalty;
             if( seq[c] != BASE_C && seq[f] != BASE_C) {
               bp_penalty += AT_PENALTY;
             }
             
             energy = 
               Fm[ pf_index( i+1, c-1, seqlength)] +
               (ALPHA_1+2*ALPHA_2+(j-f-1)*ALPHA_3 + bp_penalty) + 
               Fg[ gap_index( c, d, e, f, seqlength)] +
               DangleEnergy( f+1, j-1, seq, seqlength);
             if( Equal( energy, Renergy) == TRUE) {
               bktrF_Fm_FzN8( i+1,c-1,seq,seqlength, 
                             Fm[ pf_index(i+1,c-1,seqlength)],
                             F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
               bktrFgN8( c,d,e,f, seq, seqlength,
                        Fg[ gap_index(c,d,e,f,seqlength)],
                        F, Fb, Fm, Fp, Fz, Fg,
                        thepairs);
               
               return;
             }
           }
      }
    }
    
    
    //multiloop right + Qg
    for( c = i+1; c <= d-1; c++) {
      for( f = e+1; f <= j-6; f++) {
        if( CanPair( seq[c], seq[f]) == TRUE &&
           CanWCPair(seq[c], seq[f])) {
             
             bp_penalty = IJ_bp_penalty;
             if( seq[c] != BASE_C && seq[f] != BASE_C) {
               bp_penalty += AT_PENALTY;
             } 
             
             energy = 
               Fm[ pf_index( f+1, j-1, seqlength)] +
               (ALPHA_1+2*ALPHA_2+(c-i-1)*ALPHA_3+bp_penalty) + 
               Fg[ gap_index( c, d, e, f, seqlength)] +
               DangleEnergy( i+1, c-1, seq, seqlength);
             
             if( Equal( energy, Renergy) == TRUE) {
               bktrF_Fm_FzN8( f+1,j-1,seq,seqlength, 
                             Fm[ pf_index(f+1,j-1,seqlength)],
                             F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
               bktrFgN8( c,d,e,f, seq, seqlength,
                        Fg[ gap_index(c,d,e,f,seqlength)],
                        F, Fb, Fm, Fp, Fz, Fg,
                        thepairs);
               
               return;
             }
           }
      }
    }
    
    //multiloop both + Qg
    for( c = i+6; c <= d-1; c++) {
      for( f = e+1; f <= j-6; f++) {
        if( CanPair( seq[c], seq[f]) == TRUE &&
           CanWCPair(seq[c], seq[f])) {
             
             bp_penalty = IJ_bp_penalty;
             if( seq[c] != BASE_C && seq[f] != BASE_C) {
               bp_penalty += AT_PENALTY;
             }
             
             energy = 
               Fm[ pf_index( i+1, c-1, seqlength)] +
               Fm[ pf_index( f+1, j-1, seqlength)] +
               (ALPHA_1+2*ALPHA_2+bp_penalty) +
               Fg[ gap_index( c, d, e, f, seqlength)];
             if( Equal( energy, Renergy) == TRUE) {
               bktrF_Fm_FzN8( f+1,j-1,seq,seqlength, 
                             Fm[ pf_index(f+1,j-1,seqlength)],
                             F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
               bktrF_Fm_FzN8( i+1,c-1,seq,seqlength, 
                             Fm[ pf_index(i+1,c-1,seqlength)],
                             F, Fb, Fm, Fp, Fz, Fg, thepairs, "Fm");
               bktrFgN8( c,d,e,f, seq, seqlength,
                        Fg[ gap_index(c,d,e,f,seqlength)],
                        F, Fb, Fm, Fp, Fz, Fg,
                        thepairs);
               
               return;
             }	  
           }
      }
    }
  }
  
  
  printf("Error in bktrFgN8(%d %d).  No energy found!\n", i,j);
  exit(1);
  return;
}
#endif



#ifdef O_N4

void BacktrackN4( int seqlength, char thefold[], int thepairs[], 
                 DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[], 
                 DBL_TYPE mfe, char seq[] ) {
  // This program invokes the backtrack scheme

  bktrF_Fm_N4( 0, seqlength - 1, seq, seqlength, mfe, F, Fb, Fm, 
             thepairs, "F");

  PrintStructure( thefold, thepairs, seqlength);

  printf( "mfe = %0.2f kcal/mol\n", (DBL_TYPE) mfe);
}

/*  ***************************** */
void bktrF_Fm_N4( int i, int j, char seq[], int seqlength, DBL_TYPE Renergy, 
                 DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[], 
                 int thepairs[], char *type){
 
  int d, e; // d - e is internal basepair 
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE energy;

  if( strcmp( type, "F") == 0 && 
    Equal( DangleEnergy(i,j,seq,seqlength), Renergy) == TRUE) return; 

  for( d = i; d <= j - 4; d++) {
   for( e = d + 4; e <= j; e++) {
     if( CanPair( seq[d], seq[e]) == TRUE &&
        CanWCPair(seq[d], seq[e]) ) {
          bp_penalty = 0;
          if( seq[d] != BASE_C && seq[e] != BASE_C) {
            bp_penalty = AT_PENALTY;
          }
          
          energy = F[ pf_index(i, d-1, seqlength)] +
            Fb[ pf_index( d, e, seqlength) ] +
            bp_penalty + 
            DangleEnergy( e+1, j, seq, seqlength); 
          if( strcmp( type, "F") == 0 &&
             Equal( energy, Renergy) == TRUE) {
               
               bktrFbN4( d, e, seq, seqlength, Fb[ pf_index(d,e,seqlength)],  
                        F, Fb, Fm, thepairs);
               
               bktrF_Fm_N4( i, d-1, seq, seqlength, F[ pf_index(i,d-1,seqlength)],
                           F, Fb, Fm, thepairs, "F");
               
               return;
             }
          
          
          energy =
            (ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty) +
            Fb[ pf_index( d, e, seqlength)] +
            DangleEnergy( e+1, j, seq, seqlength) +
            DangleEnergy( i, d-1, seq, seqlength);
          
          if( strcmp( type, "Fm") == 0 && Equal( energy, Renergy) == TRUE) {
            
            bktrFbN4( d, e, seq, seqlength, Fb[ pf_index(d,e,seqlength)],
                     F, Fb, Fm, thepairs);
            return;
          }
          
          
          if( d >= i+5) {
            
            energy = Fm[ pf_index(i, d-1, seqlength)] +
              Fb[ pf_index( d, e, seqlength)] +
              (ALPHA_2 + ALPHA_3*(j-e) + bp_penalty) +
              DangleEnergy( e+1, j, seq, seqlength);
            
            if( strcmp( type, "Fm") == 0 && Equal( energy, Renergy) == TRUE) {
              
              bktrF_Fm_N4( i, d-1, seq, seqlength, 
                          F[ pf_index(i,d-1,seqlength)],
                          F, Fb, Fm, thepairs, "Fm");
              
              
              
              bktrFbN4( d, e, seq, seqlength, Fb[ pf_index(d,e,seqlength)],
                       F, Fb, Fm, thepairs);
              return;
            }
          }
        }
   }
  }

  printf("Unable to find mfe structure in bktrF_Fm_N4( %d, %d)!!!\n",
        i,j);
  exit(1);
  return;
 
}

/* ************* */
void bktrFbN4( int i, int j, char seq[], int seqlength, DBL_TYPE Renergy,
              DBL_TYPE F[], DBL_TYPE Fb[], DBL_TYPE Fm[], int thepairs[]) {
                
  int d, e;
  DBL_TYPE energy;
  DBL_TYPE bp_penalty;
  
  if( CanPair( seq[i], seq[j]) == FALSE) {
    printf( "Error in pairing for bktrFbN4( %d %d)!!!\n", i, j);
    exit(1);
  }
  
  SetPair( i, j, thepairs);
  
  if( Equal( HairpinEnergy( i, j, seq), Renergy) == TRUE) return;
  
  for( d = i+1; d <= j - 5; d++) {
    for( e = d + 4; e <= j - 1; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE) {
        bp_penalty = 0.0;
        
        energy = InteriorEnergy( i, j, d, e, seq) +
          Fb[ pf_index( d, e, seqlength) ];
        
        if( Equal( energy, Renergy) == TRUE) {
          
          bktrFbN4( d, e, seq, seqlength, Fb[ pf_index(d,e,seqlength)],
                   F, Fb, Fm, thepairs);
          return;
        }
        
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }
        if( seq[i] != BASE_C && seq[j] != BASE_C) {
          bp_penalty += AT_PENALTY;
        }
        
        if( d>= i+6 && CanWCPair(seq[d], seq[e]) &&
           CanWCPair(seq[i], seq[j])) {
             
             energy = Fm[ pf_index(i+1, d-1, seqlength)] +
               Fb[ pf_index( d, e, seqlength)] +
               (ALPHA_1 + 2*ALPHA_2 + ALPHA_3*(j-e-1) + bp_penalty) +
               DangleEnergy( e+1, j-1, seq, seqlength);
             
             if( Equal( energy, Renergy) == TRUE) {
               bktrF_Fm_N4( i+1,d-1,seq,seqlength, Fm[pf_index(i+1,d-1,seqlength)],
                           F, Fb, Fm, thepairs, "Fm");
               
               bktrFbN4( d, e, seq, seqlength, Fb[ pf_index(d,e,seqlength)],
                        F, Fb, Fm, thepairs);
               return;
             }
           }
      }
    }
  }
  
  printf("Error in bktrFbN4(%d %d).  No energy found!\n", i,j);
  exit(1);
  return;
}


#endif
/* ************************ */

