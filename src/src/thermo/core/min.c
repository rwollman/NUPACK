/*
  min.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 8/2001, Justin Bois 1/2007

  This collection of functions to calculate the energy of the minimum
  fold on a subsequence of a DNA sequence, given various restraints.
*/

#include "min.h"

/* ***************************** */ 

DBL_TYPE MinHairpin( int i, int j, int seq[], int seqlength, int **etaN) {
	
  //this version disallows nicks here
	
  int index;
  int nNicks = 0;
	
	
  index = EtaNIndex(i+0.5, j-0.5, seqlength);
  nNicks = etaN[ index][0];
	
	
  if( nNicks >= 1) return NAD_INFINITY;
	
  if( nNicks == 0 && j-i <= 3) {
    return NAD_INFINITY;
  }
  
  return HairpinEnergy( i, j, seq);	

}

/* ******************* */



DBL_TYPE MinMultiloops( int i, int j, int seq[], 
			DBL_TYPE *Fms, DBL_TYPE *Fm, int seqlength,
			int **etaN){
  // Decomposes the region inside pair i,j into multiloops, i.e.
  // and excludes the possibility of "top level" nicks
	
  DBL_TYPE min_energy = 0.0;
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE extraTerms;
  DBL_TYPE tempMin;
	
  int d; // d is the left base of a rightmost paired base between i, j.

  min_energy = NAD_INFINITY;
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
	min_energy = MIN( tempMin, min_energy);
				
      }
			
    }
  }
  
	
  return min_energy;  
}

/* ********* */
DBL_TYPE MinExteriorLoop( int i,int j, int seq[], int seqlength, 
			  DBL_TYPE *F, int *nicks, int **etaN) {
	
  DBL_TYPE min_energy = NAD_INFINITY;
  DBL_TYPE tempMin;
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
				
	min_energy = MIN( tempMin, min_energy);
				
				
      }
			
			
    }
  }
	
  return min_energy;  
}

/* ****************** */

void MinFastILoops( int i, int j, int L, int seqlength, int seq[],
		    int **etaN, DBL_TYPE *Fb, DBL_TYPE *Fx, DBL_TYPE *Fx_2,
		    DBL_TYPE *minILoopEnergyBySize) {
  
  int size;
  int pf_ij = pf_index( i, j, seqlength);
  DBL_TYPE extraTerms;
  DBL_TYPE tempMin;
		
  int isEndNicked = FALSE;
	
  if( etaN[ EtaNIndex( i-0.5,i-0.5,seqlength)][0] == 1 || 
      etaN[ EtaNIndex( j+0.5,j+0.5,seqlength)][0] == 1) 
    isEndNicked = TRUE;
	
  if( L >= 12) {
    makeNewFx( i, j, seq, seqlength, etaN, Fb, Fx);
  }
	
  //Use extensible cases              
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    for( size = 8; size <= L - 4; size++) {
			
      extraTerms = InteriorMM( seq[i], seq[j], seq[i+1], 
			       seq[j-1]);
      
      tempMin = 
	Fx[ fbixIndex( j-i, i, size, seqlength)] +
	extraTerms;
      Fb[ pf_ij] = MIN( tempMin, Fb[ pf_ij]);
      minILoopEnergyBySize[ size] = MIN( tempMin, minILoopEnergyBySize[size]);
    }
  }
  
  if( L >= 12 && i != 0 && j != seqlength -1 && isEndNicked == FALSE) {
    extendOldFx( i, j, seqlength, Fx,Fx_2);
  }
  
  /* Add in inextensible cases */  
  if( CanPair( seq[ i], seq[j]) == TRUE) {
    //first check inextensible cases
    tempMin = MinInextensibleIL( i,j, seq, seqlength, Fb,  etaN, minILoopEnergyBySize);    
    Fb[ pf_ij] = MIN( Fb[ pf_ij], tempMin);
  } 
}

/* ******************************************* */
void makeNewFx( int i, int j, int seq[], int seqlength, 
		int **etaN, DBL_TYPE Fb[], DBL_TYPE Fx[]) {
	
  /*Determine the new entries of Fx(i,j,size) that are not extended 
    versions of Fx(i+1, j-1, size-2) */
  DBL_TYPE energy;
  int d, e; //Internal pair.(d, e will be restricted to special cases)
  
  int size, L1, L2; //size parameters: L1 + L2 = size, L1 = h-i-1, L2 = j-m-1
  DBL_TYPE tempMin;
  int fbix;
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
						
      tempMin = energy + Fb[ pf_index( d,e,seqlength)];
      fbix = fbixIndex( j-i,i,size,seqlength);
      Fx[ fbix ] = MIN( tempMin, Fx[ fbix]); 
			
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
                  
      tempMin = energy + Fb[ pf_index( d,e,seqlength)];
      fbix = fbixIndex( j-i,i,size,seqlength);
      Fx[ fbix ] = MIN( tempMin, Fx[ fbix]); 
			
    }
  }
	
}

  
/* *************** */
void extendOldFx( int i, int j, int seqlength, DBL_TYPE Fx[], DBL_TYPE Fx_2[]) {
  /* Extends all entries of Fx */
	
  int size;
  DBL_TYPE oldSizeEnergy;
  DBL_TYPE newSizeEnergy;
	
  for( size = 8; size <= (j - i + 1) - 4; size++) {
    if( size <= 30) {
      oldSizeEnergy = loop37[ size - 1];
    }
    else {
      oldSizeEnergy = loop37[ 30 - 1];
      oldSizeEnergy += sizeLog (size); //1.75*kB*TEMP_K*log( size/30.0); 
    }
    if( size + 2 <= 30) {
      newSizeEnergy = loop37[ size+2 - 1];
    }
    else {
      newSizeEnergy = loop37[ 30 - 1];
      newSizeEnergy += sizeLog (size+2); //1.75*kB*TEMP_K*log( (size+2)/30.0); 
    }
		
    Fx_2[ fbixIndex( j-i+2, i-1, size+2, seqlength)] = 
      Fx[ fbixIndex( j-i, i, size, seqlength)] + newSizeEnergy - oldSizeEnergy;
  }
	
}

/* ****************** */
DBL_TYPE MinInextensibleIL( int i, int j, int seq[], int seqlength, 
			    DBL_TYPE Fb[], int **etaN, DBL_TYPE *minILoopEnergyBySize) {
  /* This finds the minimum energy IL that has a special energy 
     calculation, i.e. small loops, bulge loops or GAIL case.  None of 
     these is allowed to be nicked
  */
	
  DBL_TYPE energy;
  int d, e; //Internal pair.(h, m will be restricted to special cases)  
  int L1, L2; //size parameters: L1 + L2 = size, L1 = h-i-1, L2 = j-m-1
  int size;
	
  DBL_TYPE tempMin;
  DBL_TYPE min_energy = NAD_INFINITY;
	
  /* Consider "small" loops with special energy functions */
  
  for( L1 = 0; L1 <= 3; L1++) {
    d = i + L1 + 1;
    for( L2 = 0; L2 <= MIN( 3, j-d-2); L2++) {
      e = j - L2 - 1;
      size = L1 + L2;
			
      if( CanPair( seq[d], seq[e]) == TRUE &&
	  (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
	  (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {
				
	energy = InteriorEnergy( i, j, d, e, seq);
				
	tempMin = energy + 
	  Fb[ pf_index( d, e, seqlength)];
	min_energy = MIN( min_energy, tempMin);
	minILoopEnergyBySize[ size] = MIN( tempMin, minILoopEnergyBySize[size]);
      }
    }
  }
	
  /* Next consider large bulges or large asymmetric loops */
  // Case 2a  L1 = 0,1,2,3, L2 >= 4;
  for( L1 = 0; L1 <= 3; L1++) {
    d = i + L1 + 1;
    for( L2 = 4; L2 <= j - d - 2; L2++) {
      e = j - L2 - 1;
      size = L1 + L2;

      if( CanPair( seq[d], seq[e]) == TRUE &&
	  (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
	  (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {
				
				
	energy = InteriorEnergy( i, j, d, e, seq);	
				
	tempMin = energy +
	  Fb[ pf_index( d, e, seqlength)]; 
	min_energy = MIN( tempMin, min_energy);
	minILoopEnergyBySize[ size] = MIN( tempMin, minILoopEnergyBySize[size]);
      }
    }
  }
	
  // Case 2b L1 >= 4, L2 = 0,1,2,3;
  for( L2 = 0; L2 <= 3; L2++) {
    e = j - L2 - 1;
    for( L1 = 4;  L1 <= e - i - 2; L1++) {
      d = i + L1 + 1;
      size = L1 + L2;
			
      if( CanPair( seq[d], seq[e]) == TRUE &&
	  (etaN[ EtaNIndex(i+0.5, d-0.5,seqlength)][0] == 0) &&
	  (etaN[ EtaNIndex(e+0.5, j-0.5,seqlength)][0] == 0) ) {
				
	energy = InteriorEnergy( i, j, d, e, seq);
	tempMin = energy + 
	  Fb[ pf_index( d, e, seqlength)];
	min_energy = MIN( min_energy, tempMin);
	minILoopEnergyBySize[ size] = MIN( tempMin, minILoopEnergyBySize[size]);
      }
    }
  }    
	
  return min_energy;
}

/* ******************************************************* */
/* Fs, Fms  Recursion */

void MakeFs_Fms( int i, int j, int seq[], int seqlength, 
                  DBL_TYPE *Fs, DBL_TYPE *Fms, DBL_TYPE *Fb,
                  int *nicks, int **etaN) {
  
  int d; //base pair is i,d
  DBL_TYPE bp_penalty = 0.0;
  int pf_ij = pf_index( i, j, seqlength);
  
  DBL_TYPE extraTerms;
  int nNicks;
  int index_ij = EtaNIndex( i+0.5, j-0.5, seqlength);
  int start;
  DBL_TYPE tempMin;
  
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
         
         Fs[ pf_ij] = MIN( tempMin, Fs[ pf_ij]);
         
         // ******************** 
         
         extraTerms =  DangleEnergy( d+1, j, seq, seqlength)+
           bp_penalty + ALPHA_2 + ALPHA_3*(j-d);
         
         tempMin = Fb[ pf_index( i, d, seqlength) ] + 
           extraTerms;
         
         Fms[ pf_ij] = MIN( tempMin, Fms[ pf_ij]);
       }
  }
}


/* ******************************* */
/* F, Fm Recursions */
void MakeF_Fm_N3( int i, int j, int seq[], int seqlength, 
                 DBL_TYPE *F, DBL_TYPE *Fs, 
                 DBL_TYPE *Fms, DBL_TYPE *Fm,
                 int *nicks, int **etaN) {
 
 int d;//left base of rightmost base pair.
   int pf_ij = pf_index( i, j, seqlength);
 
 DBL_TYPE extraTerms;
 DBL_TYPE tempMin;
 
 F[ pf_ij] = NickedEmptyF( i, j, nicks, seq, seqlength, etaN);
 
 for( d = i; d <= j - 1; d++) {
   if( etaN[ EtaNIndex(d-0.5, d-0.5, seqlength)][0] == 0 || d == i ) {
     
     tempMin = F[ pf_index(i, d-1, seqlength)] + Fs[ pf_index( d, j, seqlength)];
     F[ pf_ij] = MIN( tempMin, F[ pf_ij]);
     
     extraTerms = DangleEnergy( i, d-1, seq, seqlength) +
       (ALPHA_3)*(d-i);
     
     
     if( etaN[ EtaNIndex( d-0.5, d-0.5, seqlength)][0] == 0) { 
       //otherwise Qm not possible
       
       if( etaN[ EtaNIndex(i+0.5, d-0.5, seqlength)][0] == 0 ) {
         tempMin = Fms[ pf_index( d, j, seqlength)] +
           extraTerms;
         Fm[ pf_ij] = MIN( tempMin, Fm[ pf_ij]);
       }
       
       if( d >= i+2) {
         tempMin = Fm[ pf_index( i, d - 1, seqlength) ] +
           Fms[ pf_index( d, j, seqlength) ];
         
         Fm[ pf_ij] = MIN( tempMin, Fm[ pf_ij]);
       }
     }
   }
 }
}
/* *************** */


DBL_TYPE MinInterior_Multi( int i, int j, int seq[], int seqlength, 
			    DBL_TYPE *Fm, DBL_TYPE *Fb, int *nicks, 
			    int **etaN ){
  // This finds all possible internal loops (no pseudoknots)
  // closed on the "outside" by bases i and j, as well as all 
  // multiloops
	
	
  DBL_TYPE min_energy = NAD_INFINITY;
  DBL_TYPE tempMin;
  int d, e; // d - e is internal basepair 
  DBL_TYPE bp_penalty = 0;
	
	
  for( d = i+1; d <= j - 5; d++) {
    for( e = d + 4; e <= j - 1; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE) {
        bp_penalty = 0.0;
				
	if( etaN[ EtaNIndex(e+0.5, j-0.5, seqlength)][0] == 0) {
				
	  if( etaN[ EtaNIndex(i+0.5, d-0.5, seqlength)][0] == 0) {
	    tempMin = InteriorEnergy( i, j, d, e, seq) +
	      Fb[ pf_index( d, e, seqlength) ];
					
	    min_energy = MIN( tempMin, min_energy);
	  }
					
	  if( etaN[ EtaNIndex(d-0.5, d-0.5, seqlength)][0] == 0 &&
	      etaN[ EtaNIndex(i+0.5, i+0.5, seqlength)][0] == 0 &&
	      d>= i+6 && CanWCPair(seq[d], seq[e]) &&
	      CanWCPair(seq[i], seq[j]) ) {
						
	    if( seq[d] != BASE_C && seq[e] != BASE_C) {
	      bp_penalty = AT_PENALTY;
	    }
	    if( seq[i] != BASE_C && seq[j] != BASE_C) {
	      bp_penalty += AT_PENALTY;
	    }
	    tempMin = Fm[ pf_index(i+1, d-1, seqlength)] +
	      Fb[ pf_index( d, e, seqlength)] +
	      (ALPHA_1 + 2*ALPHA_2 + ALPHA_3*(j-e-1) + bp_penalty) +
	      DangleEnergy( e+1, j-1, seq, seqlength);
	    min_energy = MIN( tempMin, min_energy);
	  }
        }
				
      }
			
    }
  }
  return min_energy;
}

/* **************** */
#ifdef O_N4

void MakeF_Fm_N4( int i, int j, int seq[], int seqlength, 
		  DBL_TYPE F[], DBL_TYPE Fm[], DBL_TYPE Fb[] ){

  int d, e; // d - e is internal basepair 
  DBL_TYPE bp_penalty = 0;
  int pf_ij = pf_index(i, j, seqlength);
  DBL_TYPE tempMin;

  F[ pf_ij] = DangleEnergy(i, j, seq, seqlength);  //Empty Graph

  for( d = i; d <= j - 4; d++) {
    for( e = d + 4; e <= j; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE &&
	        CanWCPair(seq[d], seq[e]) ) {
	bp_penalty = 0;
	if( seq[d] != BASE_C && seq[e] != BASE_C) {
	  bp_penalty = AT_PENALTY;
	}
				
	tempMin = F[ pf_index(i, d-1, seqlength)] +
	  Fb[ pf_index( d, e, seqlength) ] +
	  bp_penalty + 
	  DangleEnergy( e+1, j, seq, seqlength); 
	F[ pf_ij] = MIN( tempMin, F[ pf_ij]);
				
	tempMin =
	  (ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty) +
	  Fb[ pf_index( d, e, seqlength)] +
	  DangleEnergy( e+1, j, seq, seqlength) +
	  DangleEnergy( i, d-1, seq, seqlength);
	Fm[ pf_ij] = MIN( tempMin, Fm[ pf_ij]);
				
				
	if( d >= i+5) {
	  tempMin = Fm[ pf_index(i, d-1, seqlength)] +
	    Fb[ pf_index( d, e, seqlength)] +
	    (ALPHA_2 + ALPHA_3*(j-e) + bp_penalty) +
	    DangleEnergy( e+1, j, seq, seqlength);
	  Fm[ pf_ij] = MIN( tempMin, Fm[ pf_ij]);
					
	}
				
      }
    
    }
  }


}
#endif

/* *********************** */







