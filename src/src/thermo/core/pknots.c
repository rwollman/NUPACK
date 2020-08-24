/*
  pknots.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 8/2001, Justin Bois 1/2007

  This contains subroutines for the folding and partition function
  algorithms found in "fold.c" and "pfunc.c".  The functions in this
  file will help generate and assign energies to pseudoknots involving
  two specified bases as the pknot's 5' and 3' ends.
*/

#include "pknots.h"

/* ********************************************* */
DBL_TYPE MinFb_Pk( int i, int j, int seq[], int seqlength, 
		   DBL_TYPE Fp[], DBL_TYPE Fm[] ) {
  // Determine all possible, rightmost pseudoknots 
  // contained in i,j "loop" and calculate Sumexpl
  
	
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE min_energy = NAD_INFINITY;
  DBL_TYPE tempMin;
  int d, e; //the left and right ends of a pseudoknot
  int pf_de;
	
  if( ( seq[i]) + ( seq[j]) != 5) {
    return NAD_INFINITY;
  }
	
  if( seq[i] != BASE_C  && seq[j] != BASE_C) {
    bp_penalty = AT_PENALTY;
  }
  
  for( d = i+1; d <= j - 6; d++) {
    for( e = d + 5; e <= j - 1; e++) {
      pf_de = pf_index( d, e, seqlength);
			
      tempMin = 
	(BETA_1M + ALPHA_1 + 3*ALPHA_2 + 
	 (j-e-1 + d-i-1)*ALPHA_3 + bp_penalty) +
	DangleEnergy( i+1, d-1, seq, seqlength) +
	DangleEnergy( e+1, j-1, seq, seqlength) +
	Fp[ pf_de];
			
      min_energy = MIN( min_energy, tempMin);
			
      tempMin = Fm[ pf_index( i+1, d-1, seqlength)] + 
	Fp[ pf_de] + ( BETA_1M + ALPHA_1 + 3*ALPHA_2 + bp_penalty + 
		       (j - e - 1)*ALPHA_3 ) +
	DangleEnergy( e+1, j-1, seq, seqlength);
			
      min_energy = MIN( min_energy, tempMin);
			
    }
  }
	
  return min_energy;
}


/* ********************************** */
void MakeFg_N5( int i, int j, int seq[], int seqlength, DBL_TYPE *Fg,
		DBL_TYPE *Fm, DBL_TYPE *Fgls, DBL_TYPE *Fgrs, DBL_TYPE *FgIx,
		DBL_TYPE *FgIx_2, short *possiblePairs) {
  //  Make the gap matrix for Qg
  int c,d,e,f;
  int gap_idej;
  DBL_TYPE IJ_bp_penalty = 0.0;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE tempMin;
	
  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }
	
  if( CanPair( seq[i], seq[j]) == TRUE) {
		
    //Case 0: Only 1 pair
    Fg[ gap_index( i,i,j,j,seqlength)] = 0;
    
    //Case 1:  Terminal Inner Pair
    for( d = i+1; d <= j-5; d++) {
      for( e = d+4; e <= j-1; e++) {      
	if( CanPair( seq[d], seq[e] ) == TRUE) {
					
	  gap_idej = gap_index( i, d, e, j, seqlength);
					
	  tempMin = 
	    InteriorEnergy( i, j, d, e, seq);
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
					
	}
      }
    }
  }
	
  fastIloop_Fg(i, j, seq, seqlength, Fg, FgIx, FgIx_2, possiblePairs);
	
	
  if( CanPair( seq[i], seq[j]) == TRUE && 
      CanWCPair(seq[i], seq[j])) {
		
    //Case 2:  Multiloop Left
    for( d = i+6; d <= j-5; d++) {
      for( e = d+4; e <= j-1; e++) {      
	if( CanPair( seq[d], seq[e] ) == TRUE 
	    && CanWCPair(seq[d], seq[e])) {
	  gap_idej = gap_index(i,d,e,j,seqlength);
					
	  bp_penalty = IJ_bp_penalty;
	  if( seq[d] != BASE_C && seq[e] != BASE_C) {
	    bp_penalty += AT_PENALTY;
	  }
					
	  tempMin = 
	    Fm[ pf_index( i+1, d-1, seqlength)] +
	    (ALPHA_1 + 2*ALPHA_2 + (j-e-1)*ALPHA_3 + bp_penalty) +
	    DangleEnergy( e+1, j-1, seq, seqlength);
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
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
					
	  tempMin = 
	    Fm[ pf_index( e+1, j-1, seqlength)] +
	    (ALPHA_1 + 2*ALPHA_2 + (d-i-1)*ALPHA_3 + bp_penalty) +
	    DangleEnergy( i+1, d-1, seq, seqlength) ;
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
	}
      }
    }
		
    //Case 4: Multiloop Both Sides
    for( d = i+6; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {      
	if( CanPair( seq[d], seq[e] ) == TRUE
	    && CanWCPair(seq[d], seq[e])
	    ) {
	  gap_idej = gap_index(i,d,e,j,seqlength);
					
	  bp_penalty = IJ_bp_penalty;
	  if( seq[d] != BASE_C && seq[e] != BASE_C) {
	    bp_penalty += AT_PENALTY;
	  }	
					
	  tempMin = 
	    Fm[ pf_index( i+1, d-1, seqlength)] +
	    Fm[ pf_index( e+1, j-1, seqlength)] +
	    (ALPHA_1 + 2*ALPHA_2 + bp_penalty);
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
					
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
						
	    tempMin = 
	      (ALPHA_1+ALPHA_2 + (j-f-1)*ALPHA_3+ bp_penalty) + 
	      Fgls[ gap_index( i+1, d, e, f, seqlength)] +
	      DangleEnergy( f+1, j-1, seq, seqlength);
	    Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
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
						
	    tempMin = 
	      (ALPHA_1+ALPHA_2 + (c-i-1)*ALPHA_3 + 
	       bp_penalty) + 
	      Fgrs[ gap_index( c, d, e, j-1, seqlength)] +
	      DangleEnergy( i+1, c-1, seq, seqlength);
	    Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
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
						
	    tempMin = 
	      Fm[ pf_index( i+1, c-1, seqlength)] +
	      (ALPHA_1+ALPHA_2+bp_penalty) +
	      Fgrs[ gap_index( c, d, e, j-1, seqlength)];
	    Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
						
	  }	
	}
      }
    }
  } // Can Pair
}

/* ******************************************* */
void fastIloop_Fg(int i, int j, int seq[], int seqlength, 
		  DBL_TYPE *Fg, DBL_TYPE *FgIx, DBL_TYPE *FgIx_2,
		  short *possiblePairs) {
  
  int d, e, size;
  int L = j - i + 1;
  int pf_ij = pf_index(i, j, seqlength);
  int index;
  int gap_idej;
  DBL_TYPE tempMin;
	
  //FgIx recurions
  //expllicitly add in base cases to QgIx (smallest extensible loops) 
  
  if( L >= 17 && possiblePairs[pf_ij] == TRUE) {
    makeNewFgIx( i, j, seq, seqlength, Fg, FgIx);
  }
  
  for( d = i+1; d <= j-5; d++) {
    for( e = d+4; e <= j-1; e++) {
      if(  CanPair( seq[d], seq[e]) == TRUE) {
				
	gap_idej = gap_index( i,d,e,j, seqlength);
				
	if( L>=17 &&  CanPair( seq[i], seq[j]) == TRUE) { 
					
	  index = QgIxIndex( j-i, i, 8, d, e, seqlength);
					
	  for( size = 8; size <= L - 9; size++) {
	    tempMin = FgIx[ index] + 
	      InteriorMM( seq[i], seq[j], seq[i+1], seq[j-1]);
	    Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]); 
	    index = index + (L-2)*(L-3)/2;
	  }
					
	}
				
	if( i != 0 && j != seqlength - 1) {
	  extendOldFgIx( i, j, d, e, seq, seqlength, 
			 Fg, FgIx, FgIx_2);
	}
				
				
	tempMin = MinInextensibleIL_Fg( i, j, d, e, seq, seqlength,
					Fg);
	Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
      }
    }
  }
}

/* ************************************************* */
DBL_TYPE MinInextensibleIL_Fg( int i, int j, int d, int e, 
			       int seq[], int seqlength, 
			       DBL_TYPE *Fg) {
  /* This finds the minimum energy IL that has a special energy 
     calculation, i.e. small loops, bulge loops or GAIL case */
	
  DBL_TYPE energy;
  int c, f; //Internal pair.(c, f will be restricted to special cases)
  
  int L1, L2; //size parameters: L1 + L2 = size, L1 = c-i-1, L2 = j-f-1
	
  DBL_TYPE tempMin;
  DBL_TYPE min_energy = NAD_INFINITY;
	
  if( CanPair( seq[i], seq[j]) == FALSE) { 
    return NAD_INFINITY;
  }
	
  /* First consider small loops */
  // L1 <= 3, L2 <= 3
  for( L1 = 0; L1 <= MIN(3, d-i-2); L1++) {
    c = i + L1 + 1;
    for( L2 = 0; L2 <= MIN( 3, j-e-2); L2++) {      
      f = j - L2 - 1;
      if( CanPair( seq[c], seq[f]) == TRUE) { //c < d && f > e
				
	energy = InteriorEnergy( i, j, c, f, seq);
				
	tempMin = energy +
	  Fg[ gap_index( c, d, e, f, seqlength)];
	min_energy = MIN( tempMin, min_energy);
				
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
	energy = InteriorEnergy( i, j, c, f, seq);	
	tempMin = energy +
	  Fg[ gap_index( c, d, e, f, seqlength)];
	min_energy = MIN( tempMin, min_energy);
      }
    }
  }
	
  // Case 2b L1 >= 4, L2 = 0,1,2,3;
  for( L2 = 0; L2 <= MIN(3,j-e-2) ; L2++) {
    f = j - L2 - 1;
    for( L1 = 4; L1 <= d-i-2; L1++) {
      c = i + L1 + 1;
      if( CanPair( seq[c], seq[f]) == TRUE) { //c < d && f > e
				
	energy = InteriorEnergy( i, j, c, f, seq);
	tempMin = energy + 
	  Fg[ gap_index( c, d, e, f, seqlength)];
	min_energy = MIN( tempMin, min_energy);
				
      }
    }
  }    
	
  return min_energy;
	
}  
/* *********************************************** */


void makeNewFgIx( int i, int j, int seq[], int seqlength, 
		  DBL_TYPE *Fg, DBL_TYPE *FgIx) {
	
  /*Determine the new entries of FgIx(i,j,d,e,size) that are not extended 
    versions of FgIx(i+1, j-1, h1, e, size-2) */
  
  DBL_TYPE energy;
  int c, f; //Internal pair (c, f will be restricted to special cases);
  int d,e;  //Internal pair of gap matrix;
		
  int size, L1, L2; //size parameters: L1 + L2 = size, L1 = c-i-1, L2 = j-f-1
  int partial_index; //used to limit calls to QgIxIndex
  int len1_len2 = (j-i-1)*(j-i-2)/2;
  DBL_TYPE tempMin;
		
  //Add in all the cases that are not an extended version of an
  //extensible case.
		
		
  for( d = i + 6; d <= j - 10; d++) {
    for( e = d + 4; e <= j-6; e++) {
      if( CanPair( seq[d], seq[e]) ) {
	partial_index = QgIxIndex( j-i, i, 0, d, e, seqlength);
				
	//Case 1:  L1 = 4, L2 >= 4;
	L1 = 4;
	c= i + L1 + 1;
	for( L2 = 4; L2 <= j - e - 2; L2++) {
	  f = j - L2 - 1;
	  size = L1 + L2;
					
	  if( CanPair( seq[c], seq[f]) == TRUE) {
						
	    energy = asymmetryEfn( L1, L2, size);							
	    energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
	    /*Exclude the i-j stacking energy here, just in case i-j 
	      don't pair */
						
	    tempMin = 
	      energy + Fg[ gap_index( c, d, e, f, seqlength)];
						
	    FgIx[ partial_index+size*len1_len2 ] = 
	      MIN( FgIx[ partial_index+size*len1_len2 ], tempMin); 
						
	  }
	}    
				
	//Case 2  L1 > 4, L2 = 4
	if( d >= i+7) {
	  L2 = 4;
	  f = j-L2-1;
	  for( L1 = 5; L1 <= d-i-2; L1++) {
	    c = i + L1 + 1;
	    size = L1 + L2;
						
	    if( CanPair( seq[c], seq[f]) == TRUE) {
							
	      energy = asymmetryEfn( L1, L2, size);
	      energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
	      /*Exclude the i-j stacking energy here, just in case i-j 
		don't pair */
							
							
	      tempMin =
		energy + Fg[ gap_index( c, d,e,f, seqlength)];
							
	      FgIx[ partial_index + size*len1_len2] =
		MIN( tempMin, FgIx[ partial_index + size*len1_len2]); 
							
	    }
	  }
	}
      } 
    }
  }
		
}


/* ************ for FASTILOOP ************** */

void extendOldFgIx( int i, int j, int d, int e, int seq[], int seqlength, 
		    DBL_TYPE *Fg, DBL_TYPE *FgIx, DBL_TYPE *FgIx_2) {
  /* Extends all entries of FgIx */
	
	
  int L = j - i + 1;
  int size;
  //extern DBL_TYPE *sizeTerm; //precomputed
  
  //These following variables are to minimze calls to QgIxIndex;
  int partial_index, partial_index2;  
  int len_12 = (L-2)*(L-3)/2;
  int len10 = L*(L-1)/2;
	
  partial_index2 = QgIxIndex( j-i+2, i-1, 0, d, e, seqlength);
  partial_index = QgIxIndex( j-i, i, 0, d, e, seqlength);
  
  for( size = 8; size <= L - 9; size++) {
    FgIx_2[ partial_index2 + (size+2)*len10] = 
      FgIx[ partial_index + size*len_12] +
      sizeTerm[ size];    
  }
}


/* ************************************ */

void MakeFgls( int i, int j, int seq[], int seqlength, DBL_TYPE *Fg, 
	       DBL_TYPE *Fm,  DBL_TYPE *Fgls) {
  int c, d, e;
  int pf_ic1;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE tempMin;
  int gap_idej;
	
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
						
	    // The d-e bp_penalty already counted in MakeQgl
	    gap_idej = gap_index( i, d, e, j, seqlength);
						
	    tempMin =
	      Fm[ pf_ic1] +
	      Fg[ gap_index(c, d, e, j, seqlength)] +
	      (ALPHA_2 + bp_penalty);
						
	    Fgls[ gap_idej] = 
	      MIN( tempMin, Fgls[ gap_idej]);
						
	  }
	}
      }
    }
  }
}
/* ************************** */
void MakeFgrs( int i, int j, int seq[], int seqlength, DBL_TYPE *Fg, 
	       DBL_TYPE *Fm, DBL_TYPE *Fgrs) {
	
  int d, e, f;
  int gap_idej;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE tempMin;
	
  for( d = i + 1; d <= j-10; d++) {
    for( e = d+4; e <= j-6; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE ) {
	gap_idej = gap_index(i,d,e,j,seqlength);
	for( f = e+1; f <= j-5; f++) {
	  if( CanPair( seq[i], seq[f]) == TRUE &&
	      CanWCPair(seq[i], seq[f]) 
	      ) {
						
	    bp_penalty = 0;
	    if( seq[i] != BASE_C && seq[f] != BASE_C) {
	      bp_penalty = AT_PENALTY;
	    }
						
	    tempMin = 
	      Fm[ pf_index( f+1, j, seqlength)] + 
	      Fg[ gap_index( i, d, e, f, seqlength)] +
	      (ALPHA_2 + bp_penalty);
	    Fgrs[ gap_idej] = MIN( tempMin, Fgrs[ gap_idej]);
						
	  }
	}
      }
    }
  }
}
/* *************************************** */

void MakeFgl( int i, int j, int seq[], int seqlength, 
	      DBL_TYPE *Fg, DBL_TYPE *Fgl, DBL_TYPE *Fz) {
  //make sure to call this AFTER MakeFg and BEFORE MakeFgr
  
  int d, e, f;
  int gap_idfj, gap_iefj;
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE tempMin;
  
  for( d = i+1; d <= j-5; d++) {
    for( f = d+4; f <= j-1; f++) {
      if( CanPair( seq[d], seq[f]) == TRUE &&
	  CanWCPair(seq[d], seq[f])
	  ) {
	bp_penalty = 0.0;
	if( seq[d] != BASE_C && seq[f] != BASE_C) {
	  bp_penalty = AT_PENALTY;
	}
				
	gap_idfj =  gap_index(i, d, f, j, seqlength);
	for( e = d; e <= f-2; e++) {
					
	  gap_iefj = gap_index( i, e, f, j, seqlength);
	  tempMin = 
	    Fg[ gap_idfj]+Fz[ pf_index( d+1, e, seqlength)] +
	    (BETA_2 + bp_penalty);
	  Fgl[ gap_iefj ] =
	    MIN( tempMin, Fgl[ gap_iefj]);
	}
      }
    }
  }
}

/* ******************************** */

void MakeFgr( int i, int j, int seq[], int seqlength, 
	      DBL_TYPE *Fgr, DBL_TYPE *Fgl, DBL_TYPE *Fz) {
  
  //make sure to call this AFTER MakeQg and AFTER MakeQgl
  
  int d, e, f;
  int gap_idej;
  DBL_TYPE tempMin;
	
  for( d = i+1; d <= j-3; d++) {
    for( e = d+2; e <= j-1; e++) {
      gap_idej =  gap_index(i, d, e, j, seqlength);
      for( f = e; f <= j-1; f++) {
	tempMin =
	  Fgl[ gap_index(i,d,f,j,seqlength)]+
	  Fz[ pf_index( e, f-1, seqlength)];
	Fgr[ gap_idej] = MIN( tempMin, Fgr[ gap_idej]);
				
      }
    }
  }  
}

/* ********************************** */


DBL_TYPE MinFp_N5( int i, int j, int seq[], int seqlength, 
		   DBL_TYPE *Fgl, DBL_TYPE *Fgr, DBL_TYPE *Fg,
		   DBL_TYPE *Fz) {
  
  int d, e, f;
  DBL_TYPE bp_penalty, new_bp_penalty;
  DBL_TYPE min_energy = NAD_INFINITY;
  DBL_TYPE tempMin;
	
  int a, b, c;
  if( j - i <= 4) {
    return min_energy;
  }
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
					
	  tempMin = 
	    Fg[ gap_index( i, a, d, e, seqlength) ]+
	    Fg[ gap_index( b, c, f, j, seqlength) ]+
	    (bp_penalty + 2*BETA_2) +
	    Fz[ pf_index( e + 1, f - 1, seqlength)]+
	    Fz[ pf_index( c + 1, d - 1, seqlength)]+
	    Fz[ pf_index( a + 1, b - 1, seqlength)];
	  min_energy = MIN( tempMin, min_energy);
	}
      }
    }
  }
	
  
  if( j - i <= 6) {
    return min_energy;
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
					
	  tempMin = Fg[ gap_index(i, i, e, f, seqlength)] +
	    Fz[ pf_index( i+1, d-1, seqlength)] +
	    Fgr[ gap_index(d, e-1, f+1, j, seqlength)] +
	    new_bp_penalty + BETA_2;
	  min_energy = MIN( tempMin, min_energy);
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
						
	    tempMin = Fgl[ gap_index(i, d-1, e, f, seqlength)] +
	      Fg[ gap_index(d, d, j, j, seqlength)] +
	      Fz[ pf_index( d+1, e-1, seqlength)] +
	      Fz[ pf_index( f+1,j-1, seqlength)] +
	      new_bp_penalty+BETA_2;
	    min_energy = MIN( tempMin, min_energy);
						
	  }
	}
      }
    }
  }
  
	
  if( j - i <= 7) {// not enough room for pk
    return NAD_INFINITY;
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
						
	    tempMin = Fgl[ gap_index(i, d-1, e, f, seqlength)] +
	      Fgr[ gap_index(d, e-1, f+1, j, seqlength)] +
	      (new_bp_penalty);
	    min_energy = MIN( tempMin, min_energy);
						
	  }
	}
      }
    }
  }
	
  return min_energy;
}

/* ***************************** */

void MakeF_Fm_Fz( int i, int j, int seq[], int seqlength,
		  DBL_TYPE *F, DBL_TYPE *Fm, DBL_TYPE *Fz,
		  DBL_TYPE *Fb, DBL_TYPE *Fp) {
  
  int d, e; // d - e is internal basepair or pk boundary
  DBL_TYPE bp_penalty;
  int pf_ij = pf_index(i, j, seqlength);
  int pf_de;
  DBL_TYPE tempMin;
	
  F[ pf_ij] = DangleEnergy(i, j, seq, seqlength);  //Empty Graph
	
  
  if( i != 0 && j != seqlength - 1) { 
    Fz[ pf_ij] = DangleEnergy(i, j, seq, seqlength) +
      BETA_3*(j-i+1);
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
				
	tempMin= F[ pf_index(i, d-1, seqlength)] +
	  Fb[ pf_de ] +
	  bp_penalty +
	  DangleEnergy( e+1, j, seq, seqlength); 
	F[ pf_ij] = MIN( F[ pf_ij], tempMin);
				
	if( i != 0 && j != seqlength - 1) {
					
					
	  tempMin =
	    (ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty) +
	    Fb[ pf_de] +
	    DangleEnergy( e+1, j, seq, seqlength) +
	    DangleEnergy( i, d-1, seq, seqlength);
	  Fm[ pf_ij] = MIN( Fm[ pf_ij], tempMin);
					
					
	  if( d >= i+5) {
	    tempMin = Fm[ pf_index(i, d-1, seqlength)] +
	      Fb[ pf_de] +
	      (ALPHA_2 + ALPHA_3*(j-e) + bp_penalty) +
	      DangleEnergy( e+1, j, seq, seqlength);
	    Fm[ pf_ij] = MIN( Fm[ pf_ij], tempMin);
	  }
					
					
	  tempMin = Fz[ pf_index(i, d-1, seqlength)] + 
	    Fb[ pf_de] + (BETA_2 + BETA_3*(j-e) + bp_penalty) +
	    DangleEnergy( e+1, j, seq, seqlength);
	  Fz[ pf_ij] = MIN( Fz[ pf_ij], tempMin);
	}
				
      }
      
    }
  }
	
  for( d = i; d <= j - 5; d++) {
    for( e = d + 5; e <= j; e++) {
      
      pf_de = pf_index( d, e, seqlength);
      
      tempMin = F[ pf_index(i, d-1, seqlength)] +
	Fp[ pf_de ] +
	BETA_1 +
	DangleEnergy( e+1, j, seq, seqlength); 
      F[ pf_ij] = MIN( tempMin, F[pf_ij]);
			
      if( i != 0 && j != seqlength - 1) {
				
	tempMin =
	  (BETA_1M + 2*ALPHA_2 + ALPHA_3*(d-i + j-e)) +
	  Fp[ pf_de] +
	  DangleEnergy( e+1, j, seq, seqlength) +
	  DangleEnergy( i, d-1, seq, seqlength);
	Fm[ pf_ij] = MIN( Fm[ pf_ij], tempMin);
				
	if( d >= i+5) {
	  tempMin = Fm[ pf_index(i, d-1, seqlength)] +
	    Fp[ pf_de] + (BETA_1M + 2*ALPHA_2 + ALPHA_3*(j-e)) +
	    DangleEnergy( e+1, j, seq, seqlength);
	  Fm[ pf_ij] = MIN( tempMin, Fm[ pf_ij]);
	}
				
	tempMin = Fz[ pf_index(i, d-1, seqlength)] +
	  Fp[ pf_de] + (BETA_1P + 2*BETA_2 + BETA_3*(j-e)) +
	  DangleEnergy( e+1, j, seq, seqlength) ;	
	Fz[ pf_ij] = MIN( tempMin, Fz[ pf_ij]);
      }
    }    
  }
}



/* ******* */

#ifdef O_N8
void MakeFg_N8( int i, int j, int seq[], int seqlength, DBL_TYPE Fg[],
		DBL_TYPE Fm[]) {
  //  Make the gap matrix for Fg
  int c,d,e,f;
  int gap_idej;
  DBL_TYPE IJ_bp_penalty = 0.0;
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE tempMin;
	
  if( CanPair( seq[i], seq[j]) == FALSE) {
    return;
  }
	
  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }
	
  //Case 0: Only 1 pair
  Fg[ gap_index( i,i,j,j,seqlength)] = 0;
	
  //Case 1:  Terminal Inner Pair
  for( d = i+1; d <= j-5; d++) {
    for( e = d+4; e <= j-1; e++) {      
      if( CanPair( seq[d], seq[e] ) == TRUE) {	
				
	gap_idej = gap_index(i,d,e,j, seqlength);
				
	tempMin = InteriorEnergy( i, j, d, e, seq); 
	Fg[ gap_idej] =
	  MIN( tempMin, Fg[ gap_idej]);
				
      }
    }
  }
  
  
  //Case 2:  Simple Interior Loop
  for( d = i+2; d <= j-6; d++) {
    for( e = d+4; e <= j-2; e++) {      
      gap_idej = gap_index(i,d,e,j,seqlength);
      if( CanPair( seq[d], seq[e] ) == TRUE) { 
				
	for( c = i+1; c <= d-1; c++) {
	  for( f = e+1; f <= j-1; f++) {
	    if( CanPair( seq[c], seq[f]) == TRUE) {
							
	      tempMin = 
		InteriorEnergy( i, j, c, f, seq) + 
		Fg[ gap_index( c, d, e, f, seqlength)];
	      Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
							
							
	    }
	  }
	}
      }
    }
  }
	
  if( CanWCPair(seq[i], seq[j])) {
		
    //Case 3:  Multiloop Left
    for( d = i+6; d <= j-5; d++) {
      for( e = d+4; e <= j-1; e++) {      
	if( CanPair( seq[d], seq[e] ) == TRUE &&
	    CanWCPair(seq[d], seq[e])) {
	  gap_idej = gap_index(i,d,e,j,seqlength);
					
	  bp_penalty = IJ_bp_penalty;
	  if( seq[d] != BASE_C && seq[e] != BASE_C) {
	    bp_penalty += AT_PENALTY;
	  }
					
	  tempMin = 
	    Fm[ pf_index( i+1, d-1, seqlength)] +
	    (ALPHA_1 + 2*ALPHA_2 + (j-e-1)*ALPHA_3 + bp_penalty) +
	    DangleEnergy( e+1, j-1, seq, seqlength);
					
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
					
	}
      }
    }
    
    //Case 4:  Multiloop Right
    for( d = i+1; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {      
	if( CanPair( seq[d], seq[e] ) == TRUE &&
	    CanWCPair(seq[d], seq[e])) {
	  gap_idej = gap_index(i,d,e,j,seqlength);
					
	  bp_penalty = IJ_bp_penalty;
	  if( seq[d] != BASE_C && seq[e] != BASE_C) {
	    bp_penalty += AT_PENALTY;
	  }
					
	  tempMin = 
	    Fm[ pf_index( e+1, j-1, seqlength)] +
	    (ALPHA_1 + 2*ALPHA_2 + (d-i-1)*ALPHA_3 + bp_penalty) +
	    DangleEnergy( i+1, d-1, seq, seqlength);
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
					
	}
      }
    }
    
    
    //Case 5: Multiloop Both Sides
    for( d = i+6; d <= j-10; d++) {
      for( e = d+4; e <= j-6; e++) {      
	if( CanPair( seq[d], seq[e] ) == TRUE && 
	    CanWCPair(seq[d], seq[e])) {
	  gap_idej = gap_index(i,d,e,j,seqlength);
					
	  bp_penalty = IJ_bp_penalty;
	  if( seq[d] != BASE_C && seq[e] != BASE_C) {
	    bp_penalty += AT_PENALTY;
	  }
					
	  tempMin = 
	    Fm[ pf_index( i+1, d-1, seqlength)] +
	    Fm[ pf_index( e+1, j-1, seqlength)] +
	    (ALPHA_1 + 2*ALPHA_2 + bp_penalty);
					
	  Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
	}
      }
    }
    
    
    //Case 6: Multiloop Left + More Qg
    for( d = i+7; d <= j-6; d++) {
      for( e = d+4; e <= j-2; e++) {      
	gap_idej = gap_index(i,d,e,j,seqlength);
	if( CanPair( seq[d], seq[e] ) == TRUE ) { 
					
	  for( c = i+6; c <= d-1; c++) {
	    for( f = e+1; f <= j-1; f++) {
	      if( CanPair( seq[c], seq[f]) == TRUE &&
		  CanWCPair(seq[c], seq[f])) {
								
		bp_penalty = IJ_bp_penalty;
		if( seq[c] != BASE_C && seq[f] != BASE_C) {
		  bp_penalty += AT_PENALTY;
		}
								
		tempMin = 
		  Fm[ pf_index( i+1, c-1, seqlength)] +
		  (ALPHA_1+2*ALPHA_2+(j-f-1)*ALPHA_3 + bp_penalty) + 
		  Fg[ gap_index( c, d, e, f, seqlength)] +
		  DangleEnergy( f+1, j-1, seq, seqlength);
		Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
								
	      }
	    }
	  }
	}
      }
    }
		
    
    //Case 7: Multiloop Right + More Qg
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
								
		tempMin = 
		  Fm[ pf_index( f+1, j-1, seqlength)] +
		  (ALPHA_1+2*ALPHA_2+(c-i-1)*ALPHA_3+bp_penalty) + 
		  Fg[ gap_index( c, d, e, f, seqlength)] +
		  DangleEnergy( i+1, c-1, seq, seqlength);
								
								
		Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
								
	      }
	    }
	  }
	}
      }
    }
    
    
    //Case 8: Multiloop Both sides + More Qg
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
								
		tempMin = 
		  Fm[ pf_index( i+1, c-1, seqlength)] +
		  Fm[ pf_index( f+1, j-1, seqlength)] +
		  (ALPHA_1+2*ALPHA_2+bp_penalty) +
		  Fg[ gap_index( c, d, e, f, seqlength)];
								
		Fg[ gap_idej] = MIN( tempMin, Fg[ gap_idej]);
	      }
	    }
	  }
	}
      }
    }
  }
}


DBL_TYPE MinFp_N8( int i, int j, int seq[], int seqlength, 
		   DBL_TYPE Fg[], DBL_TYPE Fz[]) {
	
  // Add all pseudoknots with ends i, j.
  // No penalty for non-AT, terminal basepairs need to be added for the
  // i and j pairs because they are included elsewhere.  However,
  // non-AT penalty is calculated here for the two "interior" pairs.
	
	
  int a, d, c, f, e, b; // 4 pairs: i-e, a-d, b-j, c-f
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE min_energy = NAD_INFINITY;
  DBL_TYPE tempMin;
	
  if( j - i <= 4) {// not enough room for pk
    return min_energy;
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
					
	  tempMin = 
	    Fg[ gap_index( i, a, d, e, seqlength) ]+
	    Fg[ gap_index( b, c, f, j, seqlength) ]+
	    bp_penalty + 2*BETA_2  +
	    Fz[ pf_index( e + 1, f - 1, seqlength)]+
	    Fz[ pf_index( c + 1, d - 1, seqlength)]+
	    Fz[ pf_index( a + 1, b - 1, seqlength)];
					
	  min_energy = MIN( min_energy, tempMin);	  
	}
      }
    }
  }
	
  if( j - i <= 6) {// not enough room for pk
    return min_energy;
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
									
		  tempMin = 
		    Fg[ gap_index( i, a, d, e, seqlength) ]+
		    Fg[ gap_index( b, c, f, j, seqlength) ]+
		    bp_penalty + 2*BETA_2 +
		    Fz[ pf_index( e + 1, f - 1, seqlength)]+
		    Fz[ pf_index( c + 1, d - 1, seqlength)]+
		    Fz[ pf_index( a + 1, b - 1, seqlength)];
									
		  min_energy = MIN( tempMin, min_energy);
									
									
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
									
		  tempMin = 
		    Fg[ gap_index( i, a, d, e, seqlength) ]+
		    Fg[ gap_index( b, c, f, j, seqlength) ]+
		    bp_penalty + 2*BETA_2 +
		    Fz[ pf_index( e + 1, f - 1, seqlength)]+
		    Fz[ pf_index( c + 1, d - 1, seqlength)]+
		    Fz[ pf_index( a + 1, b - 1, seqlength)];
		  min_energy = MIN( tempMin, min_energy);
									
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
	
  if( j - i <= 7) {// not enough room for pk
    return min_energy;
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
											
		      tempMin = 
			Fg[ gap_index( i, a, d, e, seqlength) ]+
			Fg[ gap_index( b, c, f, j, seqlength) ]+
			(bp_penalty + 2*BETA_2) + 
			Fz[ pf_index( e + 1, f - 1, seqlength)]+
			Fz[ pf_index( c + 1, d - 1, seqlength)]+
			Fz[ pf_index( a + 1, b - 1, seqlength)];
		      min_energy = MIN( min_energy, tempMin);
											
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
	
  return min_energy;
}
#endif //O_N8


