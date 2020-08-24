/*
  sumexp_pk.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 9/2001, Justin Bois 1/2007
 
  This file contains the functions that calculate the sum of
  exponentials in order to calculate a given partition function
  matrix.  This deals with pseudoknot cases
  
  See pfuncUtilsHeader.h for more specific descriptions of each function.
*/

#include "sumexp_pk.h"

extern unsigned int seqHash;
/* ************************************** */
DBL_TYPE SumExpQb_Pk( int i, int j, int seq[], int seqlength, 
                     DBL_TYPE Qp[], DBL_TYPE Qm[] ) {
 // Determine all possible, rightmost pseudoknots 
 // contained in i,j "loop" and calculate Sumexpl
 
  DBL_TYPE bp_penalty = 0;
  DBL_TYPE sumexpl = 0.0;
  DBL_TYPE explTmp, alphaTmp;
  int d, e; //the left and right ends of a pseudoknot
  int pf_de;

  //int S1 = j-i+1;
  //int S2, S3;

  if( ( seq[i]) + ( seq[j]) != 5) {
   return 0;
  }

  if( seq[i] != BASE_C  && seq[j] != BASE_C) {
   bp_penalty = AT_PENALTY;
  }

  alphaTmp=BETA_1M + ALPHA_1 + 3*ALPHA_2 + bp_penalty;
  for( d = i+1; d <= j - 6; d++) {
    for( e = d + 5; e <= j - 1; e++) {
      pf_de = pf_index( d, e, seqlength);

      //S2 = e-d+1;
      //S3 = d-i-1;

      explTmp=ExplDangle( e+1, j-1, seq, seqlength);
      sumexpl += 
        //scale( S1 - S2) *
        EXP_FUNC( -( alphaTmp + 
               (j-e-1 + d-i-1)*ALPHA_3)/(kB*TEMP_K) ) *
        ExplDangle( i+1, d-1, seq, seqlength) * explTmp *
//        ExplDangle( e+1, j-1, seq, seqlength) *
        Qp[ pf_de];

      sumexpl += 
        //scale( S1 - S2 - S3) *
        Qm[ pf_index( i+1, d-1, seqlength)] * 
        Qp[ pf_de] * 
        EXP_FUNC( -( alphaTmp + 
               (j - e - 1)*ALPHA_3) / (kB*TEMP_K) ) * explTmp;
//        ExplDangle( e+1, j-1, seq, seqlength);
    }
  }
 return sumexpl;
}

/* *********************************************** */
void MakeQg_N5( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm, DBL_TYPE *Qgls, DBL_TYPE *Qgrs, DBL_TYPE *QgIx,
               DBL_TYPE *QgIx_2, short *possiblePairs) {

  //  Make the gap matrix for Qg
  int c,d,e,f;
  int gap_idej;
  DBL_TYPE IJ_bp_penalty = 0.0;
  DBL_TYPE bp_penalty = 0;

  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
   IJ_bp_penalty = AT_PENALTY;
  }


  if( CanPair( seq[i], seq[j]) == TRUE) {
   
   //Case 0: Only 1 pair
   Qg[ gap_index( i,i,j,j,seqlength)] = 1;
   
   //Case 1:  Terminal Inner Pair
   for( d = i+1; d <= j-5; d++) {
     for( e = d+4; e <= j-1; e++) {      
       if( CanPair( seq[d], seq[e] ) == TRUE) {
         
         Qg[ gap_index( i,d,e,j, seqlength)] += 
           EXP_FUNC( -InteriorEnergy( i, j, d, e, seq)/(kB*TEMP_K) ); 
         
       }
     }
   }
  }

  fastIloop_Qg(i, j, seq, seqlength, Qg, QgIx, QgIx_2, possiblePairs);

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
               
               Qg[ gap_idej] += 
                 Qm[ pf_index( i+1, d-1, seqlength)] *
                 EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + 
                         (j-e-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
                 ExplDangle( e+1, j-1, seq, seqlength);
               
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
               
               Qg[ gap_idej] += 
                 Qm[ pf_index( e+1, j-1, seqlength)] *
                 EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + 
                         (d-i-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
                 ExplDangle( i+1, d-1, seq, seqlength) ;
               
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
               
               Qg[ gap_idej] += 
                 Qm[ pf_index( i+1, d-1, seqlength)] *
                 Qm[ pf_index( e+1, j-1, seqlength)] *
                 EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + bp_penalty)/(kB*TEMP_K));
               
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
              
              Qg[ gap_idej] += 
                EXP_FUNC( -(ALPHA_1+ALPHA_2 + (j-f-1)*ALPHA_3+ bp_penalty)/
                     (kB*TEMP_K)) * 
                Qgls[ gap_index( i+1, d, e, f, seqlength)] *
                ExplDangle( f+1, j-1, seq, seqlength);
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
              
              Qg[ gap_idej] += 
                EXP_FUNC( -(ALPHA_1+ALPHA_2 + (c-i-1)*ALPHA_3 + 
                        bp_penalty)/(kB*TEMP_K)) * 
                Qgrs[ gap_index( c, d, e, j-1, seqlength)] *
                ExplDangle( i+1, c-1, seq, seqlength);
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
              
              Qg[ gap_idej] += 
                Qm[ pf_index( i+1, c-1, seqlength)] *
                EXP_FUNC( -(ALPHA_1+ALPHA_2+bp_penalty)/(kB*TEMP_K)) * 
                Qgrs[ gap_index( c, d, e, j-1, seqlength)];
            }
          }
        }
      }
    } // Can Pair


  }

/* ************************************** */
void fastIloop_Qg(int i, int j, int seq[], int seqlength, 
                  DBL_TYPE *Qg, DBL_TYPE *QgIx, DBL_TYPE *QgIx_2,
                  short *possiblePairs) {
                    
  int d, e, size;
  int L = j - i + 1;
  int pf_ij = pf_index(i, j, seqlength);
  int index;
  DBL_TYPE ExplInteriorMM = NAD_INFINITY;
  int gap_idej;
  
  if( CanPair( seq[i], seq[j]) == TRUE) {
    ExplInteriorMM = EXP_FUNC( -InteriorMM( seq[i], seq[j], seq[i+1], 
                                       seq[j-1])/(kB*TEMP_K));
  }
  
  //QgIx recurions
  //explicitly add in base cases to QgIx (smallest extensible loops) 
  
  if( L >= 17 && possiblePairs[pf_ij] == TRUE) {
    makeNewQgIx( i, j, seq, seqlength, Qg, QgIx);
  }

  for( d = i+1; d <= j-5; d++) {
    for( e = d+4; e <= j-1; e++) {
      if(  CanPair( seq[d], seq[e]) == TRUE) {

        gap_idej = gap_index( i,d,e,j, seqlength);

        if( L>=17 &&  CanPair( seq[i], seq[j]) == TRUE) { 

          index = QgIxIndex( j-i, i, 8, d, e, seqlength);

          for( size = 8; size <= L - 9; size++) {
            Qg[ gap_idej] += QgIx[ index] * 
              ExplInteriorMM;
            index = index + (L-2)*(L-3)/2;
          }
          
        }

        if( i != 0 && j != seqlength - 1) {
          extendOldQgIx( i, j, d, e, seq, seqlength, 
                        Qg, QgIx, QgIx_2);
        }

        Qg[ gap_idej] += SumexplInextensibleIL_Qg( i, j, d, e, seq, seqlength,
                                                  Qg);
      }
    }
  }
}




/* ************************************ */
void MakeQgls( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg, 
              DBL_TYPE *Qm,  DBL_TYPE *Qgls) {
  int c, d, e;
  int pf_ic1;
  DBL_TYPE bp_penalty = 0;

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
           Qgls[ gap_index(i, d, e, j, seqlength)] +=
             Qm[ pf_ic1] *
             Qg[ gap_index(c, d, e, j, seqlength)] *
             EXP_FUNC( -(ALPHA_2 + bp_penalty)/(kB*TEMP_K) );
         }
       }
     }
    }
  }
}

/* ************************** */
void MakeQgrs( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg, 
              DBL_TYPE *Qm, DBL_TYPE *Qgrs) {
                
  int d, e, f;
  int gap_idej;
  DBL_TYPE bp_penalty = 0;
  
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

           Qgrs[ gap_idej] +=
             Qm[ pf_index( f+1, j, seqlength)] * 
             Qg[ gap_index( i, d, e, f, seqlength)] *
             EXP_FUNC( -(ALPHA_2 + bp_penalty)/(kB*TEMP_K) );
          }
        }
      }
    }
  }
}

/* *************************************** */
void MakeQgl( int i, int j, int seq[], int seqlength, 
             DBL_TYPE *Qg, DBL_TYPE *Qgl, DBL_TYPE *Qz) {
               //make sure to call this AFTER MakeQg and BEFORE MakeQgr
               
  int d, e, f;
  int gap_idfj;
  DBL_TYPE bp_penalty = 0.0;

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
          Qgl[ gap_index( i, e, f, j, seqlength) ] +=
            Qg[ gap_idfj]*Qz[ pf_index( d+1, e, seqlength)] *
            EXP_FUNC( -(BETA_2 + bp_penalty)/(kB*TEMP_K) );
        }
      }
    }
  }
}

/* ******************************** */
void MakeQgr( int i, int j, int seq[], int seqlength, 
             DBL_TYPE *Qgr, DBL_TYPE *Qgl, DBL_TYPE *Qz) {
 //make sure to call this AFTER MakeQg and AFTER MakeQgl
 
  int d, e, f;
  int gap_idej;


  for( d = i+1; d <= j-3; d++) {
    for( e = d+2; e <= j-1; e++) {
      gap_idej =  gap_index(i, d, e, j, seqlength);
      for( f = e; f <= j-1; f++) {
       Qgr[ gap_idej ] +=
         Qgl[ gap_index(i,d,f,j,seqlength)]*
         Qz[ pf_index( e, f-1, seqlength)];
      }
    }
  }
}

/* ********************************** */
DBL_TYPE SumExpQp_N5( int i, int j, int seq[], int seqlength, 
                     DBL_TYPE *Qgl, DBL_TYPE *Qgr, DBL_TYPE *Qg, DBL_TYPE *Qz) {
                       
  int d, e, f;
  DBL_TYPE bp_penalty, new_bp_penalty;
  DBL_TYPE sumexpl = 0.0;

  int a, b, c;
  if( j - i <= 4) {
   return 0.0;
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

               sumexpl += 
                 Qg[ gap_index( i, a, d, e, seqlength) ] *
                 Qg[ gap_index( b, c, f, j, seqlength) ] *
                 EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                      (kB*TEMP_K) ) * 
                 Qz[ pf_index( e + 1, f - 1, seqlength)]*
                 Qz[ pf_index( c + 1, d - 1, seqlength)]*
                 Qz[ pf_index( a + 1, b - 1, seqlength)];
          }
        }
     }
  }


  if( j - i <= 6) {
   return sumexpl;
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

             sumexpl += Qg[ gap_index(i, i, e, f, seqlength)] *
               Qz[ pf_index( i+1, d-1, seqlength)] *
               Qgr[ gap_index(d, e-1, f+1, j, seqlength)] *
               EXP_FUNC( -(new_bp_penalty + BETA_2)/(kB*TEMP_K) );
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

               sumexpl += Qgl[ gap_index(i, d-1, e, f, seqlength)] *
                 Qg[ gap_index(d, d, j, j, seqlength)] *
                 Qz[ pf_index( d+1, e-1, seqlength)] *
                 Qz[ pf_index( f+1,j-1, seqlength)] *
                 EXP_FUNC( -(new_bp_penalty+BETA_2)/(kB*TEMP_K) );
          }
        }
      }
    }
  }


  if( j - i <= 7) {// not enough room for pk
   return sumexpl;
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

               sumexpl += Qgl[ gap_index(i, d-1, e, f, seqlength)] *
                 Qgr[ gap_index(d, e-1, f+1, j, seqlength)] *
                 EXP_FUNC( -(new_bp_penalty)/(kB*TEMP_K) );
          }
        }
      }
    }
  }

  return sumexpl;
}


/* *************************************** */
void MakeQ_Qm_Qz( int i, int j, int seq[], int seqlength,
                 DBL_TYPE *Q, DBL_TYPE *Qm, DBL_TYPE *Qz,
                 DBL_TYPE *Qb, DBL_TYPE *Qp) {
                   
  int d, e; // d - e is internal basepair or pk boundary
  DBL_TYPE bp_penalty;
  int pf_ij = pf_index(i, j, seqlength);
  int pf_de;

  //int S1 = j-i+1;
  //int S2, S3;

  Q[ pf_ij] = /*scale(S1)* */ExplDangle(i, j, seq, seqlength);  //Empty Graph


  if( i != 0 && j != seqlength - 1) { 
   Qz[ pf_ij] = 
     //scale(S1)*
     ExplDangle(i, j, seq, seqlength)*
     EXP_FUNC( -BETA_3*(j-i+1)/(kB*TEMP_K) );
  }

  for( d = i; d <= j - 4; d++) {
   //S2 = d-i;

    for( e = d + 4; e <= j; e++) {

     //S3 = e-d+1;
      if( CanPair( seq[d], seq[e]) == TRUE &&
        CanWCPair(seq[d], seq[e]) ) {
        bp_penalty = 0;
        if( seq[d] != BASE_C && seq[e] != BASE_C) {
          bp_penalty = AT_PENALTY;
        }

        pf_de = pf_index( d, e, seqlength);

        Q[ pf_ij] += 
          //	scale( S1 - S2 - S3) *
          Q[ pf_index(i, d-1, seqlength)] *
          Qb[ pf_de ] *
          EXP_FUNC( -bp_penalty/(kB*TEMP_K) ) *
          ExplDangle( e+1, j, seq, seqlength); 

        if( i != 0 && j != seqlength - 1) {

          Qm[ pf_ij] +=
            //scale( S1 - S3) *
            EXP_FUNC( -(ALPHA_2 + ALPHA_3*(d-i + j-e) + bp_penalty)/
                 (kB*TEMP_K) )*
            Qb[ pf_de] *
            ExplDangle( e+1, j, seq, seqlength) *
            ExplDangle( i, d-1, seq, seqlength);

          if( d >= i+5) {
            Qm[ pf_ij] += 
              //scale( S1 - S2 - S3) *
              Qm[ pf_index(i, d-1, seqlength)] *
              Qb[ pf_de] *
              EXP_FUNC( -(ALPHA_2 + ALPHA_3*(j-e) + bp_penalty)/
                   (kB*TEMP_K) )*
              ExplDangle( e+1, j, seq, seqlength);
          }

          Qz[ pf_ij] += 
            //scale( S1 - S2 - S3) *
            Qz[ pf_index(i, d-1, seqlength)]*
            Qb[ pf_de] * EXP_FUNC( -(BETA_2 + BETA_3*(j-e) + bp_penalty)/
                              (kB*TEMP_K) )*
            ExplDangle( e+1, j, seq, seqlength);
        }
      }
    }
  }


  for( d = i; d <= j - 5; d++) {
   //S2 = d-i;

    for( e = d + 5; e <= j; e++) {

      pf_de = pf_index( d, e, seqlength);
      //S3 = e-d+1;

      Q[ pf_ij] += 
       //scale( S1 - S2 - S3) *
       Q[ pf_index(i, d-1, seqlength)] *
       Qp[ pf_de ] *
       EXP_FUNC( -BETA_1/(kB*TEMP_K) ) *
       ExplDangle( e+1, j, seq, seqlength); 

      if( i != 0 && j != seqlength - 1) {

       Qm[ pf_ij] +=
         //scale( S1 - S3) *
         EXP_FUNC( -(BETA_1M + 2*ALPHA_2 + ALPHA_3*(d-i + j-e))/
              (kB*TEMP_K) )*
         Qp[ pf_de] *
         ExplDangle( e+1, j, seq, seqlength) *
         ExplDangle( i, d-1, seq, seqlength);

       if( d >= i+5) {
         Qm[ pf_ij] += 
           //scale( S1 - S2 - S3) *
           Qm[ pf_index(i, d-1, seqlength)] *
           Qp[ pf_de] *
           EXP_FUNC( -(BETA_1M + 2*ALPHA_2 + ALPHA_3*(j-e))/
                (kB*TEMP_K) )*
           ExplDangle( e+1, j, seq, seqlength);
       }

       Qz[ pf_ij] += 
         //scale( S1 - S2 - S3) *
         Qz[ pf_index(i, d-1, seqlength)]*
         Qp[ pf_de] * EXP_FUNC( -(BETA_1P + 2*BETA_2 + BETA_3*(j-e) )/
                           (kB*TEMP_K) )*
         ExplDangle( e+1, j, seq, seqlength);	
      }
    }
  }
}


/* ************************************************************ */
DBL_TYPE SumexplInextensibleIL_Qg( int i, int j, int d, int e, 
                                  int seq[], int seqlength, 
                                  DBL_TYPE *Qg) {
    /* This finds the minimum energy IL that has a special energy 
    calculation, i.e. small loops, bulge loops or GAIL case */
    
    /* This finds the minimum energy IL that has a special energy 
    calculation, i.e. small loops, bulge loops or GAIL case */

    DBL_TYPE energy;
    int c, f; //Internal pair.(c, f will be restricted to special cases)

    int L1, L2; //size parameters: L1 + L2 = size, L1 = c-i-1, L2 = j-f-1

    DBL_TYPE sumexpl = 0;

    if( CanPair( seq[i], seq[j]) == FALSE) { 
      return 0;
    }

    /* First consider small loops */
    // L1 <= 3, L2 <= 3
    for( L1 = 0; L1 <= MIN(3, d-i-2); L1++) {
      c = i + L1 + 1;
      for( L2 = 0; L2 <= MIN( 3, j-e-2); L2++) {      
        f = j - L2 - 1;
        if( CanPair( seq[c], seq[f]) == TRUE) {
          
          energy = InteriorEnergy( i, j, c, f, seq);
          
          sumexpl += EXP_FUNC( -energy/(kB*TEMP_K)) *
            Qg[ gap_index( c, d, e, f, seqlength)];
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
          sumexpl += EXP_FUNC( -energy/(kB*TEMP_K)) *
            Qg[ gap_index( c, d, e, f, seqlength)]; 
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
          sumexpl += EXP_FUNC( -energy/(kB*TEMP_K)) *
            Qg[ gap_index( c, d, e, f, seqlength)];
          
        }
      }
    }
    return sumexpl;
  }

/* *********************************************** */
void makeNewQgIx( int i, int j, int seq[], int seqlength, 
                 DBL_TYPE *Qg, DBL_TYPE *QgIx) {
                   
 /*Determine the new entries of QgIx(i,j,d,e,size) that are not extended 
 versions of QgIx(i+1, j-1, h1, e, size-2) */
 
 //extern DBL_TYPE loop37[];
 //extern DBL_TYPE asymmetry_penalty[];
 //extern DBL_TYPE max_asymmetry;
 
 DBL_TYPE energy;
 int c, f; 
 //Internal pair.(c, f will be restricted to special cases)
 
 int d,e; 
 // Internal pair of gap matrix (no restrictions);
 
 int size, L1, L2; //size parameters: L1 + L2 = size, L1 = c-i-1, L2 = j-f-1
 int asymmetry;
 //int asymmetry_index;
 int partial_index; //used to limit calls to QgIxIndex
 int len1_len2 = (j-i-1)*(j-i-2)/2;
 
 //Add in all the cases that are not an extended version of an
 //extensible case.
 
  for( d = i + 6; d <= j - 10; d++) {
    for( e = d + 4; e <= j-6; e++) {
      if( CanPair( seq[d], seq[e]) == TRUE) {
        partial_index = QgIxIndex( j-i, i, 0, d, e, seqlength);

        //Case 1:  L1 = 4, L2 >= 4;
        L1 = 4;
        c= i + L1 + 1;
        for( L2 = 4; L2 <= j - e - 2; L2++) {
         f = j - L2 - 1;
         size = L1 + L2;
         
         if( CanPair( seq[c], seq[f]) == TRUE)  {
           asymmetry = abs( L1 - L2);
           
           energy = asymmetryEfn( L1, L2, size);
           
           
           
           energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
           /*Exclude the i-j stacking energy here, just in case i-j 
           don't pair */
           
           QgIx[ partial_index+size*len1_len2 ] += 
             EXP_FUNC(-energy/(kB*TEMP_K))*
             Qg[ gap_index( c, d, e, f, seqlength)];
         }
        }

        //Case 2  L1 > 4, L2 = 4
        if( d >= i+7) {
          L2 = 4;
          f = j-L2-1;
          for( L1 = 5; L1 <= d-i-2; L1++) {
            c = i + L1 + 1;
            size = L1 + L2;

            if( CanPair( seq[c], seq[f]) == TRUE)  {
             asymmetry = abs( L1 - L2);
             energy = asymmetryEfn( L1, L2, size);
             
             energy += InteriorMM( seq[f], seq[c], seq[f+1], seq[c-1]);
             /*Exclude the i-j stacking energy here, just in case i-j 
             don't pair */
             
             
             QgIx[ partial_index + size*len1_len2] +=
               EXP_FUNC(-energy/(kB*TEMP_K))*
               Qg[ gap_index( c, d,e,f, seqlength)];
             
            }
          }
        }
      }
    }
  }
}


/* ************ for FASTILOOP ************** */
void extendOldQgIx( int i, int j, int d, int e, int seq[], int seqlength, 
                   DBL_TYPE *Qg, DBL_TYPE *QgIx, DBL_TYPE *QgIx_2) {
  /* Extends all entries of QgIx */

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
    QgIx_2[ partial_index2 + (size+2)*len10] = 
     QgIx[ partial_index + size*len_12]
     * sizeTerm[ size];    
  }
}
/* ************************************ */

/* *************************************** */
#ifdef O_N8

void MakeQg_N8( int i, int j, char seq[], int seqlength, DBL_TYPE Qg[],
               DBL_TYPE Qm[]) {
  //  Make the gap matrix for Qg
  int c,d,e,f;
  int gap_idej;
  DBL_TYPE IJ_bp_penalty = 0.0;
  DBL_TYPE bp_penalty = 0;

  int S1;
  int S2, S3, S4;

  if( CanPair( seq[i], seq[j]) == FALSE) {
    return;
  }

  IJ_bp_penalty = 0.0;
  if( seq[i] != BASE_C && seq[j] != BASE_C) {
    IJ_bp_penalty = AT_PENALTY;
  }

  //Case 0: Only 1 pair
  Qg[ gap_index( i,i,j,j,seqlength)] = 1*scale(2);


  //Case 1:  Terminal Inner Pair
  for( d = i+1; d <= j-5; d++) {
    for( e = d+4; e <= j-1; e++) {      
      if( CanPair( seq[d], seq[e] ) == TRUE) {	
       Qg[ gap_index( i,d,e,j, seqlength)] += 
         scale( j-e+1+d-i+1)*
         EXP_FUNC( -InteriorEnergy( i, j, d, e, seq)/(kB*TEMP_K) ); 
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

             S1 = j-e+1 + d-i+1;
             S2 = f-e+1 + d-c+1;

             Qg[ gap_idej] += 
               scale( S1 - S2) *
               EXP_FUNC( -InteriorEnergy( i, j, c, f, seq)/(kB*TEMP_K)) * 
               Qg[ gap_index( c, d, e, f, seqlength)];
            }
          }
        }
      }
    }
  }

  if( CanPair( seq[i], seq[j]) == TRUE && 
    CanWCPair(seq[i], seq[j])) {

      //Case 3:  Multiloop Left
      for( d = i+6; d <= j-5; d++) {
        S2 = d-i-1;
        
        for( e = d+4; e <= j-1; e++) {
          if( CanPair( seq[d], seq[e] ) == TRUE &&
          CanWCPair(seq[d], seq[e])) {
           gap_idej = gap_index(i,d,e,j,seqlength);
           
           bp_penalty = IJ_bp_penalty;
           if( seq[d] != BASE_C && seq[e] != BASE_C) {
             bp_penalty += AT_PENALTY;
           }

           S1 = j-e+1+d-i+1;
           Qg[ gap_idej] += 
             scale( S1 - S2) *
             Qm[ pf_index( i+1, d-1, seqlength)] *
             EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + 
                     (j-e-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
             ExplDangle( e+1, j-1, seq, seqlength);
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

           S2 = j-e-1;
           S1 = j-e+1+d-i+1;
           
           Qg[ gap_idej] += 
             scale( S1 - S2) *
             Qm[ pf_index( e+1, j-1, seqlength)] *
             EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + 
                     (d-i-1)*ALPHA_3 + bp_penalty)/(kB*TEMP_K)) *
             ExplDangle( i+1, d-1, seq, seqlength);
          }
        }
      }

      //Case 5: Multiloop Both Sides
      for( d = i+6; d <= j-10; d++) {
        S2 = d-i-1;

        for( e = d+4; e <= j-6; e++) {
          if( CanPair( seq[d], seq[e] ) == TRUE && 
             CanWCPair(seq[d], seq[e])) {
           gap_idej = gap_index(i,d,e,j,seqlength);
           
           bp_penalty = IJ_bp_penalty;
           if( seq[d] != BASE_C && seq[e] != BASE_C) {
             bp_penalty += AT_PENALTY;
           }

           S1 = j-e+1+d-i+1;
           S3 = j-e-1;
           
           Qg[ gap_idej] +=
             scale( S1 - S2 - S3) *
             Qm[ pf_index( i+1, d-1, seqlength)] *
             Qm[ pf_index( e+1, j-1, seqlength)] *
             EXP_FUNC( -(ALPHA_1 + 2*ALPHA_2 + bp_penalty)/(kB*TEMP_K));
          }
        }
      }

      //Case 6: Multiloop Left + More Qg
      for( d = i+7; d <= j-6; d++) {
        for( e = d+4; e <= j-2; e++) {
          gap_idej = gap_index(i,d,e,j,seqlength);
          if( CanPair( seq[d], seq[e] ) == TRUE) { 
            for( c = i+6; c <= d-1; c++) {
              S2 = c-i-1;
              
              for( f = e+1; f <= j-1; f++) {
                if( CanPair( seq[c], seq[f]) == TRUE &&
                   CanWCPair(seq[c], seq[f])) {

                bp_penalty = IJ_bp_penalty;
                if( seq[c] != BASE_C && seq[f] != BASE_C) {
                  bp_penalty += AT_PENALTY;
                }

                S3 = f-e+1+d-c+1;		
                S1 = j-e+1+d-i+1;
                Qg[ gap_idej] += 
                  scale( S1 - S2 - S3) *
                  Qm[ pf_index( i+1, c-1, seqlength)]*
                  EXP_FUNC( -(ALPHA_1+2*ALPHA_2+(j-f-1)*ALPHA_3 + bp_penalty)
                      /(kB*TEMP_K)) * 
                  Qg[ gap_index( c, d, e, f, seqlength)] *
                  ExplDangle( f+1, j-1, seq, seqlength);
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

                 S1 = j-e+1+d-i+1;
                 S2 = j-f-1;
                 S3 = f-e+1+d-c+1;
                 
                 Qg[ gap_idej] +=
                   scale( S1 - S2 - S3) *
                   Qm[ pf_index( f+1, j-1, seqlength)]*
                   EXP_FUNC( -(ALPHA_1+2*ALPHA_2+(c-i-1)*ALPHA_3+bp_penalty)
                        /(kB*TEMP_K)) * 
                   Qg[ gap_index( c, d, e, f, seqlength)] *
                   ExplDangle( i+1, c-1, seq, seqlength);
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
              S2 = c-i-1;
              for( f = e+1; f <= j-6; f++) {
                if( CanPair( seq[c], seq[f]) == TRUE &&
                   CanWCPair(seq[c], seq[f])) {

                  bp_penalty = IJ_bp_penalty;
                  if( seq[c] != BASE_C && seq[f] != BASE_C) {
                   bp_penalty += AT_PENALTY;
                  }

                  S1 = j-e+1+d-i+1;
                  S3 = j-f-1;
                  S4 = f-e+1+d-c+1;

                  Qg[ gap_idej] += 
                   scale( S1 - S2 - S3 - S4) *
                   Qm[ pf_index( i+1, c-1, seqlength)]*
                   Qm[ pf_index( f+1, j-1, seqlength)]*
                   EXP_FUNC( -(ALPHA_1+2*ALPHA_2+bp_penalty)
                        /(kB*TEMP_K)) * 
                   Qg[ gap_index( c, d, e, f, seqlength)];
                }
              }
            }
          }
        }
      }
    }
}


DBL_TYPE SumexplQp_N8( int i, int j, char seq[], int seqlength, DBL_TYPE Qg[], 
                       DBL_TYPE Qz[]) {
  // Add all pseudoknots with ends i, j.
  // No penalty for non-AT, terminal basepairs need to be added for the
  // h and s pairs because they are included elsewhere.  However,
  // non-AT penalty is calculated here for the two "interior" pairs,
  // i.e. i5-i3 and k5-k3
  // Dangles calculated in Qp_int_* matrixes
  
  int a, d, c, f, e, b; // 4 pairs: i-e, a-d, b-j, c-f
  DBL_TYPE bp_penalty = 0.0;
  DBL_TYPE sumexpl = 0.0;
  
  int S1 = j-i+1;
  int S2, S3, S4, S5, S6;
  
  if( j - i <= 4) {// not enough room for pk
    return 0.0;
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

                S2 = e-d+1+a-i+1;
                S3 = j-f+1+c-b+1;
                S4 = f-e-1;
                S5 = d-c-1;
                S6 = b-a-1;

                sumexpl += 

                  Qg[ gap_index( i, a, d, e, seqlength) ] *
                  Qg[ gap_index( b, c, f, j, seqlength) ] *
                  
                  EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                       (kB*TEMP_K) ) * 
                  Qz[ pf_index( e + 1, f - 1, seqlength)]*
                  Qz[ pf_index( c + 1, d - 1, seqlength)]*
                  Qz[ pf_index( a + 1, b - 1, seqlength)]*
                  scale( S1-S6-S5-S4-S3-S2);

              }
         }
       }
  }

  if( j - i <= 6) {// not enough room for pk
    return sumexpl;
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

                              S2 = e-d+1+a-i+1;
                              S3 = j-f+1+c-b+1;
                              S4 = f-e-1;
                              S5 = d-c-1;
                              S6 = b-a-1;

                              sumexpl += 
                                Qg[ gap_index( i, a, d, e, seqlength) ] *
                                Qg[ gap_index( b, c, f, j, seqlength) ] *
                                EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                                     (kB*TEMP_K) ) * 
                                Qz[ pf_index( e + 1, f - 1, seqlength)]*
                                Qz[ pf_index( c + 1, d - 1, seqlength)]*
                                Qz[ pf_index( a + 1, b - 1, seqlength)]*
                                scale( S1-S6-S5-S4-S3-S2);
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
                              
                              S2 = e-d+1+a-i+1;
                              S3 = j-f+1+c-b+1;
                              S4 = f-e-1;
                              S5 = d-c-1;
                              S6 = b-a-1;
                              
                              sumexpl += 
                                
                                Qg[ gap_index( i, a, d, e, seqlength) ] *
                                Qg[ gap_index( b, c, f, j, seqlength) ] *
                                EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                                     (kB*TEMP_K) ) * 
                                Qz[ pf_index( e + 1, f - 1, seqlength)]*
                                Qz[ pf_index( c + 1, d - 1, seqlength)]*
                                Qz[ pf_index( a + 1, b - 1, seqlength)]*
                                scale(S1-S2-S3-S4-S5-S6);
                              
                              
                            }
                       }
                     }
                }
              }
         }
    }
  }

  if( j - i <= 7) {// not enough room for pk
    return sumexpl;
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

                                  S2 = e-d+1+a-i+1;
                                  S3 = j-f+1+c-b+1;
                                  S4 = f-e-1;
                                  S5 = d-c-1;
                                  S6 = b-a-1;

                                  sumexpl += 
                                    Qg[ gap_index( i, a, d, e, seqlength) ] *
                                    Qg[ gap_index( b, c, f, j, seqlength) ] *
                                    EXP_FUNC( -(bp_penalty + 2*BETA_2)/
                                         (kB*TEMP_K) ) * 
                                    Qz[ pf_index( e + 1, f - 1, seqlength)]*
                                    Qz[ pf_index( c + 1, d - 1, seqlength)]*
                                    Qz[ pf_index( a + 1, b - 1, seqlength)]*
                                    scale(S1-S2-S3-S4-S5-S6);
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

  return sumexpl;
}
#endif //O_N8




