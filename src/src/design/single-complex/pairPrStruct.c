/*
  pairPrStruct.c
  
  Created by Robert Dirks on 6/12/06.

  Used for consensus structures.
 */
#include <stdio.h>
#include <stdlib.h>

// I'm not really sure if it needs this or not, but has useful stuff in it...
// #include "nupack.h" 
#include "design_pfunc_utils.h" 


//used to create the "greedy" consensus structure

int compare_PI( const void *a, const void *b) {
  
  double A = ((pairInfo *) a)->value;
  double B = ((pairInfo *) b)->value;

  if( A > B) return -1;
  if( A < B) return 1;
  return 0;
}

void makePairStruct( char *parens, DBL_TYPE *pairPr, int seqlength) {
  
  int pos, i,j, count;
  int size = seqlength*(seqlength+1)/2 + seqlength;
  pairInfo *pi = (pairInfo*) malloc( size*sizeof( pairInfo));
  char *tmpParens = (char*) malloc( (seqlength+1)*sizeof( char));
  int a, b;

  int nSet; //#of positions set
  
  if( pi == NULL) {
    printf("Unable to allocate pairInfo!");
    exit(-1);
  }
	
  //create array of pair info
  count = 0;
  for( i = 0; i < seqlength; i++) {
    for( j = i+1; j < seqlength; j++) {
      pi[ count].i = i;
      pi[ count].j = j;
      pi[ count++].value = pairPr[ i*(seqlength+1) + j];
    }
    pi[ count].i = i;
    pi[ count].j = -1;
    pi[ count++].value =  pairPr[ i*(seqlength+1) + seqlength];
  }      
	
  qsort( pi, count, sizeof( pairInfo), compare_PI); 
  //sort by pairPr in decreasing order
  
  //determine optimal structure
	
  //first initialize tmpParens
  tmpParens[seqlength] = '\0';
  for( i = 0; i < seqlength; i++) {
    tmpParens[ i] = '*';
  }
  
  //loop over structures in decreasing order
  nSet = 0; //# of bases fixed
  i = 0;
  while( i < count && nSet < seqlength) {
    a = pi[i].i;
    b = pi[i].j;
    if( b == -1 && tmpParens[a] == '*') {
      tmpParens[a] = '.';
      nSet++;
    }
    else if( tmpParens[a] == '*' && tmpParens[b] == '*') {
      tmpParens[a] = '(';
      tmpParens[b] = ')';
      nSet = nSet + 2;
    }
    i++;
  }
	
  
  //set parens appropriately
  pos = i = 0;
  while( pos < seqlength) {
    if( parens[i] != '+') {
      parens[i++] = tmpParens[pos++];
    }
    else {
      i++;
    }
  }
  
	
  free(pi); pi = NULL;
  free( tmpParens); tmpParens = NULL;
}


