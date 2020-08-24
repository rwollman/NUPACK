/*
  complxesUtils.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 7/2006 and Justin Bois 1/2007

   Utilities for use with the program complexes.c.
*/

#include "complexesUtils.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

//Global variables
int nStrands;


/* **************** */
void resetNicks( int nStrands, int nicks[]) {

  int i;
  for(i = 0; i <= nStrands-1; i++) {
    nicks[i] = -1;
  }
}


/* **************** */
void nextMultiset( int nStrands, int *oldCode, int *newCode, int *seqCount,
                   int maxComplexSize) {
  //this generates the next unique multiset given the current one.
  //i.e., loopless multiset generation.
  //The multisets are all possible complexes up to size maxComplexSize,
  //containing any number of each of nStrands possible strand types.

  int i;
  int state;
  /*
     state = 1 means increment as if counting
     state = 0 means the sequence is finished
     state = -1 means zero-out first non-zero term, then set state to 1
  */

  if( *seqCount < maxComplexSize) state = 1;
  else if( *seqCount == maxComplexSize) state = -1;
  else {
    printf("Error with seqCount!!!\n");
    exit(1);
  }

  for( i = 0; i <= nStrands - 1; i++) {
    newCode[i] = 0;
    if( state == 0) {
      newCode[i] = oldCode[i];
      //just copy
    }
    else if( state == 1 && oldCode[i] < maxComplexSize) {
      newCode[i] = oldCode[i]+1;
      state = 0;
      (*seqCount)++;
    }
    else if( state == -1 && oldCode[i] != 0) {
      *seqCount -= oldCode[i];
      state = 1;
      //set current position to 0 and add 1 to next available term
    }
    else if( state != -1) {
      printf("Error in generating nextMultisets!  %d %d %d!\n",
             state, oldCode[i], maxComplexSize);
      exit(1);
    }
  }

  if( state != 0) {
    printf("Error with states in generating multisets!\n");
    exit(1);
  }
}


/* ***************** */
int isCyclicP( int size, int *a, int *b) {
  //determine if a and b are cyclic permutations of one another

  int i,j;


  i = 0;
  while( i <= size-1) {
    for( j = 0; j <= size-1; j++) {
      if( a[j] != b[(i+j) % size]) {
        break;
      }
    }
    if( j == size) return TRUE;
    i++;
  }


  return FALSE;
}

/* ************************* */
void symmetryCheck( multiset *allSets, int i, permutation *pm) {
  int j, k;
  int isSym; //does the permutation have a rotational symmetry?
  int isSymJ; //is the rotational symmetry N/J?

  pm->symmetryFactor = 1;
  isSym = FALSE;
  j = 1;
  while( j <= allSets[i].nSeqs-1 && !isSym) {
    //j is distance between symmetric positions
    if( allSets[i].nSeqs % j == 0) {//only possible if j | nSeqs
      isSymJ = TRUE; //assume true until shown otherwise
      k = 0;
      while( k <= allSets[i].nSeqs - 1 && isSymJ) {
        if( (pm->code)[k] !=
            (pm->code)[ (k + j) % allSets[i].nSeqs]) {
          isSymJ = FALSE;
        }
        k++;
      }
      if( isSymJ) {
        pm->symmetryFactor = allSets[i].nSeqs/j;
        isSym = TRUE;
      }
      //else check next case
    }
    //else j is not a legal symmetry step
    j++;
  }
}


/* **************************** */
void printPerms( FILE *fp2, int id, int nStrands,
                 multiset *set) {

  int j;
  permutation *currentPerm = set->perms;
  int permId = 1;

  while( currentPerm != NULL) {
    fprintf( fp2, "%d\t%d\t", id, permId++); //complex id

    for( j = 0; j <= set->nSeqs-1; j++) {
      fprintf( fp2, "%d\t", currentPerm->code[j]);
    }

    currentPerm = currentPerm->next;
    fprintf( fp2, "\n");

  }
}


/* ********************** */
int compareMultisets( const void *p1, const void *p2) {

  //uses document wide variable, nStrands
  //used for qsort.  Sorts first by size of Complex,
  //then by lexicographic order

  const multiset* ms1 = (multiset *)p1;
  const multiset* ms2 = (multiset *)p2;
  int i = 0;

  if( ms1->nSeqs < ms2->nSeqs) return -1;
  else if( ms1->nSeqs > ms2->nSeqs) return 1;

  for( i = 0; i < nStrands; i++) {
    if( (ms1->code)[i] > (ms2->code)[i]) return -1;
    else if( (ms1->code)[i] < (ms2->code)[i]) return 1;
  }

  return 0;
}

/* *********************** */
int comparePermutations( const void * p1, const void * p2) {
  const permutation * pm1 = (permutation *)p1;
  const permutation * pm2 = (permutation *)p2;

  if(pm1->nSeqs < pm2->nSeqs) return -1;
  else if(pm1->nSeqs > pm2->nSeqs) return 1;

  int i = 0;
  int num_strands = pm1->nSeqs;
  int tot1 = 0;
  while(tot1 < num_strands && pm1->strand_sums[i] == pm2->strand_sums[i]) {
    tot1 += pm1->strand_sums[i];
    i++;
  }
  if(tot1 < num_strands) {
    if(pm1->strand_sums[i] < pm2->strand_sums[i]) {
      return 1;
    } else if(pm1->strand_sums[i] > pm2->strand_sums[i]) {
      return -1;
    }
  }


  i = 0;
  while(i < num_strands && pm1->code[i] == pm2->code[i]) {
    i++;
  }
  if(i < num_strands && pm1->code[i] < pm2->code[i]) {
    return -1;
  } else if(i < num_strands && pm1->code[i] > pm2->code[i]) {
    return 1;
  }
  return 0;
}


/* *************** */
void printMfesToFile( const dnaStructures *ds, FILE *fp,
                      const int *nicks) {

  int i, j;
  int nickIndex;
  int npairs; // Number of pairs
  int **pairlist; // Each row is i,j pair

  if( ds->nStructs <= 0) {
    return;
  }

  // Allocate memory for pairlist (this is more than we need, but be safe)
  pairlist = (int **) malloc(ds->seqlength * sizeof(int *));
  for (i = 0; i < ds->seqlength; i++) {
    pairlist[i] = (int *) malloc(2 * sizeof(int));
  }

  // Print the free energy
  fprintf(fp, "%.8Le\n", (long double) (ds->validStructs)[0].correctedEnergy);

  // Prints MFE struct on the fly and then pairs
  for( i = 0; i < ds->nStructs; i++) {
    if (i > 0) { // New record
      fprintf(fp,"\n%% %%%%%%%%%% Next degenerate MFE structure %%%%%%%%%%\n");
    }
    nickIndex = 0;
    npairs = 0;
    for( j = 0; j < ds->seqlength; j++) {
      if( (ds->validStructs)[i].theStruct[j] > j) {
        fprintf(fp, "(");
        pairlist[npairs][0] = j;
        pairlist[npairs++][1] = (ds->validStructs)[i].theStruct[j];
      }
      else if(  (ds->validStructs)[i].theStruct[j] == -1 ) {
        fprintf(fp, ".");
      }
      else fprintf(fp, ")");

      if( nicks[ nickIndex] == j) {
        fprintf( fp, "+");
        nickIndex++;
      }
    }
    fprintf( fp, "\n");

    // Print pairs
    for (j = 0; j < npairs; j++) {
      fprintf(fp,"%d\t%d\n",pairlist[j][0]+1,pairlist[j][1]+1);
    }
  }

  for (i = 0; i < ds->seqlength; i++) {
    free(pairlist[i]);
  }
  free(pairlist);

}
