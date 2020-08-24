/*
  permBG.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 10/2004, modified from the source described below.
*/

/*======================================================================*/
/* C program for distribution from the Combinatorial Object             */
/* Server. Generate permutations of a multiset in Gray code order.      */
/* This is the same version used in the book "Combinatorial Generation."*/
/* The program can be modified, translated to other languages, etc.,    */
/* so long as proper acknowledgement is given (author and source).      */  
/* Programmer: Frank Ruskey, 1995.                                      */
/* Programmer: Joe Sawada, 1997 (translation to C).                     */
/* The latest version of this program may be found at the site          */
/* http://www.theory.csc.uvic.ca/~cos/inf/mult/Multiset.html            */
/*======================================================================*/

// A direct CAT implementation of the Eades-McKay algorithm
// modified to handle multisets.                           

// This program takes as input t (the number of different elements 0..t )
//	followed by the number of each element n(0).. n(t)  

// modified by Robert Dirks 10/08/04
// for incorporation into my multistranded pf code

#include "permBG.h"
#include <stdlib.h>
#include <stdio.h>
#include "complexesUtils.h"


//global variables
int n,t,cnt,calls;
int a[100], num[20], permsum[20], off[20];
int dir[20];

int **allPerms = NULL;
int *integer = NULL;
int nPrint;

/* ********* */

void PrintIt() {
  //"print to an array"

  int i;

    
  allPerms[cnt] = (int*) malloc( permsum[0]*sizeof(int));
  
  for(i=1; i<= permsum[0]; i++) {
    allPerms[cnt][ i-1] = integer[ a[i]];
  }
  
  cnt++;

}

/* ************ */

void BigGen( int t) {

	if (dir[t]) gen(t,permsum[t-1],permsum[t]);
	else neg(t,permsum[t-1],permsum[t]);
	
	if (t>1) BigGen(t-1);
	dir[t] = (dir[t] + 1)%2;
	if (dir[t]) off[t] = off[t-1] + num[t-1];
	else off[t] = off[t-1];

}

/* ***************** */
	
void swap(int t, int x, int y) {

	int b,temp;

	if (t>1) BigGen(t-1);
	b = off[t-1];
	temp = a[x+b];
	a[x+b] = a[y+b];
	a[y+b] = temp;
	PrintIt();

}
/* **************** */
	
void gen( int t, int n, int k) { // This is the Eades-McKay algorithm

	int i;

	calls++;
	if ((1<k) && (k<n)) {
     		gen( t, n-2, k-2 );  swap( t, n-1, k-1 );
     		neg( t, n-2, k-1 );  swap( t,  n, n-1 );
     		gen( t, n-1, k );
	}
	else if (k == 1) {
		for (i=n-1; i>= 1; i--) swap(t,i,i+1);
	}
}

/* ***************** */
		
void neg( int t, int n, int k) {

	int i;

	calls++;
	if ((1<k) && (k<n)) {
     		neg( t, n-1, k );    swap( t,  n, n-1 );
     		gen( t, n-2, k-1 );  swap( t, n-1, k-1 );
     		neg( t, n-2, k-2 );  
	}
	else if (k == 1) {
		for (i=1; i<=n-1; i++) swap(t,i,i+1);
	}
}

/* ************** */

void initializeMP( int size, int *multiPerm) { 
  //computes all multiset permutations

  int i,j;
  int lastNum;
  
  nPrint = 1; //the first permutation is already set (multiPerm)
  
  t = 0;
  lastNum = multiPerm[0];
  for( i = 1; i <= size - 1; i++) {
    if( multiPerm[i] != lastNum) {
      lastNum = multiPerm[i];
      t++;
    }
  }
  
  integer = (int*) malloc( (t+1)*sizeof(int));
  
  j = 0;
  for( i = 0; i <= t; i++) {
    dir[i] = 1;
    num[i] = 0;
    integer[i] = lastNum = multiPerm[j];
    while( j <= size-1 && multiPerm[j] == lastNum) {
      a[j+1] = i;
      num[i]++;
      j++;
    }
  }
  
  //allocate memory
  allPerms = (int **) malloc( ((int) factorial(size))*
			      sizeof( int*) );
  
  off[0] = 0;
  for(i=1; i<=t; i++) off[i] = off[i-1] + num[i-1];
  dir[t+1] = 1;
  permsum[t] = num[t];
  for(i=t-1; i>=0; i--)  permsum[i] = permsum[i+1] + num[i];
  
  cnt = 0;  calls = 0;  
  PrintIt();
  BigGen( t );
}

/* ***************** */
int nextPerm() {
  
  if( nPrint == cnt) {
    freeMP();
    return 0;
  } 

  nPrint++;
  return 1;
}

// if location is NULL, return the number of permutations but save nothing.
int makeFCPermutations(permutation * location, int * content, int length, int nStrands) {
  // Generate all permutations and check for symmetry 
  int * current_code = (int *) calloc(length, sizeof(int));
  int offset = 0;
  int permi = 0;
  int strandi = 0;
  int strandj = 0;
  int codei = 0;
  int codek = 0;
  int codel = 0;

  permutation * head = NULL; // points to the head of a temporary linked list for counting
  permutation * curPerm = NULL;

  //location[offset].nSeqs = length;
  //location[offset].code = (int*) calloc(length,sizeof(int));
  //location[offset].strand_sums = (int*) calloc(nStrands,sizeof(int));
  strandj = 0;
  for(strandi = 0 ; strandi < nStrands ; strandi++) {
    // location[offset].strand_sums[strandi] = content[strandi];
    for(codei = 0 ; codei < content[strandi] ; codei++) {
      if(strandj >= length) {
        fprintf(stderr,"Internal error: too many strands in complex\n");
        exit(1);
      }
      current_code[strandj] = 1 + strandi;
      strandj++;
    }
  }

  if(location != NULL) {
    location[offset].nSeqs = length;
    location[offset].code = (int*)malloc(length * sizeof(int));
    location[offset].strand_sums = (int*) malloc(nStrands * sizeof(int));
    location[offset].symmetryFactor = 1;
    for(codei = 0; codei < length ; codei++) {
      location[offset].code[codei] = current_code[codei];
    }
    for(strandi = 0; strandi < nStrands ; strandi++) {
      location[offset].strand_sums[strandi] = content[strandi];
    }
    offset ++;
  } else {
    head = (permutation *) malloc(sizeof(permutation));
    curPerm = head;
    curPerm->nSeqs = length;
    curPerm->code = (int*)malloc(length * sizeof(int));
    for(codei = 0; codei < length ; codei++) {
      curPerm->code[codei] = current_code[codei];
    }
    curPerm->next = NULL;
    offset ++;
  }


  while(1) {
    // Generate all permutations of the set
    codek = -1;
    for(codei = 1 ; codei < length - 1; codei++) {
      if(current_code[codei] < current_code[codei + 1]) {
        codek = codei;
      }
    }
    if(codek == -1) {
      break;
    }
    codel = codek + 1;
    for(codei = codek + 1; codei < length; codei++) {
      if(current_code[codek] < current_code[codei]) {
        codel = codei;
      }
    }
    strandi = current_code[codek];
    current_code[codek] = current_code[codel];
    current_code[codel] = strandi;

    codel = length - 1;
    int distance = (length - codek - 1) / 2;
    for(codei = 0 ; codei < distance ; codei++) {
      strandi = current_code[codek + 1 + codei] ;
      current_code[codek + 1 + codei] = current_code[length - codei - 1] ;
      current_code[length - codei - 1] = strandi;
    }
    int present = 0;
    if(location != NULL) {
      for(permi = 0 ; permi < offset ; permi++) {
        if(!present && isCyclicP(length,current_code,location[permi].code)) {
          present = 1;
          break;
        }
      }
      if(!present ) {
        location[offset].nSeqs = length;
        location[offset].code = (int*)malloc(length * sizeof(int));
        location[offset].strand_sums = (int*) malloc(nStrands * sizeof(int));
        for(codei = 0; codei < length ; codei++) {
          location[offset].code[codei] = current_code[codei];
        }
        for(strandi = 0; strandi < nStrands ; strandi++) {
          location[offset].strand_sums[strandi] = content[strandi];
        }
        offset++;
      }
    } else {
      curPerm = head;
      while(NULL != curPerm) {
        if(!present && isCyclicP(length,current_code,curPerm->code)) {
          present = 1;
          break;
        }
        curPerm = curPerm->next;
      }
      if(!present) {
        curPerm = head;
        head = (permutation*)malloc(sizeof(permutation));
        head->code = (int*)malloc(length * sizeof(int));
        for(codei = 0 ; codei < length ; codei++) {
          head->code[codei] = current_code[codei];
        }
        head->next = curPerm;
        offset++;
      }
    }
  }
  if(head) {
    while(NULL != head) {
      curPerm = head;
      head = head->next;
      free(curPerm->code);
      free(curPerm);
    }
  }
  free(current_code);
  return offset;
}


int CountSets(permutation * perms,int nPerms, int nStrands) {
  int permi;
  int permj;
  int strandi;
  int nSets = 1;
  for(permi = 0 ; permi < nPerms - 1 ; permi++) {
    permj = permi + 1;
    if(perms[permi].nSeqs != perms[permj].nSeqs) {
      nSets++;
    } else {
      for(strandi = 0 ; strandi < nStrands ; strandi++) {
        if(perms[permi].strand_sums[strandi] != 
            perms[permj].strand_sums[strandi]) {
          nSets++;
          break;
        }
      }
    }
  }
  return nSets;
}

int FillSets(multiset * allSets, permutation * allPermutations, 
              int totalSets, int totalPerms, 
              int nStrands, int * seqlength) {
 
  int seti = 0; // set index
  int permi = 0;// permutation index
  int strandi = 0; // strand index (used as strand id and ordered index)
  int basei = 0; // base index
  int basej = 0;
  int maxLength = 0;
  for( seti = 0; seti < totalSets; seti++) {
    allSets[seti].code = (int*) malloc(nStrands*sizeof(int));
    allSets[seti].nSeqs = allPermutations[permi].nSeqs;
    allSets[seti].totalLength = 0;
    for(strandi = 0; strandi < nStrands ; strandi++) {
      allSets[seti].code[strandi] = allPermutations[permi].strand_sums[strandi];
      allSets[seti].totalLength += seqlength[strandi] * allSets[seti].code[strandi];
    }
    if(allSets[seti].totalLength + allSets[seti].nSeqs > maxLength) {
      maxLength = allSets[seti].totalLength + allSets[seti].nSeqs;
    }
    allSets[seti].nPerms = 0;
    allSets[seti].perms = allPermutations + permi;
    
    int reached_new_multiset = 0;
    do {
      allSets[seti].nPerms++;
      allPermutations[permi].next = NULL;
      allPermutations[permi].baseCode = (int*)malloc(allSets[seti].totalLength * 2 * sizeof(int));
      allPermutations[permi].seq = (char *)malloc((allSets[seti].totalLength 
                                  + allSets[seti].nSeqs)*sizeof(char));
      allPermutations[permi].symmetryFactor = 1;

      basej = 0;
      for(strandi = 0; strandi < allPermutations[permi].nSeqs; strandi++) {
        for(basei = 0; basei < seqlength[allPermutations[permi].code[strandi] - 1]; basei++) {
          allPermutations[permi].baseCode[2 * basej] = allPermutations[permi].code[strandi] - 1;
          allPermutations[permi].baseCode[2 * basej + 1] = basei;
          basej++;
        }
      }

      symmetryCheck(allSets,seti,allPermutations + permi);
      permi++;
      // We are out of permutations
      if(permi >= totalPerms) {
        reached_new_multiset = 1;
      }
      // There are more sequences
      if(!reached_new_multiset && allPermutations[permi].nSeqs != allSets[seti].nSeqs) {
        reached_new_multiset = 1;
      }
      // Or the sequences are different
      for(strandi = 0;!reached_new_multiset && strandi < nStrands ; strandi++) {
        if(allPermutations[permi].strand_sums[strandi] != allSets[seti].code[strandi]) {
          reached_new_multiset = 1;
        }
      }
      
      if(!reached_new_multiset) {
        allPermutations[permi - 1].next = allPermutations + permi;
      }
    } while (!reached_new_multiset);

#ifdef DEBUG
    printf("Set %d: ", seti);
    int j;
    int k;
    permutation * currentPerm;
    for( j = 0; j <= nStrands - 1; j++) {
      printf( "%d",allSets[seti].code[j]);
    }
    printf("\nTotal permutations: %d\n", allSets[seti].nPerms);
    printf("Cumulative perms: %d\n", totalPerms);

    currentPerm = allSets[seti].perms;
    while( currentPerm != NULL) {
      printf("Symmetry Factor for following code: %d\n",
             currentPerm->symmetryFactor);
      for( k = 0; k <= allSets[seti].nSeqs - 1; k++) {
        printf("%d", (currentPerm->code)[k]);
      }
      printf("\n");
      /*
      for( k = 0; k < allSets[i].totalLength; k++) {
        printf("%d %d\n", currentPerm->baseCode[2*k],
        currentPerm->baseCode[2*k+1]);
      }
      */
      currentPerm = currentPerm->next;
    }
#endif
  }
  return maxLength;

}

int GetMaxComplexSize(multiset * allSets, int totalSets) {
  int seti;
  int maxSize = 0;
  for(seti = 0 ; seti < totalSets ; seti++) {
    if(allSets[seti].nSeqs > maxSize) {
      maxSize = allSets[seti].nSeqs;
    }
  }
  return maxSize;
}

/* ************* */
// return the number of permutations generated
int makePermutations(permutation * location, int length, int nStrands) {
  int * current_code = (int *) calloc(length + 1, sizeof(int));
  int offset = 0;
  int code_counter = 0;
  int i = length;
  int j = 0;
  int test = nStrands > 1;

  location[offset].nSeqs = length;
  location[offset].code = (int*) calloc(length,sizeof(int));
  location[offset].strand_sums = (int*) calloc(nStrands,sizeof(int));
  for(code_counter = 0 ; code_counter < length ; code_counter++) {
    int current_strand = current_code[code_counter + 1] + 1;
    location[offset].code[code_counter] = current_strand;
    location[offset].strand_sums[current_strand - 1] ++;
  }
  offset ++;
  
  while(test) {
    current_code[i] = current_code[i] + 1;
    for(j = 1 ; j <= length - i ; j++) {
      current_code[i + j] = current_code[j];
    }
    if(length % i == 0) {
      location[offset].nSeqs = length;
      location[offset].code = (int *) calloc(length, sizeof(int));
      location[offset].strand_sums = (int *)calloc(nStrands,sizeof(int));
      for(code_counter = 0; code_counter < length ; code_counter++) {
        int current_strand = current_code[code_counter + 1] + 1;
        location[offset].code[code_counter] = current_strand;
        location[offset].strand_sums[current_strand - 1] ++;
      }
      offset++;
    }
    i = length;
    while(current_code[i] == nStrands - 1) {
      i = i - 1;
    }
    test = i != 0;
  }
  free(current_code);
  return offset;
}

/* ************* */

void setPerm( int *perm) {
  
  int i;

  if( nPrint <= 0 || nPrint >= cnt+1) {
    printf("Error, setPerm called illegally!\n");
    exit(1);
  }

  //printf("%d %d\n", nPrint-1, permsum[0]-1);
  for(i=0; i<= permsum[0]-1; i++) {
    perm[i] = allPerms[nPrint-1][ i];
  }

}

/* ************** */
void freeMP() {
  
  int i;

  free( integer); integer = NULL;
  for( i = 0; i <= cnt - 1; i++) {
    free( allPerms[i]); allPerms[i] = NULL;
  }
  free( allPerms); allPerms = NULL;

}
  

