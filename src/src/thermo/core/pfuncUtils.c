/*
  pfuncUtils.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2001, Justin Bois 1/2007

  Generic utility functions used in partition function calculations.
*/

#include "pfuncUtils.h"
/* ********************************************** */
//Added 08/21/2001

int gap_index( int h, int r, int m, int s, int seqlength) {
  // index variable should actually be h..h1..m1..m, but i will use
  // the above indices for simplicity. 

  int n = seqlength;
  long int indexG;

  long int h2 = h*h;
  long int h3 = h2*h;
  long int h4 = h3*h;
  long int m2 = m*m;
  long int n2= n*n;
  long int n3 = n2*n;
  long int r2 = r*r;
  long int r3 = r2*r;

  extern long int maxGapIndex;

  if( h == r && m == s) { //new case for only 1 bp in gap matrix
    return maxGapIndex - 1;
  }

  if( h >= r || r >= m || m >= s || s >= seqlength) {
    fprintf(stderr, "Illegal call to Gap_index! %d, %d, %d, %d\n", h, r, m, s);
    exit(1);
  }

  // the following indexing formula imitates 4 nested for loops
  
  indexG = (-24 - 50*h - 35*h2 - 10*h3 - h4 - 36*m -12*m2 +
            12*n + 70*h*n + 30*h2*n + 4*h3*n + 24*m*n - 12*n2 -30*h*n2 -
            6*h2*n2 + 4*h*n3 + 44*r - 48*n*r + 12*n2*r + 
            24*r2 - 12*n*r2 +  4*r3 + 24*s)/24 ;

  if( indexG >= maxGapIndex) {
    fprintf(stderr, "Gap Index too large! %ld %ld %d %d %d %d %d\n", 
           indexG, maxGapIndex,
           h,r,m,s,n);
    exit(1);
  }
  return indexG;
}



/* ******************************************** */
int pf_index_old( int i, int j, int N) {
  /* Calculate index for partition function array 
  N = sequence length */

  static int seqlength  = -1;
  static int storeValues[(MAXSEQLENGTH)*(MAXSEQLENGTH)];
  static int ss2;

  int ind,value;

  if( seqlength == N) {
    ind = i*N+j;
    if( j != i-1) {
      value=storeValues[ind];
      if( value != -1) {
        return value;
      }
      else storeValues[ind] = value = ss2 - (N-i)*(N-i+1)/2+(j-i);
      return value;
    }
    else {
      return ss2 + i;
    }
  }
  else {
    seqlength = N;
    ss2 = N*(N+1)/2;
    int maxN=N*N;
    for( ind = 0; ind < maxN; ind++) storeValues[ind] = -1;
    if( j != i-1) {
      ind =i*N+j;
      storeValues[ ind] = ss2 - (N-i)*(N-i+1)/2+(j-i);
      return storeValues[ ind];
    }
    else {
      return ss2 + i;
    }
  }

  /*
  if( j == i - 1) {
  return N*(N+1)/2 + i;
  }
  
  if( i < 0 || j > N - 1 || j < i) {
  printf("Illegal partition function index!\n i = %d, j = %d, ind = %d\n",
  i, j, N*(N+1)/2 - (N - i)*( N - i + 1)/2 + (j - i));
  exit(1);
  }
  
  return N*(N+1)/2 - (N - i)*( N - i + 1)/2 + (j - i);
  //return (N+1-i)*(N-i)/ 2 + N+1-j;
  */
} 


/* **************************************** */
int QgIxIndex( int d, int i, int size, int h1, int m1, int N) {

  int indexG;
  int d1d2 = (d-1)*(d-2);
  int d5 = d-5;
  int h1_i_1 = h1-i-1;

  if( i < 0 || d < 0 || i + d > N-1 || size > d-6 || h1 <= i ||
     m1 >= i+d) {
       fprintf(stderr, "Error in QgIxIndex %d %d %d %d %d %d\n", 
              d, i, size, h1, m1, N);
       exit(1);
  }

  indexG = i*d5*d1d2/2 + size*d1d2/2 + 
    h1_i_1*(d-1) - h1_i_1*(h1 - i)/2 + m1 - h1 - 1;

  if( indexG >= (N-d)*d5*(d1d2/2) ) {
    fprintf(stderr,  "%d >= %d\n", indexG, (N-d)*(d-5)*(d1d2/2) );
    exit(1);
  }

  return indexG;
}

/* **************************************************** */
int Base2int( char base) {
  /* Convert A, C, G, T, U to 1,2,3,4 respectively */
  if( base == 'A') {
    return BASE_A;
  }
  else if ( base == 'C') {
    return BASE_C;
  }
  else if (base == 'G') {
    return BASE_G;
  }
  else if (base == 'T' || base =='U') {
    return BASE_U;
  }
  else if (base == '+') {
    return STRAND_PLUS;
  }
  else {
    fprintf(stderr, "Error in Converting base %c!\n", base);
    exit(1);
    return NAD_INFINITY; // never returns this
  }
}

/* **************************************************** */
char Int2base( int base) {
  /* Convert 1,2,3,4,5 to A, C, G, T, U respectively */
  if( base == BASE_A) {
    return 'A';
  }
  else if ( base == BASE_C) {
    return 'C';
  }
  else if (base == BASE_G) {
    return 'G';
  }
  else if (base == BASE_T) {
    return 'T';
  }
  else if (base == BASE_U) {
    return 'U';
  }
  else if (base == STRAND_PLUS) {
    return '+';
  }
  else {
    fprintf(stderr, "Error in Converting base %d!\n", base);
    exit(1);
    return '\0'; // never returns this
  }
}


/* *********************************************************** */ 
int convertSeq(char seqchar[], int seqnum[], int seqlength){
  int i;
  for (i=0; i<seqlength; i++)
    seqnum[i]=Base2int(toupper(seqchar[i]));
  seqnum[i]=-1;
  return -1;
}

/* *********************************************************** */ 
int printSeqNum(int seqnum[]){
  int i;
  for (i=0; seqnum[i]>0 && seqnum[i]<6; i++)
    printf("%c",Int2base(seqnum[i]));
  printf("\n");
  return -1;
}

/* *********************************************************** */ 
int GetMismatchShift( int base1, int base2) {
  /* base1 and base2 are basepaired. the returned value is needed to 
  index energy arrays
  */
  int shift;
  
  if( base1 == BASE_A) { /* this is for retrieving 
  proper mismatch energy */
    shift = 0;
  }
  else if( base1 == BASE_C) {
    shift = 1;
  }
  else if( base1 == BASE_G && base2 == BASE_C) {
    shift = 2;
  }
  else if( base1 == BASE_G && (base2 == BASE_T|| base2 == BASE_U) ) {
    shift = 4;
  }
  else if( (base1 == BASE_T|| base1 == BASE_U) && base2 == BASE_A) {
    shift = 3;
  }
  else if( (base1 == BASE_T|| base1 == BASE_U) && base2 == BASE_G) {
    shift = 5;
  }
  else {
//    printf("Error in GetMismatchShift. %c and %c don't pair!\n", 
//           Int2base (base1), Int2base (base2));
    shift=-1; // No longer exit...
  }
  return shift;
}

/* ************************************************** */
int GetPairType( int b) { //assume pair of b is the watson crick pair
  int shift;
  
  if( b == BASE_A) { 
    shift = 0;
  }
  else if( b == BASE_C) {
    shift = 1;
  }
  else if( b == BASE_G) {
    shift = 2;
  }
  else if( b == BASE_T|| b == BASE_U) {
    shift = 3;
  }
  else {
    fprintf(stderr, "Error in GetPairType: %d!\n",b);
    return -1;
  }
  return shift;
}


/* **************************************************** */
void CheckPossiblePairs( short **possiblePairs, int seqlength, int seq[]) {
  int i, j, b, pfi;
  int outerPair;
  *possiblePairs = (short *) calloc( seqlength*(seqlength+1)/2, 
                                    sizeof( short));

  if( *possiblePairs == NULL) {
    fprintf(stderr, "Error in calloc for possiblePairs!\n");
  }

  for( i = 0; i < seqlength - 4; i++) {
    for( j = i+4; j < seqlength; j++) {
      if( (*possiblePairs)[ pf_index(i, j, seqlength)] == 0) {
        outerPair = -1;
        for( b = 0; b <= MIN( i, seqlength - 1 - j); b++) {
          if( CanPair( seq[i-b], seq[j+b]) == TRUE  
             //&& CanWCPair(seq[i-b], seq[j+b]) 
             ) {
               outerPair = b;
          }
        }
        for( b = 0; b <= MIN( i, seqlength - 1 - j); b++) {
          pfi=pf_index(i-b, j+b, seqlength);
          if( b <= outerPair) {
            (*possiblePairs)[pfi] = TRUE;
          }
          else {
            (*possiblePairs)[pfi] = -1;
          }
        }
      }
    }
  }
}

/* ******************************************** */
int fbixIndexOld( int d, int i, int size, int N ) {
  if( d < 0 || i < 0 || i + d >= N) {
    fprintf(stderr,  "Error in fbixIndex %d %d %d %d\n", d, i, size, N);
  }
  if( i*(d-1) + size > (N-d)*(d-1) ) {
    fprintf(stderr,   "Error in fbixIndex! %d > %d, %d %d %d %d,\n", i*(d-1) + size,
           (N-1-d)*(d-1), d,i,size,N);
  }

  return i*(d-1) + size;
}


/* *********************************** */
void LoadFold( fold *thefold, char filename[]) {

  FILE *fp;
  char line[MAXLINE];
  int i; //position in line
  int linenumber; //line number
  int mode;  // mode = 0 means sequence being read, mode = 1 means graph
  int init; // loop index for initializations
  int done; // has all the necessary info been read?
  int line_done; //has a given line been completely read?

  int paira, pairb, helixlength; // graph property variables
  int pk1, pk2;
  int start_fix;
  char *tmp_fix=NULL;// fixed bases
  int indx; // index for some loops
  int foldInitialized = FALSE;

  int nRead;
  int array[MAXLINE];
  char *token;
  int j;
  int nicks[MAXSEQLENGTH];

  fp = fopen(filename, "r");
  if( fp == NULL) {  // Make sure input file exits 
    fprintf(stderr, "Error opening file %s!\n", filename);
    exit(1);  
  }

  //init thefold
  thefold->seqlength = 0;
  thefold->nStrands = 1;

  linenumber = 0;
  mode = 0; //read in seqlength
  done = FALSE;
  while( fgets(line, MAXLINE, fp)!= NULL && done == FALSE) {  // Read lines
    linenumber++;
    line_done = FALSE;
    i = 0;

    while( line[i] != '\0' && line[0] != '>' && line[i] != '\n'
    && line_done == FALSE) {
      if( mode == 0) {  //read seqlength
        if( line[i] == '.') {
          mode = 1; // go to next line and start reading
          line_done = TRUE; //needed to exit loop
        }
        else {
          nRead = 0;
          token = strtok( line, " ");
          while( token != NULL) {
            if( sscanf( token, "%d", &(array[ nRead]) ) == 1) {
              //set thefold->nicks
              if( nRead == 0) nicks[0] = array[0] - 1;
              else nicks[nRead] = nicks[nRead-1] + 
                array[ nRead];
              thefold->seqlength += array[nRead];
              nRead++;
            }
            token = strtok( NULL, " ");
          }

          //allocate thefold->isNicked
          thefold->isNicked = (int *) calloc( thefold->seqlength, sizeof( int));

          //set thefold->isNicked
          for( j = 0; j <= nRead - 2; j++) {
            thefold->isNicked[ nicks[j] ] = 1;
          }
          thefold->nStrands = nRead; //set number of strands
          line_done = TRUE;
        }
      }
      else if(mode == 1) {  //Read in graph properties
        if( foldInitialized == FALSE) {
          thefold->pairs = 
            (int*) calloc( thefold->seqlength+1, sizeof(int));
          thefold->pknots = 
            (int*) calloc( thefold->seqlength+1, sizeof(int));

          thefold->fixedBases = 
            (int*) calloc( (thefold->seqlength)+1, sizeof(int));

          tmp_fix = (char*) calloc( (thefold->seqlength)+1, sizeof(char));

          for( init = 0; init <= thefold->seqlength; init++) {
            tmp_fix[init] = '0';
            thefold->pairs[init] = -1;
            thefold->pknots[init] = -1;
          }
          foldInitialized = TRUE;
        }

        if( sscanf( line, "p %d %d", &pk1, &pk2) == 2) {
          pk1--;
          pk2--;
          if( thefold->pknots[pk1] == -1 && thefold->pknots[pk2] == -1) {
            thefold->pknots[ pk1] = pk2;
            thefold->pknots[ pk2] = pk1;
            if( pk2 - pk1 < 8) {
              fprintf(stderr, "Error! Pseudoknot region too small! %d %d\n", 
                     pk1, pk2);
              exit(1);
            }
            line_done = TRUE;
          }
          else {
            fprintf(stderr, "Error! base %d or %d is the end of multiple pknots\n",
                   pk1, pk2);
            exit(1);
          }
        }
        else if( sscanf( line, "%d %d %d", &paira, &pairb, &helixlength) < 3
                && line[0] != '.') {
                  fprintf(stderr, "error: 3 integers expected at line %d\n", linenumber);
                  exit(1);   //error if format is incorrect
                }
        else if( line[0] != '.') {
          if( (paira > 0) && (pairb > paira) && (helixlength > 0) && 
             (pairb <= thefold->seqlength) ) {
               paira--;
               pairb--;
               if( thefold->pairs[ paira] == -1 
                  && thefold->pairs[ pairb] == -1 && pairb - paira >= 4) {

                    thefold->pairs[ paira] = pairb; 
                    thefold->pairs[ pairb] = paira;

                    for( indx = 1; indx <= helixlength - 1; indx++) {
                      if( thefold->pairs[ paira + indx] == -1 && 
                         thefold->pairs[ pairb -indx] == -1) {
                           
                           thefold->pairs[ paira + indx] = pairb - indx;
                           thefold->pairs[ pairb - indx] = paira + indx;	
                       }

                      else {
                        fprintf(stderr, "Base %d or %d assigned twice!\n",
                               paira + indx, pairb-indx);
                        exit(1);
                      }
                    }
                  }
               else {
                 fprintf(stderr, "Base %d or %d assigned twice!, or too close\n", 
                        paira, pairb);
                 exit(1);
               }
             }
          line_done = TRUE;
        }
        else {
          line_done = TRUE; //needed to exit loop
          mode = 2;
        }
      }
      else if( mode == 2) { //fixed base data
        if( line[0] == '.') {
          line_done = TRUE;
          done = TRUE;
        }
        else if( sscanf( line, "%d%s", &start_fix, tmp_fix) == 2) {
          indx = 0;
          start_fix -= 1;

          while( toupper( tmp_fix[ indx]) == 'A' ||
                toupper( tmp_fix[ indx]) == 'C' ||
                toupper( tmp_fix[ indx]) == 'T' ||
                toupper( tmp_fix[ indx]) == 'U' ||
                toupper( tmp_fix[ indx]) == 'G' ||
                toupper( tmp_fix[ indx]) == 'R' ||
                toupper( tmp_fix[ indx]) == 'Y' ||
                toupper( tmp_fix[ indx]) == 'M' ||
                toupper( tmp_fix[ indx]) == 'K' ||
                toupper( tmp_fix[ indx]) == 'S' ||
                toupper( tmp_fix[ indx]) == 'W' ||
                toupper( tmp_fix[ indx]) == 'V' ||
                toupper( tmp_fix[ indx]) == 'H' ||
                toupper( tmp_fix[ indx]) == 'B' ||
                toupper( tmp_fix[ indx]) == 'D' ||
                toupper( tmp_fix[ indx]) == 'N') {
                  if( start_fix + indx > thefold->seqlength) {
                    fprintf(stderr, "Base specification exceed sequence length!\n");
                    exit(1);
                  }
                  if( thefold->fixedBases[ start_fix + indx] != 0) {
                    fprintf(stderr,  "Attemping to specify base %d more than once!",
                           start_fix + indx);
                    exit(1);
                  }
                  if(  toupper( tmp_fix[ indx]) == 'A') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_A;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'C') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_C;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'G') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_G;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'T' ||
                          toupper( tmp_fix[ indx]) == 'U' ) {
                            thefold->fixedBases[ start_fix + indx] =  BASE_T;
                          }
                  else if(  toupper( tmp_fix[ indx]) == 'R') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_R;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'Y') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_Y;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'M') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_M;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'K') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_K;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'S') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_S;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'W') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_W;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'V') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_V;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'H') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_H;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'D') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_D;
                  }
                  else if(  toupper( tmp_fix[ indx]) == 'B') {
                    thefold->fixedBases[ start_fix + indx] =  BASE_B;
                  }
                  indx++;
                }
          line_done = TRUE;
        }
        else {
          fprintf(stderr, "Error in reading in fixed bases!");
          exit(1);
        }
      }
      else { //this should never happen
        fprintf(stderr, "Illegal mode value at line %d\n", linenumber);
        exit(1);
      }
    }
  }
  
  
  for( i = 0; i <= thefold->seqlength - 1; i++) {
    if( thefold->fixedBases[i] != 0 && thefold->pairs[i] != -1) {

      if( thefold->fixedBases[ i] <= 4) {
        if( thefold->fixedBases[ thefold->pairs[i] ] == 0) {
          thefold->fixedBases[ thefold->pairs[i] ] =
            5 - thefold->fixedBases[i];
        }
        else if( thefold->fixedBases[ thefold->pairs[i] ] !=
                5 - thefold->fixedBases[i] &&
                thefold->fixedBases[ thefold->pairs[i] ] !=
                7 - thefold->fixedBases[i]) {
                  fprintf(stderr, "Fixed bases %d %d at position %d and %d cannot pair!\n",
                         thefold->fixedBases[ thefold->pairs[i] ],
                         thefold->fixedBases[ i], i+1, thefold->pairs[i]+1);
                  exit(1);
                }
        //wobble pair OK?...sometimes
      }
      else if( thefold->fixedBases[i] == 5) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 10;
        } 
        else if( thefold->fixedBases[ thefold->pairs[i]] != 10) {
          fprintf(stderr, "1Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 6) {
        if( thefold->fixedBases[ thefold->pairs[i] ] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 9;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 9) {
          fprintf(stderr, "1aError with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 7) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 7;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 7) {
          fprintf(stderr, "1Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1); 
        }
      }
      else if( thefold->fixedBases[i] == 11) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 14;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 14) {
          fprintf(stderr, "2Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 10) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 5;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 5) {
          fprintf(stderr, "3Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 9) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 6;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 6) {
          fprintf(stderr, "4Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 8) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 8;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 8) {
          fprintf(stderr, "5Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 11) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 14;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 14) {
          fprintf(stderr, "6Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 12) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 13;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 13) {
          fprintf(stderr, "7Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 14) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 11;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 11) {
          fprintf(stderr, "8Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
      else if( thefold->fixedBases[i] == 13) {
        if( thefold->fixedBases[thefold->pairs[i]] == 0) {
          thefold->fixedBases[thefold->pairs[i]] = 12;
        } 
        else if( thefold->fixedBases[thefold->pairs[i]] != 12) {
          fprintf(stderr, "9Error with fixed bases %d and %d!\n", i, 
                 thefold->pairs[i]);
          exit(1);
        }
      }
    }
  }

  fclose( fp);
  if( done == FALSE) { //make sure file ended appropriately
    fprintf(stderr, "File terminated before final \".\" reached\n");
    exit(1);
  }

  free( tmp_fix);
}



/* **************** */
DBL_TYPE expectedCorrectBases( int *structPairs, int seqlength) {

  int i;
  DBL_TYPE value = 0;
  int pair;
  extern DBL_TYPE *pairPr;

  value = 0;
  for( i = 0; i< seqlength; i++) {
    pair = structPairs[i];
    if( pair == -1) {
      value += pairPr[ i*(seqlength+1) + seqlength]; //unpaired probability
    }
    else {
      value += pairPr[ i*(seqlength+1) + pair];
    }
  }

  return value;
}

/* ************ */
int checkSymmetry( const int *thepairs, int seqlength, const int *nicks, int possibleSymmetry, int nStrands) {
  int i, j;
  int nSeqsPerUnit, unitLength;
  int sym = 1;
  int isSym;
  int symPos;
  for( i = 2; i <= possibleSymmetry; i++) {
    if( possibleSymmetry %  i != 0) continue;
    nSeqsPerUnit = nStrands/i;
    unitLength = nicks[ nSeqsPerUnit-1] + 1;
    isSym = TRUE;
    for( j = 0; j < seqlength - unitLength; j++) {
      symPos = (j + unitLength);
      if( (thepairs[j] == -1 && thepairs[ symPos] != -1) ||
         (thepairs[j] != -1 && thepairs[ symPos] == -1) ) isSym = FALSE;
      if( thepairs[j] >= 0 && thepairs[symPos] >= 0 && 
         thepairs[j] != (thepairs[ symPos] - symPos + j + seqlength) % seqlength) isSym = FALSE;
      
      if( isSym == FALSE) break;
    }
    if( isSym) sym = i;
  }

  return sym;
}

/* ******************** */
void findUniqueMins( dnaStructures *ds, const int *nicks, int symmetry, 
                     int nStrands, DBL_TYPE minDev) {
  dnaStructures newDs = {NULL, 0, 0, 0, NAD_INFINITY};
  int i, j;
  int unique;
  int less_than;
  int nSeqsPerUnit = nStrands/symmetry;
  int unitLength = nicks[ nSeqsPerUnit-1] + 1;
  DBL_TYPE minE = (ds->validStructs)[0].correctedEnergy;

  newDs.validStructs = (oneDnaStruct*) malloc( (ds->nStructs)*sizeof(oneDnaStruct));
  newDs.nAlloc = ds->nStructs;
  newDs.seqlength = ds->seqlength;
  newDs.minError = ds->minError;
  newDs.nStructs = 0;

  for( i = 0; i < ds->nStructs; i++) {
    if( (ds->validStructs)[i].correctedEnergy > minE + ENERGY_TOLERANCE + minDev) break;
    unique = TRUE;
    less_than = -1;
    if( symmetry != 1) {
      j = 0;
      while( j < newDs.nStructs && unique) {
        
        if( comparePairs( (ds->validStructs)[i].theStruct, newDs.validStructs[j].theStruct, 
                         ds->seqlength, unitLength, symmetry) == 0) {
          unique = FALSE;
          if(compareDnaStructsOutput(&((ds->validStructs)[i]),&((newDs.validStructs)[j])) < 0) {
            less_than = j;
          } 
        }
        j++;
      }
    }
    if( unique) {
      newDs.validStructs[ newDs.nStructs].theStruct = (int*) malloc( sizeof( int)*(newDs.seqlength));
      for( j = 0; j < newDs.seqlength; j++) newDs.validStructs[ newDs.nStructs].theStruct[j] =
        ds->validStructs[i].theStruct[j];
      newDs.validStructs[ newDs.nStructs].error = (ds->validStructs)[i].error;
      newDs.validStructs[ newDs.nStructs].correctedEnergy = 
        (ds->validStructs)[i].correctedEnergy;
      newDs.nStructs++;
    } 
    if(less_than >= 0) {
      for( j = 0 ; j < newDs.seqlength; j++) {
        newDs.validStructs[less_than].theStruct[j] = ds->validStructs[i].theStruct[j];
      }
    }
  }

  qsort(newDs.validStructs, newDs.nStructs,sizeof(oneDnaStruct),&compareDnaStructsOutput);

  copyDnaStructures( ds, &newDs);
  clearDnaStructures( &newDs);
}

/* ********* */
int comparePairs( const int *p1, const int *p2, const int seqlength, 
                  const int unitLength, const int possibleSymmetry) {

  int i, j;
  int isSame;
  int symPos;
  int symPair;

  for( i = 1; i < possibleSymmetry; i++) {
    isSame = TRUE;
    j = 0;
    while( j < seqlength && isSame) {
      symPos = (j + i*unitLength) % seqlength;
      if( p1[ j] == -1)
        symPair = -1;
      else {
        symPair = (p1[ j] + i*unitLength) % seqlength;
      }

      //printf("%d %d %d %d %d %d\n", i, j, symPos, p1[j], symPair);
      
      if( p2[ symPos] != symPair)
        isSame = FALSE;
      j++;
    }

    if( isSame) {
      //printf("repeat\n");
      return 0;
    }
  }
  //printf("new\n");
  return 1;
}

/***************** */
void getStructureFromParens( char *line, int *pairs, int seqlength) {
  int i, j;
  int braces[seqlength];
  int leftParenIndex;
  
  char pairSymbols[] = { '(', ')', '{','}', '[', ']', '<', '>' };
  int type = 0;
  int nTypes = 4;
  
  for( i = 0; i <= seqlength-1; i++) {
    pairs[i] = -1;
  }

  for( type = 0; type < nTypes; type++) {
    leftParenIndex = 0;		
    for( i = 0; i <= seqlength-1; i++) {
      braces[i] = -5;
    }

    i = 0;
    j = 0;

    while( i <= seqlength - 1) {
      if( leftParenIndex < 0 || leftParenIndex >= seqlength) {
        fprintf(stderr, "Too many %c, not enough %c!\n", pairSymbols[2*type+1], pairSymbols[2*type]);
        exit(1);
      }
      if( line[j] == pairSymbols[ 2*type]) {
        braces[ leftParenIndex++] = i;
      } 
      else if( line[j] == pairSymbols[ 2*type+1]) {
        pairs[ braces[ --leftParenIndex]] = i;
        pairs[ i] = braces[ leftParenIndex];
      }
      j++;
      if( line[j] != '+') {
        i++;
      }
    }
  }
}

/******** */
void PrintStructure( char *thefold, const int *thepairs,  int **etaN, int seqlength, char *filename) {
  /*
  This prints the structure of the fold using a '.' for 
  unpaired bases, and ( ), { }, [ ], < > for pairs. 

  The file specified by file name is appended.

  Initially, thefold[i] = '.' for all i

  If this ever becomes the slow step, it can be optimized to run faster
  */

  int i,j, pos;
  char pairSymbols[] = { '(', ')', '{','}', '[', ']', '<', '>' };
  int type = 0;
  int nTypes = 4;
  int nStrands = etaN[ EtaNIndex( 0.5, seqlength-0.5, seqlength)][0] + 1;
  int **pairlist; // Each row is i,j pair
  int npairs; // number of pairs in structure
  FILE *fp; // The output file

  char *parensString;
  parensString = ( char*) malloc( (seqlength+1)*sizeof( char) );
  int lastL, lastR;

  if ((fp = fopen(filename,"a")) == NULL){
    fprintf(stderr, "Error opening file %s!\n",filename);
    exit(1);
  }

  // Allocate memory for pairlist (this is more than we need, but be safe)
  pairlist = (int **) malloc(seqlength * sizeof(int *));
  for (i = 0; i < seqlength; i++) {
    pairlist[i] = (int *) malloc(2 * sizeof(int));
  }

  // Create pairlist from thepairs
  npairs = 0;
  for( j = 0; j < seqlength; j++) {
    if(thepairs[j] > j) {
      pairlist[npairs][0] = j;
      pairlist[npairs++][1] = thepairs[j];
    }
  }

  // Creat dot-paren structure
  for( i = 0; i < seqlength+1; i++) {
    parensString[i] = '.';
  }

  //offSet = 0;
  lastL = -1; 
  lastR = seqlength;
  for( i = 0; i < seqlength; i++) {
    if( thepairs[i] != -1 && thepairs[i] > i) {
      if( i > lastR || thepairs[i] > lastR) {
        for( j = 0; j < i; j++) {
          if( thepairs[j] > i && thepairs[j] < thepairs[i]) {
            type = (type + 1) % nTypes;
            break;
          }
        }
      }

      parensString[i] = pairSymbols[ 2*type];
      parensString[ thepairs[i]] = pairSymbols[2*type + 1];
      lastL = i;
      lastR = thepairs[i];
    }
  }

  for( i = 0; i < seqlength+nStrands-1; i++) {
    thefold[i] = '.';
  }

  pos = 0;
  for( i = 0; i < seqlength; i++) {
    thefold[ pos++] = parensString[ i];
    if( etaN[ EtaNIndex_same(i+0.5, seqlength)][0] == 1) {
      thefold[ pos++] = '+';
    }
  }

  // Print the structure
  for( i = 0; i < seqlength+nStrands-1; i++) {
    fprintf(fp,"%c", thefold[ i]);
  }

  fprintf(fp,"\n");

  for (j = 0; j < npairs; j++) {
    fprintf(fp,"%d\t%d\n",pairlist[j][0]+1,pairlist[j][1]+1);
  }

  fclose(fp);
  free( parensString); parensString = NULL;

  for (i = 0; i < seqlength; i++) {
    free(pairlist[i]);
  }
  free(pairlist);

}

/* *************** */
void PrintDnaStructures( const dnaStructures *ds, int **etaN, const int *nicks, int symmetry, char *filename) {
  int i,j;
  int nStrands = etaN[ EtaNIndex( 0.5, ds->seqlength-0.5, ds->seqlength)][0]+1;
  char *foldParens = (char*) malloc( (ds->seqlength + nStrands) * sizeof(char));
  int nPercent = 40; // Number of percent signs in a line with all comments
  FILE *fp; // File to write output to
  
  if ((fp = fopen(filename,"a")) == NULL){
    fprintf(stderr, "Error opening file %s!\n",filename);
    exit(1);
  }
  fclose(fp);

  for( i = 0; i < ds->nStructs; i++) {
    // Print a comment line for separation
    fp = fopen(filename,"a");
    fprintf(fp,"%% ");
    for (j = 0; j < nPercent; j++) {
      fprintf(fp,"%%");
    }

    fprintf(fp," %%\n");

    fprintf(fp,"%d\n", ds->seqlength);
    if(!NUPACK_VALIDATE) {
      fprintf(fp,"%.3Lf\n", (long double) ds->validStructs[i].correctedEnergy);
    } else {
      fprintf(fp,"%.14Le\n", (long double) ds->validStructs[i].correctedEnergy);
    }
    fclose(fp);

    PrintStructure( foldParens, (ds->validStructs)[i].theStruct, etaN, ds->seqlength, filename);

    // Print a comment line for separation
    fp = fopen(filename,"a");
    fprintf(fp,"%% ");
    for (j = 0; j < nPercent; j++) {
      fprintf(fp,"%%");
    }

    fprintf(fp," %%\n\n");
    fclose(fp);
  }
  free( foldParens);
}

/* ******* */
void PrintS( const dnaStructures *ds) {
  //print all the structures currently held in ds, ignoring strand breaks
  int i, j;

  for( i = 0; i < ds->nStructs; i++) {
    for( j = 0; j < ds->seqlength; j++) {
      if( (ds->validStructs)[i].theStruct[j] > j) {
        printf("(");
      }
      else if(  (ds->validStructs)[i].theStruct[j] == -1 ) {
        printf(".");
      }
      else printf( ")");
    }
    printf(" %.2f  %.2f\n", (float) (ds->validStructs)[i].error, (float) (ds->validStructs)[i].correctedEnergy);
  }
}


/* ************* */
void copyDnaStructures( dnaStructures *target, const dnaStructures *source) {
  //clear target and make a copy of source
  int i, j;
  
  clearDnaStructures( target);
  target->nAlloc = source->nAlloc;
  target->nStructs = source->nStructs;
  target->minError = source->minError;
  target->seqlength = source->seqlength;
  
  target->validStructs = (oneDnaStruct *) malloc( sizeof( oneDnaStruct)*(target->nAlloc));
  for( i = 0; i < target->nStructs; i++) {
    (target->validStructs)[i].theStruct = (int *) malloc( (target->seqlength)*sizeof( int));
    for( j = 0; j < target->seqlength; j++) {
      (target->validStructs)[i].theStruct[j] = (source->validStructs)[i].theStruct[j]; 
    }

    (target->validStructs)[i].error = (source->validStructs)[i].error;
    (target->validStructs)[i].correctedEnergy = (source->validStructs)[i].correctedEnergy;
  }
}

/* ************ */
void clearDnaStructures( dnaStructures *ds) {
  int i;

  for( i = 0; i < ds->nStructs; i++) {
    free( (ds->validStructs)[i].theStruct);
    (ds->validStructs)[i].theStruct = NULL;
  }

  free( ds->validStructs);
  ds->validStructs = NULL;
  ds->nStructs = 0;
  ds->nAlloc = 0;
  ds->minError = NAD_INFINITY;
  //ds->seqlength = 0;
}

/* ************* */
void addDnaStructures( dnaStructures *target, 
                      const dnaStructures *source, 
                      DBL_TYPE newErr, 
                      DBL_TYPE mfeEpsilon, 
                      int add_min) {

  int i, j;
  int maxStructs = target->nStructs + source->nStructs;
  DBL_TYPE addedError;
  int si;
  oneDnaStruct * temporary;

  if (add_min) {
    if (target->nAlloc == 0) {
      target->validStructs = (oneDnaStruct*) malloc( 1 * sizeof(oneDnaStruct));
      target->validStructs[0].theStruct = (int*) malloc((target->seqlength)*sizeof(int));
      target->nAlloc = 1;
    } 

    for (i = 0; i < source->nStructs; i++) {
      addedError = newErr + (source->validStructs)[i].error;
      if ( addedError < mfeEpsilon + ENERGY_TOLERANCE) {
        if (target->nStructs == 0 || addedError < target->minError) {
          target->validStructs[0].error = addedError;
          target->validStructs[0].correctedEnergy = 0;

          for( j = 0; j < target->seqlength; j++) {
            (target->validStructs)[0].theStruct[j] = (source->validStructs)[i].theStruct[j]; 
          }
          target->minError = addedError;
          target->nStructs = 1;
        }
      }
    }
  } else {

    if( target->nAlloc == 0) {
      target->validStructs = (oneDnaStruct*) malloc( maxStructs*sizeof( oneDnaStruct));	
      target->nAlloc = maxStructs;
    }
    else if( target->nAlloc < maxStructs) {
      temporary = (oneDnaStruct*) calloc(maxStructs,sizeof(oneDnaStruct));
      memcpy(temporary,target->validStructs,sizeof(oneDnaStruct)*target->nAlloc);
      target->validStructs = temporary;
      target->nAlloc = maxStructs;
    }

    for( i = 0; i < source->nStructs; i++) {
      addedError = newErr + (source->validStructs)[i].error;
      if( addedError < mfeEpsilon + ENERGY_TOLERANCE) {
        si = (target->nStructs)++;
        (target->validStructs)[ si ].error = addedError;
        (target->validStructs)[ si ].correctedEnergy = 0;
        (target->validStructs)[ si ].theStruct = (int *) malloc( (target->seqlength)*sizeof( int));
        
        for( j = 0; j < target->seqlength; j++) {
          (target->validStructs)[ si].theStruct[j] = (source->validStructs)[i].theStruct[j]; 
        }
        if( addedError < target->minError) target->minError = addedError;
      }
    }
  }
}


/* ******************************************** */
int WithinEps( DBL_TYPE x, DBL_TYPE y, DBL_TYPE tol1, DBL_TYPE tol2) {
  // Check if y is in (x-tol1, x+tol2);

  if( y > x + tol2 || y < x-tol1) {
    return FALSE;
  }
  else {
    return TRUE;
  }
}

void print_pf_mat(DBL_TYPE *mat, int len, char *name, FILE *stream) {
  fprintf(stream, "Matrix name: %s\n", name);
  int i, j;
  double val;
  for (i = 0; i < len; i++) {
    for (j = 0; j < len; j++) {
      while (j < i-1) {
        fprintf(stream, "%7.1f ", 0.0);
        j++;
      }
      val = (double) mat[pf_index(i,j,len)];
      if (val < 99999.9) {
        fprintf(stream, "%7.1f ", val);
      } else {
        fprintf(stream, "%2.1e ", val);
      }
    }
    fprintf(stream, "\n");
  }
}


