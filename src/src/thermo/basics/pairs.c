/*  
    pairs.c is part of the NUPACK software suite
    Copyright (c) 2007 Caltech. All rights reserved.
    Coded by: Robert Dirks, 6/2006 and Justin Bois 1/2007
    
    This program will print out the pair probabilities for every pair
    of bases in a nucleic acid ordered complex.  In addition, for a
    sequence of length N, the "pair" probability between base i and
    base N+1 is the probability that base i is unpaired.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <shared.h>
#include <thermo/core.h>


/* ************************************************ */
extern DBL_TYPE *pairPrPbg;  //for pseudoknots
extern DBL_TYPE *pairPrPb;  //for pseudoknots

extern double CUTOFF;
extern int Multistranded;
extern int perm[MAXSTRANDS];
extern int seqlengthArray[MAXSTRANDS];
extern int nUniqueSequences; // Number of unique sequences entered

int main( int argc, char *argv[] ) {

  char seqChar[MAXSEQLENGTH];
  int seqNum[MAXSEQLENGTH+1];  

  DBL_TYPE pf;
  int stoichiometry[MAXSTRANDS]; // stoichiometry[i] = # of strands of sequence i in complex
  int complexity = 3;
  int length, tmpLength;
  int i, j, q, r;
  char inputFile[ MAXLINE];
  int vs;
  int nStrands = 1; // Number of strands in complex
  int baseNumber = 0; // Number of a base in a strand
  int *baseCode = NULL; // Code for base/strand ID, as in the code complexes
  int **indexCx = NULL;
  int totalStrandsLength = 0;
  char ppairsFile[MAXLINE];
  char epairsFile[MAXLINE];
  //permAvgPairs stores the expected value of each 
  //class of base pairs, grouped by permutation or complex, respectively
  long double *permAvgPairs = NULL; 
  long double ppf;
  long double expectedUnpaired;
  long double **permPr = NULL;
  long double **avgBp = NULL;
  int indQ;
  int pf_qr;
  int strandId, strandPos, strandId2, strandPos2;
  int row, column;

  int inputFileSpecified;
  FILE *F_permPr = NULL; // ppairs file
  FILE *F_permAvg = NULL; // epairs file

  strcpy( inputFile, "");
  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);

  if(NupackShowHelp) {
    printf("Usage: pairs [OPTIONS] PREFIX\n");
    printf("Calculate the pair probabilities of the input sequence.\n");
    printf("Example: pairs -multi -T 25 -material dna example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
    printf("Output control:\n");
    printf(" -cutoff CUTOFFVALUE    only probabilities and expected values\n");
    printf("                        at or above CUTOFFVALUE are saved in the\n");
    printf("                        output file(s)\n");
    exit(1);
  }

  if( !inputFileSpecified ) {
    printf("Enter output file prefix: ");
    scanf("%s", inputFile);
    strcat(inputFile,".in"); // Here, .in is just a placeholder
  }

  if(!inputFileSpecified ||
     !ReadInputFile( inputFile, seqChar, &vs, NULL, NULL, NULL) ) {
       if (inputFileSpecified==0) getUserInput( seqChar, &vs, NULL, NULL);
       else abort();
  }

  strncpy(ppairsFile,inputFile,strlen(inputFile)-3);
  ppairsFile[strlen(inputFile)-3] = '\0';
  strcat(ppairsFile,".ppairs");
  strncpy(epairsFile,inputFile,strlen(inputFile)-3);
  epairsFile[strlen(inputFile)-3] = '\0';
  strcat(epairsFile,".epairs");

  header( argc, argv, "pairs", ppairsFile);
  printInputs( argc, argv, seqChar, vs, NULL, NULL,ppairsFile);

  if( !DO_PSEUDOKNOTS ) {
    complexity = 3;
  }
  else {
    complexity = 5;
  }

  tmpLength = length = strlen( seqChar);
  convertSeq(seqChar, seqNum, tmpLength);
  int ns1,ns2;
  getSequenceLength(seqChar, &ns1);
  getSequenceLengthInt(seqNum, &ns2);

  pairPr = (DBL_TYPE*) calloc( (length+1)*(length+1), sizeof(DBL_TYPE));
  pairPrPbg = (DBL_TYPE*) calloc( (length+1)*(length+1), sizeof(DBL_TYPE));
  pairPrPb = (DBL_TYPE*) calloc( (length+1)*(length+1), sizeof(DBL_TYPE));

  pf = pfuncFullWithSym(seqNum, complexity, DNARNACOUNT, DANGLETYPE, 
      TEMP_K - ZERO_C_IN_KELVIN, 1,vs,
      SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);



  
  if ((F_permPr = fopen(ppairsFile,"a")) == NULL) {
    printf("Error opening file %s!\n",ppairsFile);
    exit(1);
  }

  // Print the free energy to the output file
  if(!NUPACK_VALIDATE) {
    fprintf(F_permPr,"%s Free energy: %.8Le kcal/mol\n",
          COMMENT_STRING,-kB*TEMP_K*logl(pf));
  } else {
    fprintf(F_permPr,"%s Free energy: %.14Le kcal/mol\n",
          COMMENT_STRING,-kB*TEMP_K*logl(pf));
  }

  // Put newline for stylistic reasons
  fprintf(F_permPr,"\n");

  if (Multistranded) {
    header( argc, argv, "pairs", epairsFile);
    printInputs( argc, argv, seqChar, vs, NULL, NULL,epairsFile);

    if ((F_permAvg = fopen(epairsFile,"a")) == NULL) {
      printf("Error opening file %s!\n",epairsFile);
      exit(1);
    }
    fprintf(F_permAvg,"\n");

    // Get base code and total number of bases and number of strands
    //initialize expectation of pairs for complex
    baseCode = (int *) malloc(2*tmpLength * sizeof(int));
    nStrands = 1;
    j = 0;
    for ( i = 0; i < tmpLength; i++) {
      if ( seqChar[i] == '+') {
        length--;
        nStrands++;
        baseNumber = 0;
      }
      else {
        baseCode[2*j] = perm[nStrands-1]-1;
        baseCode[2*j+1] = baseNumber;
        j++;
        baseNumber++;
      }
    }

    // complex index identifiers
    indexCx = (int **) malloc( nUniqueSequences * sizeof( int*));
    totalStrandsLength = 0;
    for( i = 0; i < nUniqueSequences; i++) {
      indexCx[i] = (int *) malloc( seqlengthArray[i]*sizeof( int) );
      for( j = 0; j < seqlengthArray[i]; j++) {
        indexCx[i][j] = totalStrandsLength++;
      }
    }

    // Print strand lengths
    fprintf(F_permPr,"%d\n",length);
    fprintf(F_permAvg,"%d\n",totalStrandsLength);

    // Determine strand multiplicy
    for (i = 0; i < nUniqueSequences; i++) {
      stoichiometry[i] = 0;
    }

    for (i = 0; i < nStrands; i++) {
      stoichiometry[perm[i]-1]++;
    }

    //these will store the average number of pairs for a particular base type
    permAvgPairs = (long double *) malloc (totalStrandsLength*
                                           totalStrandsLength*
                                           sizeof( long double));
    permPr = (long double **) malloc( nUniqueSequences*sizeof( long double*));

    for( j = 0; j < nUniqueSequences; j++) { //calloc initialize to zero
      permPr[j] = (long double *) calloc( seqlengthArray[j], sizeof( long double));
    }

    //initialize avgBp
    avgBp = (long double **) malloc( nUniqueSequences*sizeof( long double*));
    for( j = 0; j < nUniqueSequences; j++) { //calloc initialize to zero
      avgBp[j] = (long double*) calloc( seqlengthArray[j], sizeof( long double));
    }
    
    //initialize expectation of pairs for complex
    for( j = 0; j < totalStrandsLength*totalStrandsLength; j++) {
      permAvgPairs[ j] = 0;
    }

    for( q = 0; q <= length-1; q++) {
      //compute the strand type and relative position for each base in the 
      //concatenation of sequences
      strandId = baseCode[2*q];
      strandPos = baseCode[2*q+1];
      indQ = (1 + length)*q;
      for( r = 0; r <= length-1; r++) {
        strandId2 = baseCode[2*r];
        strandPos2 = baseCode[2*r+1];

        pf_qr = indQ + r;
        ppf = pairPr[pf_qr]*pf;
        avgBp[strandId][strandPos] += ppf;
        permPr[strandId][strandPos] += pairPr[ pf_qr];

        permAvgPairs[ indexCx[strandId][strandPos]*totalStrandsLength +
                     indexCx[strandId2][strandPos2] ] += 
          pairPr[pf_qr];

        if(r > q && pairPr[ pf_qr] >= CUTOFF) {
          if(!NUPACK_VALIDATE) {
            fprintf( F_permPr, "%d\t%d\t%.4Le\n", q+1, r+1, (long double)pairPr[ pf_qr]);
          } else {
            fprintf( F_permPr, "%d\t%d\t%.14Le\n", q+1, r+1, (long double)pairPr[ pf_qr]);
          }
        }
      }
    }

    //add in unpaired probabilities to F_permPr (.ppairs)
    for( q = 0; q <= length-1; q++) {
      r = length;
      pf_qr = (1+length)*q+r;
      if(pairPr[pf_qr] >= CUTOFF) {
        if(!NUPACK_VALIDATE) {
          fprintf( F_permPr, "%d\t%d\t%.4Le\n",q+1, r+1, (long double)pairPr[ pf_qr]);
        } else {
          fprintf( F_permPr, "%d\t%d\t%.14Le\n",q+1, r+1, (long double)pairPr[ pf_qr]);
        }
      }
    }

    //store average number of each type of pair (.epairs)
    for( q = 0; q < totalStrandsLength*totalStrandsLength; q++) {
      row = (q/totalStrandsLength) + 1;
      column = (q % totalStrandsLength) + 1;
      
      if( permAvgPairs[q] >= CUTOFF && row < column ) {
        if(!NUPACK_VALIDATE) {
          fprintf( F_permAvg, "%d\t%d\t%.4Le\n", row, column, permAvgPairs[q]);
        } else {
          fprintf( F_permAvg, "%d\t%d\t%.14Le\n", row, column, permAvgPairs[q]);
        }
      }
    }

    //Also store average unpaired for each base (.epairs)
    for( q = 0; q < nUniqueSequences; q++) {
      for( r = 0; r < seqlengthArray[q]; r++) {
        row = indexCx[q][r] + 1;
        column = totalStrandsLength + 1;
        expectedUnpaired = stoichiometry[q] - permPr[q][r];
        if (expectedUnpaired >= CUTOFF) {
          if(!NUPACK_VALIDATE) {
            fprintf( F_permAvg, "%d\t%d\t%.4Le\n", row, column, expectedUnpaired); 
          } else {
            fprintf( F_permAvg, "%d\t%d\t%.14Le\n", row, column, expectedUnpaired);
          }
        }
      }
    }

    fclose(F_permAvg);
  }
  
  else { // Single strand
    fprintf(F_permPr,"%d\n",length);
    for( i = 0; i < length; i++) {
      for( j = i+1; j < length; j++) {
        if (pairPr[(length+1)*i + j] >= CUTOFF) {
          if( !DO_PSEUDOKNOTS) {
            if(!NUPACK_VALIDATE) {
              fprintf(F_permPr,"%d\t%d\t%.4Le\n",
                    i+1, j+1, (long double) pairPr[ (length+1)*i + j]);
            } else {
              fprintf(F_permPr,"%d\t%d\t%.14Le\n",
                    i+1, j+1, (long double) pairPr[ (length+1)*i + j]);
            }
          } else {
            if(!NUPACK_VALIDATE) {
              fprintf(F_permPr,"%d %d %.4Le %.4Le %.4Le\n", i+1, j+1, 
                    (long double) pairPr[ (length+1)*i + j],
                    (long double) pairPrPb[ (length+1)*i + j],
                    (long double) pairPrPbg[ (length+1)*i + j]);
            } else {
              fprintf(F_permPr,"%d %d %.14Le %.14Le %.14Le\n", i+1, j+1, 
                    (long double) pairPr[ (length+1)*i + j],
                    (long double) pairPrPb[ (length+1)*i + j],
                    (long double) pairPrPbg[ (length+1)*i + j]);
            }
          }
        }
      }
    }
    for(i = 0; i < length ; i++) {
      if (pairPr[(length+1)*i + j] >= CUTOFF) {
        if( !DO_PSEUDOKNOTS) {
          if(!NUPACK_VALIDATE) {
            fprintf(F_permPr,"%d\t%d\t%.4Le\n",
                  i+1, j+1, (long double) pairPr[ (length+1)*i + j]);
          } else {
            fprintf(F_permPr,"%d\t%d\t%.14Le\n",
                  i+1, j+1, (long double) pairPr[ (length+1)*i + j]);
          }
        } else {
          if(!NUPACK_VALIDATE) {
            fprintf(F_permPr,"%d %d %.4Le %.4Le %.4Le\n", i+1, j+1, 
                  (long double) pairPr[ (length+1)*i + j],
                  (long double) pairPrPb[ (length+1)*i + j],
                  (long double) pairPrPbg[ (length+1)*i + j]);
          } else {
            fprintf(F_permPr,"%d %d %.14Le %.14Le %.14Le\n", i+1, j+1, 
                  (long double) pairPr[ (length+1)*i + j],
                  (long double) pairPrPb[ (length+1)*i + j],
                  (long double) pairPrPbg[ (length+1)*i + j]);
          }
        }
      }
    }
  }

  fclose(F_permPr);

  if (Multistranded) {
    free(baseCode);
    free(permAvgPairs);
    for (j = 0; j < nUniqueSequences; j++) {
      free(indexCx[j]);
    }

    free(indexCx);
    for (j = 0; j < nUniqueSequences; j++) {
      free(permPr[j]);
      free(avgBp[j]);
    }

    free(permPr);
    free(avgBp);
  }

  free( pairPr);
  free( pairPrPbg);
  free( pairPrPb);
#ifdef GC_DEBUG
  CHECK_LEAKS();
#endif
  return 0;
}
/* ****** */


