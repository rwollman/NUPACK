#ifdef NUPACK_SAMPLE
/*  
    sample.c is part of the NUPACK software suite
    Copyright (c) 2007 Caltech. All rights reserved.
    Coded by: Robert Dirks, 6/2006 and Justin Bois 1/2007 
    Brian Wolfe 3/2009
    
    This program will print out N structure samples based 
    on a nucleic acid input sequence
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
    
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
  int complexity = 3;
  int length, tmpLength;
  char inputFile[ MAXLINE];
  int vs;
  char sampleFile[MAXLINE];
  //permAvgPairs stores the expected value of each 
  //class of base pairs, grouped by permutation or complex, respectively

  int inputFileSpecified;
  FILE *F_sample = NULL; // ppairs file
  int index;

  strcpy( inputFile, "");
  nupack_sample = 1;
  nupack_num_samples = 10;
  struct timeval rand_time;

  gettimeofday(&rand_time,0);
  nupack_random_seed = (rand_time.tv_sec)*1000000 + rand_time.tv_usec;

  inputFileSpecified = ReadCommandLineNPK( argc, argv, inputFile);

  if(NupackShowHelp) {
    printf("Usage: sample [OPTIONS] PREFIX\n");
    printf("Randomly sample unpseudoknotted structures from the equilibrium distribution\n");
    printf("Example: sample -multi -T 25 -material dna -samples 100 example\n");
    PrintNupackThermoHelp();
    PrintNupackUtilitiesHelp();
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

  strncpy(sampleFile,inputFile,strlen(inputFile)-3);
  sampleFile[strlen(inputFile)-3] = '\0';
  strcat(sampleFile,".sample");

  header( argc, argv, "sample", sampleFile);
  printInputs( argc, argv, seqChar, vs, NULL, NULL,sampleFile);


  tmpLength = length = strlen( seqChar);
  convertSeq(seqChar, seqNum, tmpLength);
  int ns1,ns2;
  getSequenceLength(seqChar, &ns1);
  getSequenceLengthInt(seqNum, &ns2);

  init_genrand(nupack_random_seed);

  pairPr = NULL;
  if (complexity != 3) {
    printf("Sampling supported only for complexity = 3. Exiting\n");
    exit(1);
  }

  nupack_sample_list = (char **)calloc(nupack_num_samples, sizeof(char *));
  printf("Number of Samples = %i\n",nupack_num_samples);
  for(index = 0 ; index < nupack_num_samples ; index++) {
    nupack_sample_list[index] = (char *) calloc(tmpLength+1,sizeof(char));
  }

  printf("Started Calculation\n");
  pf = pfuncFull(seqNum, complexity, DNARNACOUNT, DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, 0,
      SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_HELIX_FOR_SALT_CORRECTION);
  printf("Finished Calculation\n");


  

  
  if ((F_sample = fopen(sampleFile,"a")) == NULL) {
    printf("Error opening file %s!\n",sampleFile);
    exit(1);
  }

  // Print the free energy to the output file
  fprintf(F_sample,"%s Free energy: %.8Le kcal/mol\n",
          COMMENT_STRING,-kB*TEMP_K*logl(pf));
  fprintf(F_sample,"%s Number of Samples: %i\n",COMMENT_STRING,nupack_num_samples);

  // Put newline for stylistic reasons
  fprintf(F_sample,"\n");


  for(index = 0 ; index < nupack_num_samples ; index++) {
    fprintf(F_sample, "%s\n",nupack_sample_list[index]);
    free(nupack_sample_list[index]);
    nupack_sample_list[index] = NULL;
  }
  free(nupack_sample_list);

#ifdef GC_DEBUG
  CHECK_LEAKS();
#endif
  return 0;
}
/* ****** */



#endif // NUPACK_SAMPLE
