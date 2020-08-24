
#include "design_engine.h"
#include <shared.h>
#include <thermo.h>

int main(int argc, char **argv) {

int myrank = 0;
  
  SetExecutionPath(argc, argv);
  if (myrank == 0) {
  
  char *inputPrefix = (char*)calloc(MAX_FILENAME_SIZE,sizeof(char));
  char *psFile = (char*)calloc(MAX_FILENAME_SIZE,sizeof(char));
  
  char *foldFile;
  
  int loadSeeds = 0;
  int output_init = 0;
  int output_seed = 0;
  int output_ppairs = 0;
  int output_json = 0;
  readCommandLine(argc,argv, inputPrefix, psFile, &initMode, &bypassDesign, 
      &bypassHierarchy, &bypassGuidance, &designMode, &loadSeeds, &mReopt, 
      &mLeafopt, &nRatio, &loadInit, &hMin, &ppairsCutoff, &quickFlag, 
      &output_init,&output_seed, &output_ppairs, &output_json);

  initEngine(psFile);
  free(psFile);
  
  char *fprefix = inputPrefix;

  char *seedFilename = (char *)calloc(MAX_FILENAME_SIZE, sizeof(char));
  sprintf(seedFilename, "%s.seed", fprefix);
  
  #ifdef EARLY_TERMINATE
    timerFilename = (char *)calloc(MAX_FILENAME_SIZE, sizeof(char));
    sprintf(timerFilename, "%s.timer", fprefix);
  #endif
  
  FILE *seedFile = NULL;
  if (loadSeeds) {
    seedFile = fopen(seedFilename, "r");
  }
  else if(output_seed) {
    seedFile = fopen(seedFilename, "w");
  }
  
  free(seedFilename);
    
  foldFile = (char *)calloc(MAX_FILENAME_SIZE,sizeof(char));
  
  sprintf(foldFile, "%s", inputPrefix);
  char *initFilename = (char *)calloc(MAX_FILENAME_SIZE, sizeof(char));
  sprintf(initFilename, "%s.init", foldFile);
  
  FILE *initFile = NULL;
  
  if (loadInit) {
    initFile = fopen(initFilename, "r");    
  }
  else if(output_init) {
    initFile = fopen(initFilename, "w");    
  }
  
  free(initFilename);

  // char *seqFile = (char *)calloc(MAX_FILENAME_SIZE, sizeof(char));
  // sprintf(seqFile, "%s.seq", foldFile);
  // FILE *seqOut = fopen(seqFile, "w");
  // free(seqFile);  
  
  #ifdef DEBUG_TEST
    unitTester();
  #endif
    
  int stringSum = 0;
  int z;
  for (z = 0; z < strlen(inputPrefix); z++) {
    stringSum += (int)inputPrefix[z];
  }
  
  unsigned rand_seed;
  if (loadSeeds) {
    fscanf(seedFile, "%u\n", &rand_seed);
  }
  else {
    /* get microsecond resolution time for random seed */
    struct timeval rand_time;
    gettimeofday(&rand_time, 0);
    rand_seed = (unsigned)(rand_time.tv_sec) + stringSum;
    if(output_seed) {
      fprintf(seedFile, "%u\n",rand_seed);
      fflush(seedFile);
    }
  }
  init_genrand(rand_seed);
  
  secStruc *myStruc = (secStruc *)calloc(1, sizeof(secStruc));
  time_t start_time_t;
  struct tm * start_s_tm;
  start_time_t = time(NULL);
  start_s_tm = localtime(&start_time_t);
  
  struct timeval start_timeval;
  struct timeval end_timeval;
  struct timeval elapsed_timeval;
  
  #ifdef EARLY_TERMINATE
    global_timeval = &start_timeval;
  #endif
  
  DBL_TYPE elapsed;
  
  gettimeofday(&start_timeval, 0);
  
  parseStructure(myStruc, foldFile);
  #ifdef FINDBUG
  return 0;
  #endif

  #ifdef DEBUG
  char *treeFile = (char *)calloc(MAX_FILENAME_SIZE, sizeof(char));
  sprintf(treeFile, "%s.dot", foldFile);
  FILE *treeOut = fopen(treeFile, "w");
  free(treeFile);  
  fprintf(treeOut, "digraph decompositionTree {\n");
  fprintf(treeOut,"ordering = out\n");
  printTree(myStruc, 0, treeOut);
  fprintf(treeOut, "}\n");
  fclose(treeOut);
  #endif

  if (myStruc->length > MAX_SEQ_LENGTH || myStruc->length < MIN_SEQ_LENGTH) {
    fprintf(stderr, "STRUCTURE %s NOT RIGHT SIZE!\n", foldFile);
    exit(1);
  }
  
  if (loadInit) {
    loadStringFromInitFile(initFile, myStruc->length, myStruc->seq);
    if (violatesSkeleton(myStruc->length, myStruc->seq, myStruc->skel)) {
      
      fprintf(stderr,"Given initial sequence violates constraints:\n");
      char *tempSeq = (char *)calloc(myStruc->length + 1, sizeof(char));
      convertIntsToBases(myStruc->seq, tempSeq, myStruc->length);
      fprintf(stderr,"GIVEN     : %s\n", tempSeq);
      convertIntsToBases(myStruc->skel, tempSeq, myStruc->length);
      fprintf(stderr,"CONSTRAINT: %s\n", tempSeq);
      exit(1);
    }
    fixPatternViolations(myStruc);  

  }
  else {
    initSeq(myStruc);
    char *initSeq = (char *)calloc(myStruc->length + 1, sizeof(char));
    convertIntsToBases(myStruc->seq, initSeq, myStruc->length);
    if(output_init) {
      fprintf(initFile, "%s\n",initSeq);
    }
    free(initSeq);
  }
  
  if (!bypassDesign) {
    myStruc->foldFile = foldFile;
    designSeq(myStruc, NULL);

  }
  else {
    gettimeofday(&end_timeval, 0);

  }
  
  DBL_TYPE nVal = (DBL_TYPE)myStruc->length;
  DBL_TYPE mfeVal = (DBL_TYPE)myStruc->length;
  DBL_TYPE pVal = 0.0;
  
  if (myStruc->checked) {
    /* Value was calculated */
    if (designMode == N_OPTIMIZATION) {
      nVal = myStruc->optVal;
      #ifdef DEBUG
        printf("%Lf:%Lf\n",myStruc->optVal, myStruc->nVal);
        assert(floatsEqual(myStruc->optVal,myStruc->nVal));
      #endif
    }
    else if (designMode == MFE_OPTIMIZATION) {
      mfeVal = myStruc->optVal;
      #ifdef DEBUG
        assert(floatsEqual(myStruc->optVal, myStruc->bestMFEDifference));
      #endif
    }
    else {
      #ifdef DEBUG
        assert(designMode == P_OPTIMIZATION);
      #endif
      pVal = myStruc->optVal;
      #ifdef DEBUG
        assert(floatsEqual(myStruc->optVal, myStruc->bestPVal));
      #endif
    }
    gettimeofday(&end_timeval, 0);

  }
  
  else {
    /* We skipped calculations on the highest level */
    if (designMode == N_OPTIMIZATION) {
      computeNSstar(myStruc);
      nVal = myStruc->optVal;
    }
    else if (designMode == MFE_OPTIMIZATION) {
      computeMFE(myStruc);
      //mfeVal = myStruc->optVal;
    }
    else {
      #ifdef DEBUG
        assert(designMode == P_OPTIMIZATION);
      #endif
      computeP(myStruc);
      //pVal = myStruc->optVal;
    }
    
    if (myStruc->optVal < myStruc->bestOptVal) {
      printf("SETTING: %Lf,%Lf\n", myStruc->optVal, myStruc->bestOptVal);
      setNewBestValue(myStruc);
    }
    
    
    //if (mReopt == 0) {
    gettimeofday(&end_timeval, 0);
  
    //}
  }    

  revertToOldValue(myStruc);
  
  if (designMode == N_OPTIMIZATION) {
    nVal = myStruc->optVal;
  }
  else if (designMode == MFE_OPTIMIZATION) {
    mfeVal = myStruc->optVal;
  }
  else {
    #ifdef DEBUG
      assert(designMode == P_OPTIMIZATION);
    #endif
    pVal = myStruc->optVal;
  }
  
  
  if (timeDifference (&elapsed_timeval, &end_timeval, &start_timeval)) {
    fprintf(stderr,"Timing function returned negative number\n");
    exit(1);
  }
  
  elapsed =  (DBL_TYPE)elapsed_timeval.tv_sec + 1e-6 * (DBL_TYPE)elapsed_timeval.tv_usec;
  
  printf("ELAPSED TIME: %Lf\n", elapsed);

  #ifdef DEBUG  
    assert(seqCmp(myStruc->seq, myStruc->bestSeq, myStruc->length) == 0);
    
    if (designMode == N_OPTIMIZATION) {
      computeNSstar(myStruc);
      printf("NVAL: %Lf,%Lf\n", (long double)nVal, (long double)myStruc->optVal);
      assert(floatsEqual(nVal, myStruc->optVal));
      assert(floatsEqual(nVal, myStruc->bestOptVal));
  
      
    }
    else if (designMode == MFE_OPTIMIZATION) {
      assert(floatsEqual(myStruc->bestOptVal, myStruc->bestMFEDifference));
      computeMFE(myStruc);
      assert(floatsEqual(mfeVal, myStruc->optVal));
      assert(floatsEqual(mfeVal, myStruc->bestMFEDifference));
    }
    else {
      assert(designMode == P_OPTIMIZATION);
      assert(floatsEqual(myStruc->bestOptVal,myStruc->bestPVal));
      computeP(myStruc);
      assert(floatsEqual(pVal,myStruc->optVal));
      assert(floatsEqual(pVal,myStruc->bestPVal));
    }
    
    assertInSync(myStruc);
  #endif

  if (!quickFlag) {
    if (designMode == N_OPTIMIZATION) {
      computeMFE(myStruc);
      computeP(myStruc);
    }
    else if (designMode == MFE_OPTIMIZATION) {
      computeNSstar(myStruc);
      computeP(myStruc);
    }
    else {
      computeMFE(myStruc);
      computeNSstar(myStruc);
    }
  }
  
  setAllBestValues(myStruc);
  nVal = myStruc->nVal;
  mfeVal = myStruc->bestMFEDifference;
  pVal = myStruc->bestPVal;
  
  char *tempSeq = (char *)calloc(myStruc->length + 1, sizeof(char));
  int *tempSeqExpInts = (int *)calloc(myStruc->expandedLength + 1, sizeof(int));
  tempSeqExpInts[myStruc->expandedLength] = -1;
  
  char *tempSeqExp = (char *)calloc(myStruc->expandedLength + 1, sizeof(char));

  createExpandedSeq(tempSeqExpInts, myStruc->seq, myStruc->length, myStruc->numStrands, myStruc->strandBreaks);

  
  DBL_TYPE freeEnergy = (-1.0 * (TEMP_K * kB)) * (logl(myStruc->pfVal));
  
  
  DBL_TYPE freeEnergyTarget = naEnergyPairsOrParensFullWithSym( myStruc->tempPairing, NULL, tempSeqExpInts, 
         DNARNACOUNT, DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, 1, SODIUM_CONC, 
         MAGNESIUM_CONC, USE_LONG_SALT);
  DBL_TYPE *base_dist_t = (DBL_TYPE *)calloc(NUM_BASES + 1, sizeof(DBL_TYPE));
  
  calculateBaseDist(base_dist_t, myStruc->seq, myStruc->length);

  DBL_TYPE CGTot = base_dist_t[BASE_C]+base_dist_t[BASE_G];  
  
  char *structure = (char *)calloc(myStruc->length + 1, sizeof(char));
  
  getStructure( myStruc->length, myStruc->bestPairs, structure); 
   
  char *outFile = (char *)calloc(MAX_FILENAME_SIZE, sizeof(char));
  FILE * out;
  sprintf(outFile, "%s.summary", foldFile);
  designHeader(argc, argv, outFile,start_s_tm,myStruc->expFolding);
  out = fopen(outFile,"a");
     
  #ifdef LEAFCORRECTION
  DBL_TYPE leafSum2 = (DBL_TYPE)sumLeafIDs(myStruc);
  #else
  DBL_TYPE leafSum2 = -1.0;
  #endif

  convertIntsToBases(myStruc->seq, tempSeq, myStruc->length);
  convertIntsToBases(tempSeqExpInts, tempSeqExp, myStruc->expandedLength);
  
  fprintf(out, "%s\n",COMMENT_STRING);
  // Timing
  fprintf(out, "%s Elapsed time: %Lf\n",COMMENT_STRING, (long double)elapsed);
  if (designMode == N_OPTIMIZATION) {
    fprintf(out, "%s Root evaluation time: %Lf\n", COMMENT_STRING, myStruc->pfTimeSpent);
  }
  else if (designMode == MFE_OPTIMIZATION) {
    fprintf(out, "%s Root evaluation time: %Lf\n", COMMENT_STRING, myStruc->mfeTimeSpent);
  }
  else {
    fprintf(out, "%s Root evaluation time: %Lf\n", COMMENT_STRING, myStruc->pTimeSpent);
  }

  fprintf(out,"%s\n",COMMENT_STRING);
  // Properties of the sequence
  //fprintf(out, "%s Length: %d\n", COMMENT_STRING, myStruc->length);
  //fprintf(out, "%s Independent bases: %d\n", COMMENT_STRING, myStruc->indyLength);
  fprintf(out, "%s CG fraction: %Lf\n", COMMENT_STRING, (long double)CGTot);
  fprintf(out, "%s Pattern violations: %d\n", COMMENT_STRING, violatesPattern(myStruc->seq, myStruc->skel,myStruc->length*2));

  // Energetics of the sequence
  fprintf(out, "%s Target free energy: %Lf\n",COMMENT_STRING, (long double)freeEnergyTarget);
  fprintf(out, "%s Ensemble free energy: %Lf\n", COMMENT_STRING, (long double)freeEnergy);
  fprintf(out,"%s\n",COMMENT_STRING);

  if (mReopt <= 0) {
    DBL_TYPE leafSum = (DBL_TYPE)calculateLeafSum(myStruc);
    #ifdef DEBUG
      floatsEqual(leafSum,leafSum2);
    #endif
     fprintf(out, "%s Leaf sum: %Lf\n", COMMENT_STRING,leafSum);

  }
  else {
    fprintf(out, "%s Leaf sum: %Lf\n", COMMENT_STRING, leafSum2);
  }

  fprintf(out, "%s MFE defect: %Lf\n",COMMENT_STRING,(long double)mfeVal);
  fprintf(out, "%s Probability: %Lf\n" , COMMENT_STRING,(long double)(pVal));
  fprintf(out, "%s Ensemble defect: %Lf\n", COMMENT_STRING, (long double)(nVal));
  fprintf(out, "%s\n", tempSeqExp);

  // if (!verifyStructure(structure, myStruc->tempFolding)) {
  //   fprintf(out,"0,");
  // } 
  // else {
  //   fprintf(out,"1,");
  // }
  
   //printf("Leaf Sum = %Lf\n", (long double)leafSum);
  
  // fprintf(seqOut, "%s\n", tempSeqExp);
  printf("Length : %d\n",myStruc->length);
  printf("N(S*) : %Lf\n",(long double)nVal);
  printf("P : %Lf\n",(long double)(pVal));
  printf("MFE : %Lf\n",(long double)mfeVal);
  printf("Sequence : %s\n", tempSeqExp);
  
  if (!quickFlag) {
    char * baseProbabilityFilename;
    FILE * bpf = NULL;
    char * probMatrixFilename;
    FILE * pmf = NULL;

    if(output_json) {
      baseProbabilityFilename = (char*)calloc(MAX_FILENAME_SIZE,sizeof(char));
      sprintf(baseProbabilityFilename, "%s.json", foldFile);
      bpf = fopen(baseProbabilityFilename, "w");
      free(baseProbabilityFilename);
    }
  
    if(output_ppairs) {
      probMatrixFilename = (char*)calloc(MAX_FILENAME_SIZE,sizeof(char));
      sprintf(probMatrixFilename, "%s.ppairs", foldFile);

      header(argc,argv,"design",probMatrixFilename);
      printInputs(argc,argv,tempSeq,1,NULL,NULL,probMatrixFilename);

      pmf = fopen(probMatrixFilename, "a");
      free(probMatrixFilename);
      fprintf(pmf, "%s Free energy: %Le\n",COMMENT_STRING,(long double) freeEnergy);
    
      fprintf(pmf, "\n");
      fprintf(pmf, "%d\n",myStruc->length);
      
    }
    
    
    int p;
    
    if(output_json) {
     
      struct tm * loctime;
      loctime = localtime(&start_time_t);

      fprintf(bpf, "{\n");
      fprintf(bpf, "\t\"NUPACK version\":\"%s\",\n", NUPACK_VERSION);
      fprintf(bpf, "\t\"program\":\"design\",\n");
      char * time_src = asctime(loctime);
      int time_len = strlen(time_src);
      char * time_str = (char *)calloc(time_len,sizeof(char));
      if(time_str != NULL) {
        int strpos;
        for(strpos = 0 ; strpos < time_len && time_src[strpos] != '\n'; strpos++) {
          time_str[strpos] = time_src[strpos];
        }
        fprintf(bpf, "\t\"start time\":\"%s\",\n",time_str);
      }
      fprintf(bpf, "\t\"command\":\"");
      int argi;
      for(argi = 0 ; argi < argc; argi++) {
        fprintf(bpf,"%s ", argv[argi]);
      }
      fprintf(bpf, "\",\n");
      fprintf(bpf, "\t\"parameters\":\"");
      if( DNARNACOUNT == DNA) {
        fprintf(bpf, "DNA, (Mfold 2.3)\",\n");
      }
      else if(  DNARNACOUNT == RNA) {
        fprintf(bpf, "RNA, (Mfold 2.3)\",\n");
      }
      else if( DNARNACOUNT == RNA37) {
        fprintf(bpf, "RNA, (Mfold 3.0)\",\n");
      }
      else if( DNARNACOUNT == USE_SPECIFIED_PARAMETERS_FILE) {
        fprintf(bpf, "Custom, (%s)\",\n", PARAM_FILE);
      }

      fprintf(bpf, "\t\"dangles setting\": %d,\n", DANGLETYPE);
      fprintf(bpf, "\t\"temperature (C)\": %.1f,\n",
              (float) (TEMP_K - ZERO_C_IN_KELVIN) );

      // Say what salt concentrations were used in the calculation
      fprintf(bpf,"\t\"sodium concentration\": %.4f,\n", (float) SODIUM_CONC);
      fprintf(bpf,"\t\"magnesium concentration\": %.4f,\n",(float) MAGNESIUM_CONC);


      fprintf(bpf, "\t\"structure\":\"%s\",\n",myStruc->expFolding);
      fprintf(bpf, "\t\"sequence\":\"%s\",\n",tempSeqExp);
      fprintf(bpf, "\t\"energypreamble\":\"Free energy of secondary structure:\",\n");
      fprintf(bpf, "\t\"target free energy\": %Lf,\n", (long double)freeEnergyTarget);
      fprintf(bpf, "\t\"ensemble free energy\": %Lf,\n", (long double)freeEnergy);
      fprintf(bpf, "\t\"elapsed time\": %Lf,\n", (long double)elapsed);
      if(designMode == N_OPTIMIZATION) {
        fprintf(bpf, "\t\"root evaluation time\": %Lf,\n", (long double)myStruc->pfTimeSpent);
      } else if(designMode == MFE_OPTIMIZATION) {
        fprintf(bpf, "\t\"root evaluation time\": %Lf,\n", (long double)myStruc->mfeTimeSpent);
      } else {
        fprintf(bpf, "\t\"root evaluation time\": %Lf,\n", (long double)myStruc->pTimeSpent);
      }
      fprintf(bpf, "\t\"length\": %d,\n", myStruc->length);
      fprintf(bpf, "\t\"independent bases\": %d,\n", myStruc->indyLength);
      fprintf(bpf, "\t\"CG fraction\": %Lf,\n", (long double)CGTot);
      fprintf(bpf, "\t\"pattern violations\": %d,\n", violatesPattern(myStruc->seq, myStruc->skel, myStruc->length * 2));
      if(mReopt <= 0) {
        DBL_TYPE leafSum = (DBL_TYPE)calculateLeafSum(myStruc);
        fprintf(bpf, "\t\"leaf sum\": %Lf,\n", (long double)leafSum);
      } else {
        fprintf(bpf, "\t\"leaf sum\": %Lf,\n", (long double)leafSum2);
      }
      fprintf(bpf, "\t\"mfe defect\": %Lf,\n", (long double)mfeVal);
      fprintf(bpf, "\t\"probability\": %Lf,\n", (long double)pVal);
      fprintf(bpf, "\t\"ensemble defect\": %Lf,\n", (long double)nVal);
      fprintf(bpf, "\t\"nucleotide probabilities\":[\n\t\t");
    }
    for (p = 0; p < myStruc->length; p++) {
      if (myStruc->tempPairing[p] >= 0) {
        if(output_json) {
          fprintf(bpf, "%Lf", myStruc->myPairPr[p*(myStruc->length+1) +myStruc->tempPairing[p]]);
        }
      }
      else {
        if(output_json) {
          fprintf(bpf, "%Lf", myStruc->myPairPr[p*(myStruc->length+1) + myStruc->length ]);
        }
      }
      if (p < myStruc->length - 1) {
        if(output_json) {
          fprintf(bpf, ",\n\t\t");
        }
      }
      int q;
      if(output_ppairs) {
        for (q = p; q < myStruc->length + 1; q++) {
          if (myStruc->myPairPr[p*(myStruc->length+1) +q] >= ppairsCutoff) {
            fprintf(pmf, "%d %d %Lf\n", p+1, q+1, myStruc->myPairPr[p*(myStruc->length+1) +q]);
          }
        }
      }
    }
    if(output_json) {
      fprintf(bpf, "\n\t]\n}\n");
      fclose(bpf);
    }
  }
  free(tempSeq);
  free(tempSeqExpInts);
  free(tempSeqExp);


  
  printf("********************\n");
  
  fflush(out);
  // fflush(seqOut);

  free(base_dist_t);
  
  free(structure);
  freeStruc(myStruc);
  free(myStruc);

  
  if(loadInit || output_init) {
    fclose(initFile);
  }
  if(loadSeeds || output_seed) {
    fclose(seedFile);
  }
  fclose(out);
  // fclose(seqOut);
  free(inputPrefix);
  free(foldFile);

  freeEngine();
  }
  else {
  }
  return 0;
}
