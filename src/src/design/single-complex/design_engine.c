
#include "design_engine.h"
#include "design_pfunc_utils.h"

#ifdef LEAFCORRECTION
int nodeCounter = 0;
#endif


void designHeader(int argc, char ** argv, 
                  char * filename, struct tm * start_s_tm,
                  char * structure) {
  
  FILE * out = fopen(filename,"w");
  
  fprintf(out,"%s NUPACK %s\n", COMMENT_STRING, NUPACK_VERSION);
  fprintf(out,"%s Program: %s\n",COMMENT_STRING,"design");
  fprintf(out,"%s Start time: %s", COMMENT_STRING, asctime( start_s_tm));
  fprintf(out,"%s Command: ", COMMENT_STRING);
  int argi;
  for( argi = 0; argi < argc; argi++) {
    fprintf( out,"%s ", argv[ argi]);
  }
  fprintf( out, "%s\n",COMMENT_STRING);
  fclose(out);
  printInputs(argc,argv,NULL,1,NULL,structure,filename);
  
  out = fopen(filename, "a");
  if(loadInit) {
    fprintf(out, "%s Init mode: LOAD\n",COMMENT_STRING);
  } else {
    switch(initMode) {
      case RANDOM_INIT:
        fprintf(out,"%s Init mode: RANDOM\n",COMMENT_STRING);
        break;
      case AU_INIT:
        fprintf(out,"%s Init mode: AU\n",COMMENT_STRING);
        break;
      case CG_INIT:
        fprintf(out,"%s Init mode: CG\n",COMMENT_STRING);
        break;
      case SSM_INIT:
        fprintf(out,"%s Init mode: SSM\n",COMMENT_STRING);
        break;
      default:
        fprintf(out,"%s Init mode: unknown\n",COMMENT_STRING);
        break;
    }
  }
  fprintf(out,"%s M_UNFAVORABLE: %lf\n",COMMENT_STRING,M_UNFAVORABLE);
  fprintf(out,"%s M_REOPT: %i\n",COMMENT_STRING,mReopt);
  fprintf(out,"%s M_LEAFOPT: %i\n",COMMENT_STRING,mLeafopt);
  fprintf(out,"%s f_STOP: %Lf\n",COMMENT_STRING,nRatio);
  fprintf(out,"%s Minimum split helix: %i\n",COMMENT_STRING,2*hMin);
  fprintf(out,"%s Minimum leaf size: %i\n",COMMENT_STRING,U_MIN);
  if(designMode == N_OPTIMIZATION) {
    fprintf(out,"%s Design mode: N\n",COMMENT_STRING);
  } else if(designMode == P_OPTIMIZATION) {
    fprintf(out,"%s Design mode: P\n",COMMENT_STRING);
  } else if(designMode == MFE_OPTIMIZATION) {
    fprintf(out,"%s Design mode: MFE\n",COMMENT_STRING);
  } else {
    fprintf(out,"%s Design mode: unsupported\n",COMMENT_STRING);
  }
  fclose(out);
  
}

void outputConstantParameters(char *filename) {

	FILE *paramOut = fopen(filename,"w");

	#ifdef DEBUG
		fprintf(paramOut, "*Running in DEBUG mode*\n");

	#endif
	
	if (loadInit) {
		fprintf(paramOut, "*LOADING INITIAL SEED*\n");	
	}
	else {
		fprintf(paramOut, "INIT_MODE\t%d\n", initMode);	
	}
	
	fprintf(paramOut, "NA_TYPE\t\t%d\n", DNARNACOUNT);
	fprintf(paramOut, "DANGLES\t\t%d\n", DANGLETYPE);
	
	fprintf(paramOut,"M_UNFAVORABLE\t%f\n", M_UNFAVORABLE);
	fprintf(paramOut,"M_LEAFOPT\t%d\n", mLeafopt);
	
	
	if (designMode == N_OPTIMIZATION) {
		fprintf(paramOut, "*N OPTIMIZATION*\n");
		fprintf(paramOut,"N_RATIO\t\t%Lf\n", nRatio);
	}
	else if (designMode == P_OPTIMIZATION) {
		fprintf(paramOut, "*P OPTIMIZATION*\n");
	}
	else if (designMode == MFE_OPTIMIZATION) {
		fprintf(paramOut, "*MFE OPTIMIZATION*\n");

		fprintf(paramOut,"M_EVAL\t\t%d\n", M_EVAL);
		fprintf(paramOut,"B_ACCEPT\t%f\n", B_ACCEPT);
		fprintf(paramOut,"MFE_EPS\t%f\n", MFE_EPS);
	}
	
	if (!bypassHierarchy) {
		fprintf(paramOut, "*Using Hierarchical Design*\n");
		fprintf(paramOut,"H_MIN\t\t%d\n", hMin);
		fprintf(paramOut,"U_MIN\t\t%d\n", U_MIN);
		fprintf(paramOut,"M_REOPT_DEFAULT\t%d\n", M_REOPT_DEFAULT);
	}
	else {
		fprintf(paramOut, "*Using Single-scale Design*\n");
	}
	
	if (!bypassGuidance) {
		fprintf(paramOut, "*Using Guidance*\n");

	}
	else {
		fprintf(paramOut, "*Bypassing Guidance*\n");
	}
  
  if (quickFlag) {
   fprintf(paramOut, "*Running in Qiick Mode*\n"); 
  }

	fclose(paramOut);
}


void initEngine(char *psFile) {
	loadNucNumbers();
	if (psFile[0]) {
		loadBadStrings(psFile);
	}	
	
	if (designMode == P_OPTIMIZATION && !bypassGuidance) {
		bypassGuidance = 1;
		printf("Bypassing guidance for P optimization\n");
	}
	
	if (designMode == P_OPTIMIZATION && !bypassHierarchy) {
		bypassHierarchy = 1;
		printf("Bypassing decomposition for P optimization\n");
	}
  
}


void freeEngine(void) {
	freeLoadedStrings();	
}


void resetOptVal(secStruc *myStruc) {
	myStruc->checked = 0;
	if (designMode == N_OPTIMIZATION || designMode == MFE_OPTIMIZATION) {
		myStruc->optVal = (DBL_TYPE)myStruc->length;
		if (designMode == N_OPTIMIZATION) {
			myStruc->nVal = (DBL_TYPE)myStruc->length;
			myStruc->nDifference = myStruc->length;
			myStruc->bestNDifference = myStruc->length;


		}
		else if (designMode == MFE_OPTIMIZATION) {
			myStruc->mfeDifference = myStruc->length;
			myStruc->bestMFEDifference = myStruc->length;

		}
	}
	else {
		#ifdef DEBUG
			assert(designMode == P_OPTIMIZATION);
		#endif
		myStruc->optVal = 0.0;
	}
}


void resetNodeVals(secStruc *myStruc, DBL_TYPE optVal) {
	myStruc->optVal = optVal;
	if (designMode == N_OPTIMIZATION) {
		myStruc->bestOptVal = myStruc->optVal;
	}
	else if (designMode == MFE_OPTIMIZATION) {
		myStruc->bestMFEDifference = (int)myStruc->optVal;
		myStruc->bestOptVal = myStruc->optVal;
	}
	else {
		#ifdef DEBUG
			assert(designMode == P_OPTIMIZATION);
		#endif
		myStruc->bestPVal = myStruc->optVal;
		myStruc->bestOptVal = myStruc->optVal;
	}
	
	#ifdef DEBUG
		assertInSync(myStruc);
	#endif
}

#ifdef LEAFCORRECTION
void overwriteChildrenLeafArrays(secStruc *myStruc) {
  if(hasChildren(myStruc)) {
    #ifdef DEBUG
    printf("overwriting children\n");
    #endif
    memcpy(myStruc->leftChild->leafArray, myStruc->leafArray, MAX_NODES * sizeof(DBL_TYPE));
    memcpy(myStruc->rightChild->leafArray, myStruc->leafArray, MAX_NODES * sizeof(DBL_TYPE));
    memcpy(myStruc->leftChild->leafNative, myStruc->leafNative, MAX_NODES * sizeof(DBL_TYPE));
    memcpy(myStruc->rightChild->leafNative, myStruc->leafNative, MAX_NODES * sizeof(DBL_TYPE));
    //memcpy(myStruc->leftChild->leafFull, myStruc->leafFull, MAX_NODES * sizeof(DBL_TYPE));
    //memcpy(myStruc->rightChild->leafFull, myStruc->leafFull, MAX_NODES * sizeof(DBL_TYPE));
    memcpy(myStruc->leftChild->leafNativeCheck, myStruc->leafNativeCheck, MAX_NODES * sizeof(int));
    memcpy(myStruc->rightChild->leafNativeCheck, myStruc->leafNativeCheck, MAX_NODES * sizeof(int));
    memcpy(myStruc->leftChild->leafCheck, myStruc->leafCheck, MAX_NODES * sizeof(int));
    memcpy(myStruc->rightChild->leafCheck, myStruc->leafCheck, MAX_NODES * sizeof(int));
    myStruc->leftChild->modified = 0;
    myStruc->rightChild->modified = 0;
  }
}

void copyLeafArrayFromChild(secStruc *parent, secStruc *child) {
  //assert child modified?
  int i = 0;
  for (i = 0; i < MAX_NODES; i++) {
    if (child->leafCheck[i] != -1) {
      parent->leafCheck[i] = child->leafCheck[i];
      parent->leafArray[i] = child->leafArray[i];
    }
    if (child->leafNativeCheck[i] != -1) {
      parent->leafNative[i] = child->leafNative[i];
      //parent->leafFull[i] = child->leafFull[i];
      parent->leafNativeCheck[i] = child->leafNativeCheck[i]; 
    }
  }
}

#endif

void setNewBestN(secStruc *myStruc) {
	int matrixLength = (myStruc->length + 1) * (myStruc->length + 1);
	memcpy(myStruc->bestPr, myStruc->myPairPr, matrixLength *sizeof(DBL_TYPE));
	memcpy(myStruc->bestProbs, myStruc->myProbs, myStruc->length *sizeof(DBL_TYPE));
	memcpy(myStruc->bestConflictArray, myStruc->conflictArray, myStruc->length * sizeof(int));
  #ifdef LEAFCORRECTION
  if (hasChildren(myStruc)) {
    if (myStruc->leftChild->modified) {
      // secStruc *leftChild = myStruc->leftChild;
      copyLeafArrayFromChild(myStruc,myStruc->leftChild);

    }
    if (myStruc->rightChild->modified) {

      copyLeafArrayFromChild(myStruc,myStruc->rightChild);
    }
    // secStruc *leftChild = myStruc->leftChild;

    
    int i = 0;
    DBL_TYPE *tempValues = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
    getNValues(myStruc, tempValues);
    DBL_TYPE sum = 0.0;
    int checksum = 0;
    
    for (i = 0; i < myStruc->length; i++) {
      if (myStruc->parentMapping[i] >= 0) {
        sum += tempValues[i];
        checksum += myStruc->seq[i];
      }
      //fullSum += tempValues[i];
    }
    myStruc->leafNative[myStruc->id] = sum;
    myStruc->leafNativeCheck[myStruc->id] = checksum;
    //myStruc->leafFull[myStruc->id] = fullSum;
    assertNativeChecksum(myStruc);
  }
  #endif
}


void setNewBestMFE(secStruc *myStruc) {
	memcpy(myStruc->bestPairs, myStruc->mfePairs, myStruc->length * sizeof(int));
	memcpy(myStruc->bestConflictArray, myStruc->conflictArray, myStruc->length * sizeof(int));
	memcpy(myStruc->bestProbs, myStruc->myProbs, myStruc->length *sizeof(DBL_TYPE));	
	myStruc->bestMFEDifference = myStruc->mfeDifference;
  #ifdef LEAFCORRECTION
  if (hasChildren(myStruc)) {
    if (myStruc->leftChild->modified) {
      // secStruc *leftChild = myStruc->leftChild;
      copyLeafArrayFromChild(myStruc,myStruc->leftChild);

    }
    if (myStruc->rightChild->modified) {

      copyLeafArrayFromChild(myStruc,myStruc->rightChild);
    }
    // secStruc *leftChild = myStruc->leftChild;

    
    int i = 0;
    DBL_TYPE *tempValues = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
    getNValues(myStruc, tempValues);
    DBL_TYPE sum = 0.0;
    int checksum = 0;
    
    for (i = 0; i < myStruc->length; i++) {
      if (myStruc->parentMapping[i] >= 0) {
        sum += tempValues[i];
        checksum += myStruc->seq[i];
      }
      //fullSum += tempValues[i];
    }
    myStruc->leafNative[myStruc->id] = sum;
    myStruc->leafNativeCheck[myStruc->id] = checksum;
    //myStruc->leafFull[myStruc->id] = fullSum;
    assertNativeChecksum(myStruc);
  }
  #endif	
}


void setNewBestP(secStruc *myStruc) {
	myStruc->bestPVal = myStruc->pVal;	
}


#ifdef DEBUG
void assertInSync(secStruc *myStruc) {
	if (mReopt > 0) {
		if (!hasChildren(myStruc)) {
				assert(seqCmp(myStruc->seq, myStruc->bestSeq, myStruc->length) == 0);
				assert(floatsEqual(myStruc->optVal,myStruc->bestOptVal));
		}
		
		if (designMode == N_OPTIMIZATION) {
			
		}
		else if (designMode == MFE_OPTIMIZATION) {
			
		}
		else {
			assert(designMode == P_OPTIMIZATION);
		}
	}
}
#endif

#ifdef LEAFCORRECTION
DBL_TYPE sumLeafIDs(secStruc *myStruc) {
  int i = 0;
  DBL_TYPE sum = 0.0;
  for (i = 0; i < MAX_NODES; i++) {
    sum += myStruc->leafArray[i];
  }
  assertChecksum(myStruc);
  return sum;
}


//DBL_TYPE sumNativeIDs(secStruc *myStruc) {
//  int i = 0;
//  DBL_TYPE sum = 0.0;
//  for (i = 0; i < MAX_NODES; i++) {
//    sum += myStruc->leafNative[i];
//  }
// assertNativeChecksum(myStruc);
//  return sum;
//}


void assertChecksum(secStruc *myStruc) {
  // don't do this if pattern prevention is on
  if (maxBadStrings > 0) {
    return; 
  }
  int i = 0;
  int checksum = 0;
  for (i = 0; i < MAX_NODES; i++) {
    if (myStruc->leafCheck[i] >= 0) {
      checksum += myStruc->leafCheck[i];
    }
  }
  int myCheck = 0;
  if (checksum > 0) {
    for (i = 0; i < myStruc->length; i++) {
     myCheck += myStruc->seq[i] * myStruc->purelyNativeBases[i];
    }

    assert(myCheck == checksum);
  }

}


void assertNativeChecksum(secStruc *myStruc) {

  
  // don't do this if pattern prevention is on
  if (maxBadStrings > 0) {
    return; 
  }
  int i = 0;
  int checksum = 0;
  
  if (hasChildren(myStruc)) {
     checksum = myStruc->leafNativeCheck[myStruc->leftChild->id] + myStruc->leafNativeCheck[myStruc->rightChild->id];
     #ifdef DEBUG
     printf("CHILDCHECK: %d:%d,%d\n",myStruc->length, myStruc->leafNativeCheck[myStruc->leftChild->id], myStruc->leafNativeCheck[myStruc->rightChild->id]);
     #endif     
  }
  else {
     checksum = myStruc->leafNativeCheck[myStruc->id];
  }
  
  int myCheck = 0;
  if (checksum > 0) {
    for (i = 0; i < myStruc->length; i++) {
      if (hasChildren(myStruc)) {
        myCheck += myStruc->seq[i];
      }
      else {
       if (myStruc->parentMapping[i] >= 0) {
        myCheck += myStruc->seq[i];  
       }
      }
      
    }
    #ifdef DEBUG
    printf("Native CHECKSUM TEST: %d:%d\n",myCheck,checksum);
    #endif
    assert(myCheck == checksum);
  }

}



#endif

void extractBestValues(secStruc *myStruc, DBL_TYPE *bestVal, int *bestSeq, DBL_TYPE *bestPr, DBL_TYPE *bestProbs) {
	memcpy(bestVal, &myStruc->bestOptVal, sizeof(DBL_TYPE));
	memcpy(bestSeq, myStruc->bestSeq, sizeof(int)*myStruc->length);
	memcpy(bestPr, myStruc->bestPr, sizeof(DBL_TYPE)*(myStruc->length + 1)*(myStruc->length+1));
	memcpy(bestProbs, myStruc->bestProbs, sizeof(DBL_TYPE)*(myStruc->length));
  
  #ifdef LEAFCORRECTION
  int i;
  DBL_TYPE *values = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
  getNValues(myStruc,values);
  DBL_TYPE sum = 0.0;
  int checksum = 0;
  int checksumNative = 0;
  DBL_TYPE nativeSum = 0.0;
  for (i = 0; i < myStruc->length; i++) {
    sum += values[i] * myStruc->purelyNativeBases[i];
    checksum += myStruc->seq[i] * myStruc->purelyNativeBases[i];
    if (myStruc->parentMapping[i] >= 0) {
      nativeSum += values[i];
      checksumNative += myStruc->bestSeq[i];
    }
  }
  myStruc->leafArray[myStruc->id] = sum;
  myStruc->leafNative[myStruc->id] = nativeSum;
  myStruc->leafCheck[myStruc->id] = checksum;
  myStruc->leafNativeCheck[myStruc->id] = checksumNative;
  assertNativeChecksum(myStruc);
  free(values);
  #endif
	
}


void setAndRevertToOldValue(secStruc *myStruc, DBL_TYPE bestVal, int *bestSeq, DBL_TYPE *bestPr, DBL_TYPE *bestProbs) {
	int matrixLength = (myStruc->length + 1) * (myStruc->length + 1);
	
	memcpy(myStruc->bestPr, bestPr, matrixLength *sizeof(DBL_TYPE));
	memcpy(myStruc->bestProbs, bestProbs, myStruc->length *sizeof(DBL_TYPE));

	memcpy(myStruc->bestSeq, bestSeq,sizeof(int)*myStruc->length);
 
	myStruc->bestOptVal = bestVal;
	
	revertToOldValue(myStruc);

	#ifdef DEBUG
		assertInSync(myStruc);
	#endif
	
	
}

void revertToOldValue(secStruc *myStruc) {
	int matrixLength = (myStruc->length + 1) * (myStruc->length + 1);
	memcpy(myStruc->myPairPr, myStruc->bestPr, matrixLength *sizeof(DBL_TYPE));
	memcpy(myStruc->myProbs, myStruc->bestProbs, myStruc->length *sizeof(DBL_TYPE));
  memcpy(myStruc->mfePairs, myStruc->bestPairs, myStruc->length *sizeof(int));
	memcpy(myStruc->seq, myStruc->bestSeq,sizeof(int)*myStruc->length);
 
	myStruc->optVal = myStruc->bestOptVal;
	
	#ifdef DEBUG
		assert(myStruc->optVal < (DBL_TYPE)myStruc->length);
	#endif
	memcpy(myStruc->conflictArray, myStruc->bestConflictArray, sizeof(int)*myStruc->length);
	
	if (designMode == N_OPTIMIZATION) {
		myStruc->nVal = myStruc->optVal;
	}
	else if (designMode == MFE_OPTIMIZATION) {
		myStruc->mfeDifference = myStruc->bestMFEDifference;
	}
	else {
		#ifdef DEBUG
			assert(designMode == P_OPTIMIZATION);
		#endif
		myStruc->pVal = myStruc->bestPVal;
	}
	myStruc->checked = 1;
	
  #ifdef LEAFCORRECTION
  overwriteChildrenLeafArrays(myStruc);
  #endif
  
	#ifdef DEBUG
		assertInSync(myStruc);
	#endif
	
	
}


void setNewBestValue(secStruc *myStruc) {
	if (designMode == N_OPTIMIZATION) {
		setNewBestN(myStruc);
	}
	else if (designMode == MFE_OPTIMIZATION) {
		setNewBestMFE(myStruc);
		myStruc->bestOptVal = myStruc->bestMFEDifference;

	}
	else {
		#ifdef DEBUG
			assert(designMode == P_OPTIMIZATION);
		#endif
		setNewBestP(myStruc);
		myStruc->bestOptVal = myStruc->bestPVal;
	}
	myStruc->bestOptVal = myStruc->optVal;

	memcpy(myStruc->bestSeq, myStruc->seq, sizeof(int)*(myStruc->length + 1));
	#ifdef DEBUG
		assertInSync(myStruc);
	#endif
	
}


/* 	this function ensures that all ways of evaluating a sequence are synchronized
the bestSeq and seq should already be in sync
*/
void setAllBestValues(secStruc *myStruc) {

	setNewBestN(myStruc);
	setNewBestP(myStruc);
	setNewBestMFE(myStruc);
	
	revertToOldValue(myStruc);
	
	#ifdef DEBUG
		assertInSync(myStruc);	
	#endif
}


void freeStruc(secStruc *myStruc) {
	if (myStruc->leftChild) {
		freeStruc(myStruc->leftChild);
		free(myStruc->leftChild);
		free(myStruc->leftChildMapping);	
	}
	if (myStruc->rightChild) {
		freeStruc(myStruc->rightChild);
		free(myStruc->rightChild);
		free(myStruc->rightChildMapping);
	}
	free(myStruc->tempPairing);
	free(myStruc->tempFolding);	
	free(myStruc->expFolding);

	free(myStruc->seq);
	free(myStruc->skel);
	free(myStruc->bestSeq);
	free(myStruc->fakeBases);
  free(myStruc->purelyNativeBases);


	free(myStruc->parentMapping); 
	free(myStruc->bestPr);
	free(myStruc->bestProbs);
	free(myStruc->myPairPr);
	free(myStruc->myProbs);
	free(myStruc->bestPairs);
	free(myStruc->mfePairs);
	free(myStruc->conflictArray);
	free(myStruc->bestConflictArray);
	free(myStruc->strandBreaks);
	
	
	if (!myStruc->parent) {
			free(myStruc->initNucs);	
			free(myStruc->initComps);
			free(myStruc->initMap);	
	}
}


void fixPatternViolations(secStruc *myStruc) {
	if (maxBadStrings < 1) {
		return;
	}
	
	int currentViolations = violatesPattern(myStruc->seq, myStruc->skel, myStruc->length * maxBadStringLength);

	while (currentViolations > 0) {

		mutation candidateMutation;
		mutateRandomBaseOrPair(myStruc, &candidateMutation);
		if (!mutationProhibited(myStruc, &candidateMutation, &currentViolations, currentViolations)) {
			applyMutation(myStruc, &candidateMutation);
		}
	}
	
	#ifdef DEBUG
		assert(!violatesPattern(myStruc->seq, myStruc->skel,0));
	#endif
}


int nearStrandBreak(int pos, int length, int *pairing, int numStrands, int *breaks, int hval) {
	if (numStrands <= 1) {
		#ifdef DEBUG
			assert(numStrands == 1);
		#endif
		return 0;
	}
	#ifdef DEBUG
		assert(pairing[pos] > pos);
		assert(numStrands > 0);
	#endif
	
	int i;
	for (i = 0; i < numStrands - 1; i++) {
		int startPoint = pos - hval + 1;
		int endPoint = pos + hval;
		#ifdef DEBUG
			assert(startPoint >= 0);
			assert(endPoint < length);
			assert((breaks[i] - i) < length && (breaks[i] - i) >= hval);
		#endif
		
		int breakPoint = breaks[i] - i;
		if (breakPoint >= startPoint && breakPoint <= endPoint) {
			return 1;	
		}
		
		startPoint = pairing[pos] - hval;
		endPoint = pairing[pos] + hval - 1;
		#ifdef DEBUG
			assert(startPoint >= 0);
			assert(endPoint < length);
		#endif
		if (breakPoint >= startPoint && breakPoint <= endPoint) {
			return 1;	
		}
	}
	return 0;
}


int findSplitPoint(int *pairing, int length, int numStrands, int *breaks) {

	
	if (bypassHierarchy) {
		return -1;
	}

	int smallestPoint = -9999;
	#ifdef DEBUG
		int sizeSatisfied = 0;
	#endif
	
	int i;
	int smallestDistance = length * length;

	for (i = hMin; i < length - hMin; i++) {
		if (pairing[i] > i && pairing[i] < (length - hMin) &&
			(pairing[i] - i - 1 > U_MIN) && (length - pairing[i] + i + 1 > U_MIN) && 
			!nearStrandBreak(i, length, pairing, numStrands, breaks, hMin)) {
			#ifdef DEBUG
				assert(pairing[i] >= 0);
			#endif
			int j;
			int found = 1;
			for (j = 1; j < hMin; j++) {
				#ifdef DEBUG
					assert(i - j >= 0);
					assert(i + j < length);
				#endif
				if (pairing[i - j] < 0 || pairing[i + j] < 0 || pairing[i] + j < 0 || pairing[i] - j < 0) {
					found = 0;
				}
				if (pairing[i - j] != (pairing[i] + j) || pairing[i + j] != (pairing[i] - j)) {
					found = 0;
				}
			}
			#ifdef DEBUG
				assert(i + hMin < length);
			#endif
			if (pairing[i + hMin] != (pairing[i] - hMin)) {
				found = 0;
			}
			#ifdef DEBUG
				if (found) {
					assert(pairing[i + hMin] >= 0);
					assert(pairing[i] - hMin >= 0);
				}
					
			#endif
			if (found) {
				int possSplitPointA = i;
				int possSplitPointB = pairing[i];
				
				int rightSplitLength = possSplitPointB - possSplitPointA - 1;
				int leftSplitLength = length - rightSplitLength;

				#ifdef DEBUG
					int z;
					int checkLeftLength = 0;
					int checkRightLength = 0;
					for (z = 0; z < possSplitPointA + 1; z++) {
						checkLeftLength++;	
					}
					for (z = possSplitPointA + 1; z < possSplitPointB; z++) {
						checkRightLength++;	
					}
					for (z = possSplitPointB; z < length; z++) {
						checkLeftLength++;	
					}

					assert(checkLeftLength == leftSplitLength);
					assert(checkRightLength == rightSplitLength);
				
				#endif
				
				int diff = (rightSplitLength - leftSplitLength);
				int diffSq = diff * diff;
				
				
				#ifdef DEBUG
					assert(rightSplitLength > U_MIN && leftSplitLength > U_MIN);
				#endif

				if (diffSq < smallestDistance) {
					smallestDistance = diffSq;
					smallestPoint = i;
				}
			}
		}
	}
	

	if (smallestPoint >= 0) {
		#ifdef DEBUG
			if (!sizeSatisfied) {
				assert(smallestPoint > 0 && smallestPoint < (length * length));
			}
			sizeSatisfied = 1;
		#endif

	}
	return smallestPoint + 1;
}


void initializeMapping(int *mapping, int length) {
	int i;
	for (i = 0; i < length; i++) {
		mapping[i] = -99;
	}
}



int mapFakeBases(secStruc *parent, secStruc *child) {
	int realBaseCount = child->length;
	int childLength = child->length;
	
	int *childToParentMapping = child->parentMapping;
	int *childFakeBases = child->fakeBases;
	int *parentFakeBases = parent->fakeBases;
	
	int i;
	for (i = 0; i < childLength; i++) {
		
		if (childToParentMapping[i] < 0) {
			childFakeBases[i] = 1;
			
		}
		else if (!child->amILeft && (i < hMin || i > childLength - hMin - 1)){
			childFakeBases[i] = 1;
			
		}
		else if (child->amILeft && (i >= child->loopStart - hMin && i < child->loopStart + hMin)) {
			childFakeBases[i] = 1;
			
		}
		else {
			childFakeBases[i] = parentFakeBases[childToParentMapping[i]];
			
		}
		realBaseCount -= childFakeBases[i];
	}

  
  #ifdef DEBUG
  int pncounter = 0;
  #endif
  for (i = 0; i < childLength; i++) {
    child->purelyNativeBases[i] = 1;
    if (childFakeBases[i] == 1) {
      child->purelyNativeBases[i] = 0;
    }
    if (child->parentMapping[i] >=0 && parent->purelyNativeBases[child->parentMapping[i]] == 0) {
      child->purelyNativeBases[i] = 0; 
    }
    #ifdef DEBUG
    if (child->purelyNativeBases[i] == 1) pncounter++;
    #endif

  }
  
  #ifdef DEBUG
  assert(pncounter > 0);
  #endif
  
	return realBaseCount;
}

#ifdef DEBUG
void assertValidStrandBreaks(secStruc *myStruc) {
		int i;
		for (i = 0; i < myStruc->numStrands - 1; i++) {
			assert(myStruc->strandBreaks[i] > hMin);
			assert(myStruc->strandBreaks[i] < myStruc->length - hMin);
		}
}
#endif


void initStruc(secStruc *myStruc, int length, int numStrands) {
	myStruc->parent = 0;
	myStruc->checked = 0;
	myStruc->timesOptimized = 0;
	myStruc->totalTimeSpent = 0;
	myStruc->pfTimeSpent = 0.0;
	myStruc->mfeTimeSpent = 0.0;
	myStruc->pTimeSpent = 0.0;

	myStruc->nsCalls = 0;
	myStruc->failures = 0;
	myStruc->length = length;
	if (numStrands > 0) {
		myStruc->expandedLength = length + numStrands - 1;
	}
	else {
		myStruc->expandedLength = length;	
	}
	myStruc->numStrands = numStrands;
	
	
	myStruc->bestSeq = (int *)calloc(myStruc->length + 1, sizeof(int));
	myStruc->bestSeq[myStruc->length] = -1;
	
	myStruc->seq = (int *)calloc(myStruc->length + 1, sizeof(int));
	myStruc->seq[myStruc->length] = -1;
	
	
	if (!myStruc->skel) {
		myStruc->skel = (int *)calloc(myStruc->expandedLength + 1, sizeof(int));
		myStruc->skel[myStruc->expandedLength] = -1;
		initSeqToN(myStruc->skel, myStruc->expandedLength);
	}
	
	if (!myStruc->tempFolding) {
		myStruc->tempFolding = (char *)calloc(myStruc->expandedLength + 1, sizeof(char));
		#ifdef DEBUG
			assert(!myStruc->tempPairing);
		#endif
		
		myStruc->tempPairing = (int *)calloc(myStruc->expandedLength + 1, sizeof(int));
	}
	
	initSeqToN(myStruc->seq, myStruc->length);
	
	myStruc->leftChild = 0;
	myStruc->rightChild = 0;
	
	myStruc->parentMapping = (int *)calloc(myStruc->length, sizeof(int));
	myStruc->bestPr = (DBL_TYPE *) calloc( (myStruc->length + 1) * (myStruc->length + 1),sizeof( DBL_TYPE) );
	myStruc->bestProbs = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
	myStruc->myPairPr = (DBL_TYPE *)calloc((myStruc->length + 1) * (myStruc->length + 1), sizeof(DBL_TYPE));
	myStruc->myProbs = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
	myStruc->bestPairs = (int *)calloc(myStruc->length, sizeof(int));
	
	myStruc->mfePairs = (int *)calloc(myStruc->length, sizeof(int));
	myStruc->conflictArray = (int *)calloc(myStruc->length, sizeof(int));
	myStruc->bestConflictArray = (int *)calloc(myStruc->length, sizeof(int));
	resetOptVal(myStruc);
	myStruc->mfeDifference = myStruc->length;
	myStruc->bestOptVal = (DBL_TYPE)myStruc->length;
	
	myStruc->nVal = (DBL_TYPE)myStruc->length;
	myStruc->bestNVal = (DBL_TYPE)myStruc->length;

	myStruc->bestMFEDifference = myStruc->length;
	myStruc->bestPVal = 0.0;

	
	if (designMode == N_OPTIMIZATION || designMode == MFE_OPTIMIZATION) {
		myStruc->targetObjective = myStruc->length * nRatio;
	}
	else {
		#ifdef DEBUG
			assert(designMode == P_OPTIMIZATION);
		#endif
		myStruc->targetObjective = 1.0 - nRatio;
	}
  
  #ifdef LEAFCORRECTION
  myStruc->leafArray = (DBL_TYPE *)calloc(MAX_NODES,sizeof(DBL_TYPE)); 
  myStruc->leafCheck = (int *)calloc(MAX_NODES,sizeof(int)); 
  myStruc->leafNative = (DBL_TYPE *)calloc(MAX_NODES,sizeof(DBL_TYPE)); 
  myStruc->leafNativeCheck = (int *)calloc(MAX_NODES,sizeof(int)); 
  int k;
  for (k = 0; k < MAX_NODES; k++) {
    myStruc->leafCheck[k] = -1;
    myStruc->leafNativeCheck[k] = -1;
  }
  #endif
}


void createTree(secStruc *parentStruc) {
  
	int *pairing = parentStruc->tempPairing;
	int length;
	
	length = parentStruc->length;
	
	int bestSplitPointA = 0;
	int bestSplitPointB = length - 1;
	
	
	
	bestSplitPointA = findSplitPoint(pairing, length, parentStruc->numStrands, parentStruc->strandBreaks);

	if (bestSplitPointA >= 0) {
		bestSplitPointB = pairing[bestSplitPointA];
		
		#ifdef DEBUG
			assert(bestSplitPointB > bestSplitPointA);

		#endif
		
		secStruc *rightStruc = (secStruc *)calloc(1,sizeof(secStruc));
		
		int fakeBaseCount = 2 * hMin; 

		int rightLength = bestSplitPointB - bestSplitPointA + fakeBaseCount + 1;
		
		#ifdef DEBUG
			int c;
			int checkRight = hMin * 2;
			for (c = bestSplitPointA; c < bestSplitPointB + 1; c++) {
				checkRight++;
			}
			assert(checkRight == rightLength);
		
		#endif
		
		initStruc(rightStruc, rightLength, -1);
		rightStruc->parent = parentStruc;
		rightStruc->initMap = parentStruc->initMap;
		rightStruc->initComps = parentStruc->initComps;
		rightStruc->ssmWordLength = parentStruc->ssmWordLength;
		
		
		rightStruc->amILeft=0;
		rightStruc->level = parentStruc->level + 1;
		
		rightStruc->k = parentStruc->k + 1;
		rightStruc->l = parentStruc->l * 2 + 1;
    #ifdef LEAFCORRECTION
    rightStruc->id = ++nodeCounter;
    #else
    rightStruc->id = 2 * parentStruc->id + 2;
    #endif
		#ifdef DEBUG
    printf("NODE: %d\n",rightStruc->id);
    assert(rightStruc->id < MAX_NODES);
    #endif
    
    
		rightStruc->initNucs = parentStruc->initNucs;
		rightStruc->numInitBases = parentStruc->numInitBases;
		
		parentStruc->rightChildMapping = (int *)calloc(length, sizeof(int));
		
		initializeMapping(rightStruc->parentMapping, rightStruc->length);
		initializeMapping(parentStruc->rightChildMapping, length);
		
		parentStruc->rightChildMapping[0] = -1000;
		
		rightStruc->fakeBases = (int *)calloc(rightStruc->length, sizeof(int));
    rightStruc->purelyNativeBases = (int *)calloc(rightStruc->length, sizeof(int));

		
		int l;
		
		int startPoint = bestSplitPointA - (fakeBaseCount / 2);
		
		for (l = 0; l < rightStruc->length; l++) {
			
			parentStruc->rightChildMapping[startPoint + l] = l;
			
							
			if (l > hMin - 1 && l < rightStruc->length - hMin) {
				rightStruc->parentMapping[l] = startPoint + l;
			}
			else {
				
				rightStruc->parentMapping[l] = -1;		
			}
			
			if (pairing[startPoint + l] >= 0) {
				rightStruc->tempPairing[l] = pairing[startPoint + l] - startPoint;
			}
			else {
				rightStruc->tempPairing[l] = -1;
			}
			
		}
		
		mapSequenceDown(parentStruc, rightStruc, 1);
		
		int z;
		for (z = 0; z < rightStruc->length; z++) {
			if (rightStruc->tempPairing[z] >= 0) {
				if (rightStruc->tempPairing[z] > z) {
					rightStruc->tempFolding[z] = '(';
				}
				else {
					rightStruc->tempFolding[z] = ')';
				}
			}
			else {
				rightStruc->tempFolding[z] = '.';
			}
		}
		
		rightStruc->indyLength = 0;
		for (z = 0; z < rightStruc->length; z++) { 
			if (rightStruc->tempPairing[z] >= z || rightStruc->tempPairing[z] < 0) { 
				rightStruc->indyLength++;
			}
		}
		
		rightStruc->relevantLength = mapFakeBases(parentStruc, rightStruc);

		secStruc *leftStruc = (secStruc *)calloc(1,sizeof(secStruc));
		
		
		int leftLength = length - bestSplitPointB + bestSplitPointA - 1 + fakeBaseCount;
		
	
		#ifdef DEBUG
			int checkLeft = hMin * 2;
			
			for (c = 0; c < bestSplitPointA; c++) {
				checkLeft++;
			}
			for (c = bestSplitPointB + 1; c < length; c++) {
				checkLeft++;
			}
			assert(checkLeft == leftLength);
		
		#endif

		initStruc(leftStruc, leftLength, -1);
		leftStruc->loopStart = bestSplitPointA + hMin;
		leftStruc->amILeft = 1;
		leftStruc->parent = parentStruc;
		leftStruc->initMap = parentStruc->initMap;
		leftStruc->initComps = parentStruc->initComps;
		leftStruc->ssmWordLength = parentStruc->ssmWordLength;

		leftStruc->level = parentStruc->level + 1;
		parentStruc->leftChildMapping = (int *)calloc(length, sizeof(int));
		
		leftStruc->k = parentStruc->k + 1;
		leftStruc->l = parentStruc->l * 2;
    #ifdef LEAFCORRECTION
    leftStruc->id = ++nodeCounter;
    #else
    leftStruc->id = 2 * parentStruc->id + 1;
    #endif
    
    #ifdef DEBUG
    printf("NODE: %d\n",leftStruc->id);

    assert(leftStruc->id < MAX_NODES);
    #endif
    
		
		initializeMapping(leftStruc->parentMapping, leftStruc->length);
		initializeMapping(parentStruc->leftChildMapping, length);
		
		leftStruc->fakeBases = (int *)calloc(leftStruc->length, sizeof(int));
 		leftStruc->purelyNativeBases = (int *)calloc(leftStruc->length, sizeof(int));

		
		leftStruc->initNucs = parentStruc->initNucs;
		leftStruc->numInitBases = parentStruc->numInitBases;
		
		memcpy(leftStruc->tempPairing, parentStruc->tempPairing, (bestSplitPointA + (fakeBaseCount / 2)) * sizeof(int));
		memcpy(&leftStruc->tempPairing[bestSplitPointA + (fakeBaseCount / 2)], &parentStruc->tempPairing[bestSplitPointB + 1 - (fakeBaseCount / 2)], (length - bestSplitPointB + (fakeBaseCount / 2) - 1) * sizeof(int));
		
		for (l = 0; l < leftStruc->length; l++) {
			if (l < leftStruc->loopStart) {

				if (l < leftStruc->loopStart - hMin) {
					leftStruc->parentMapping[l] = l;
				}
				else {
					leftStruc->parentMapping[l] = -1;	
				}
				parentStruc->leftChildMapping[l] = l;
			}
			else {
				if (l >= leftStruc->loopStart + hMin) {
					leftStruc->parentMapping[l] = (bestSplitPointB - hMin)  + (l - (bestSplitPointA + hMin - 1));

				}
				else {
					leftStruc->parentMapping[l] = -1;
				}
				parentStruc->leftChildMapping[(bestSplitPointB - hMin)  + (l - (bestSplitPointA + hMin - 1))] = l; 
			}
			
			if (leftStruc->tempPairing[l] > bestSplitPointA + (fakeBaseCount / 2)) {
				leftStruc->tempPairing[l] = leftStruc->tempPairing[l] - (bestSplitPointB - bestSplitPointA) + fakeBaseCount - 1;
			}
			
		}
		
		for (z = 0; z < leftStruc->length; z++) {
			if (leftStruc->tempPairing[z] >= 0) {
				if (leftStruc->tempPairing[z] > z) {
					leftStruc->tempFolding[z] = '(';
				}
				else {
					leftStruc->tempFolding[z] = ')';
				}
			}
			else {
				leftStruc->tempFolding[z] = '.';
			}
		}
		
		
		#ifdef DEBUG
		printf("%s\n",leftStruc->tempFolding);
		for (l = 0; l < leftStruc->length; l++) {
			printf("%d ",leftStruc->parentMapping[l]);			
		}
		printf("\n");
		#endif
		
		leftStruc->indyLength = 0;
		for (z = 0; z < leftStruc->length; z++) { 
			if (leftStruc->tempPairing[z] >= z || leftStruc->tempPairing[z] < 0) { 
				leftStruc->indyLength++;
			}
		}
		
		leftStruc->relevantLength = mapFakeBases(parentStruc, leftStruc);
	
		mapSequenceDown(parentStruc, leftStruc, 1);
		
		parentStruc->leftChild = leftStruc;
		parentStruc->rightChild = rightStruc;
		
		mapStrandBreaksDown(parentStruc, leftStruc);
		mapStrandBreaksDown(parentStruc, rightStruc);

		leftStruc->expandedSeq = (int *)calloc(leftStruc->length + leftStruc->numStrands, sizeof(int));
		rightStruc->expandedSeq = (int *)calloc(rightStruc->length + rightStruc->numStrands, sizeof(int));

		
		#ifdef DEBUG
			assert(leftStruc->numStrands > 0);
			assert(rightStruc->numStrands > 0);
			assertValidStrandBreaks(leftStruc);
			assertValidStrandBreaks(rightStruc);
      int b;
      int pcount = 0;
      int lcount = 0;
      int rcount = 0;
      for (b = 0; b < parentStruc->length; b++) {
        pcount += parentStruc->purelyNativeBases[b];
      }
      for (b = 0; b < leftStruc->length; b++) {
        lcount += leftStruc->purelyNativeBases[b];
      }
      for (b = 0; b < rightStruc->length; b++) {
        rcount += rightStruc->purelyNativeBases[b];
      }
      assert(pcount == lcount + rcount);

		#endif
		
		if (leftStruc->length >= 2 * U_MIN) {
			createTree(leftStruc);
		}
		
		if (rightStruc->length >= 2 * U_MIN) {
			createTree(rightStruc);
		}
	}
}


void parseStructure(secStruc *myStruc, char *foldFile) {

	char *filename = (char *)calloc(1000,sizeof(char));
	sprintf(filename,"%s.fold",foldFile);
	
	FILE *in = fopen(filename, "r");
	if(in==NULL) {
		fprintf(stderr, "Error: can't open file: %s.\n", filename);
		exit(1);
	}
	
	free(filename);
	
	char *buffer = (char *)calloc(MAX_SEQ_LENGTH,sizeof(char));
	
	fscanf(in, "%s\n", buffer);
	int buffLength = strlen(buffer);
	
	
	int i;
	for (i = 0; i < buffLength; i++) {
		if (buffer[i] != '(' &&
			buffer[i] != ')' &&
			buffer[i] != '.' &&
			buffer[i] != '+') {
		fprintf(stderr, "INVALID STRUCTURE SPECIFIED IN FILE %s:\n ", foldFile);
		exit(1);
			}
	}
	
	myStruc->initNucs = (int*)calloc(NUM_BASES,sizeof(int));
	myStruc->initComps = (int*)calloc(NUM_BASES,sizeof(int));
	myStruc->initMap = (int *)calloc(NUM_BASES, sizeof(int));
	myStruc->strandBreaks = (int *)calloc(buffLength, sizeof(int));
	
	myStruc->tempFolding = (char *)calloc(buffLength + 1, sizeof(char));
	sprintf(myStruc->tempFolding,"%s", buffer);

	myStruc->expFolding = (char *)calloc(buffLength + 1, sizeof(char));
	sprintf(myStruc->expFolding,"%s", buffer);

	
	myStruc->tempPairing = (int *)calloc(buffLength, sizeof(int));
	int numStrands = convertParensToPairs(myStruc->tempFolding, myStruc->tempPairing, buffLength, myStruc->strandBreaks);
	
	#ifdef DEBUG
		assert(numStrands > 0);
	#endif
	
	int length = strlen(myStruc->tempFolding);

	initStruc(myStruc, length, numStrands);
	myStruc->relevantLength = length;
	myStruc->expandedSeq = (int *)calloc(length + numStrands, sizeof(int));

	
	#ifdef DEBUG
		assert(buffLength == myStruc->length + myStruc->numStrands - 1);
	#endif

	fscanf(in, "%s\n", buffer);
	fclose(in);
	
	buffLength = strlen(buffer);
	
	#ifdef DEBUG
		assert(buffLength == myStruc->length + myStruc->numStrands - 1); 
	#endif
	convertBasesToInts(buffer, myStruc->skel);
	

	if (!fixSkeleton(myStruc->skel, myStruc->tempPairing, buffLength)){
		fprintf(stderr, "ERROR IN SEQUENCE CONSTRAINTS\n");
		exit(1);
	}

	myStruc->level = 1;	
	myStruc->fakeBases = (int *)calloc(myStruc->length, sizeof(int));
	myStruc->purelyNativeBases = (int *)calloc(myStruc->length, sizeof(int));

  for (i = 0; i < myStruc->length; i++) {
    myStruc->purelyNativeBases[i] = 1;
  }
  
	
	myStruc->indyLength = getIndependentLength(myStruc->tempPairing, myStruc->length);
	
	setupInitNucs(myStruc);


	myStruc->k = 0;
	myStruc->l = 0;
	myStruc->id = 0;
  
	if (initMode == SSM_INIT) {
		myStruc->ssmWordLength = ssmWordLength(myStruc->length);
	}
	
  
  
	createTree(myStruc);
	
	free(buffer);
}  


long int q2d(int *entry, int length){
	long int decVal = 0;
	int i;
	for (i = 0; i < length; i++) {
		decVal += (long int)pow(4, i) * (long int)(entry[i] - 1);
	}
	return decVal;
}


void clearHash(long int *hash, int wordLength) {
	int i;
	for (i = 0; i < (long int)pow(4,wordLength); i++) {
		hash[i] = -1;
	}
}


long int inHash(long int *hash,int *entry, int wordLength) {
	return hash[q2d(entry, wordLength)];
}


void removeFromHash(long int *hash, int *entry, int wordLength) {
	hash[q2d(entry, wordLength)] = -1;;
}


int invalidComplementInHash(long int *hash, int *word, int pos, int *pairing, int wordLength) {
	int retVal = 0;
	
	int *reverse = reverseStr(word, wordLength);
	
	int *comp = (int *)calloc(wordLength + 1, sizeof(int));
	getComplement(reverse, comp, wordLength);

	int cpos = inHash(hash, comp, wordLength);

	if (cpos >= 0) {
		int i;

		for (i = 0; i < wordLength; i++) {
			if (pairing[pos + i] < 0) {
				retVal = 1;
			}
			else if (pairing[pos+i] != (cpos + wordLength - 1 - i)) {
				retVal = 1;
			}
		}
	}

	free(comp);
	free(reverse);
	
	return retVal;	
}


void addToHash(long int *hash, int *entry, int pos, int wordLength) {
	#ifdef DEBUG
		int i;
		for (i = 0; i < wordLength; i++) {
			assert(entry[i] != BASE_N);	
		}
	#endif
	
	hash[q2d(entry, wordLength)] = pos;
}


void resetTriedOptions(int *tried, int numberOfOptions, int *seq, int *skel, int pos) {
	int j;
	
	for (j = 0; j < numberOfOptions; j++) {
		tried[j] = 0;	
	}
	if (!isAFamily(skel[pos])) {
		seq[pos] = skel[pos];
	}
	else {
		seq[pos] = BASE_N;	
	}
}


int optionsAtThisPosition(int *tried, int numberOfPossibleOptions) {

	int i;
	for (i = 0; i < numberOfPossibleOptions; i++) {
		if (!tried[i]) {
			return 1;	
		}
	}
	return 0;
}


int violatesSSM(secStruc *myStruc, int wordLength) {
	int *pairing = myStruc->tempPairing;
	
	long int *hash = (long int *)malloc((long int)pow(4, wordLength) * sizeof(long int));
	clearHash(hash, wordLength);

	int i;
	
	int *word = (int *)calloc(wordLength + 1, sizeof(int));

		
	#ifdef DEBUG
	if (wordLength == 6) {
		
		word[0] = BASE_A;
		word[1] = BASE_C;
		word[2] = BASE_G;
		word[3] = BASE_U;
		word[4] = BASE_A;
		word[5] = BASE_C;
		
		addToHash(hash, word, 10, wordLength);
		assert(inHash(hash, word, wordLength) >= 0);

	}
	#endif
	
	for (i = 0; i < myStruc->length - wordLength; i++) {
		
		int cont = 1;
		int j;
		for (j = i; j < i + wordLength; j++) {
			if (myStruc->seq[j] == BASE_N) {
				cont = 0;	
			}
		}

		if (cont) {
			int badflag = 0;
			memcpy(word, &myStruc->seq[i], wordLength * sizeof(int));
			word[wordLength] = -1;
			if (inHash(hash, word, wordLength) >= 0) {
				badflag = 1;
			}
			if (invalidComplementInHash(hash, word, i, pairing, wordLength)) {
				badflag = 1;
					
			}
	
			if (badflag == 1) {
				free(word);
				free(hash);
				return 1;
			}
			else {
				addToHash(hash, word, i, wordLength);
			}
		}
	}

  	free(word);	
	free(hash);

	return 0;
}


void setupInitNucs(secStruc *myStruc) {
	int *initNucs = myStruc->initNucs;
	int *initComps = myStruc->initComps;
	int *initMap = myStruc->initMap;
	
	if (initMode == RANDOM_INIT || initMode == SSM_INIT) {
		myStruc->numInitBases = 4;
		initNucs[0] = BASE_A;
		initNucs[1] = BASE_C;
		initNucs[2] = BASE_G;
		initNucs[3] = BASE_U;
		initComps[0] = BASE_U;
		initComps[1] = BASE_G;
		initComps[2] = BASE_C;
		initComps[3] = BASE_A;
		initMap[0] = 3;
		initMap[1] = 2;
		initMap[2] = 1;
		initMap[3] = 0;
	}
	else if (initMode == AU_INIT) {
		myStruc->numInitBases = 2;
		initNucs[0] = BASE_A;
		initNucs[1] = BASE_U;
		initComps[0] = BASE_U;
		initComps[1] = BASE_A;
		initMap[0] = 1;
		initMap[1] = 0; 		
	}
	else if (initMode == CG_INIT) {
		myStruc->numInitBases = 2;
		initNucs[0] = BASE_C;
		initNucs[1] = BASE_G;
		initComps[0] = BASE_G;
		initComps[1] = BASE_C;
		initMap[0] = 1;
		initMap[1] = 0; 		
	}
	else {
		fprintf(stderr,"INVALID INITIALIZATION\n");
		exit(1);
	}	
}


void initSeq(secStruc *myStruc) {
	
	
	
	initHardBases(myStruc);
	
	int wordLength = myStruc->ssmWordLength;
	
	int *initNucs = myStruc->initNucs;
	int *initComps = myStruc->initComps;
	int *initMap = myStruc->initMap;
	
	int numberOfOptions = myStruc->numInitBases;

	#ifdef DEBUG_FULL
		int rnumTracker[numberOfOptions];
		int rnumCount = 0;

		int z;
		for (z = 0; z < numberOfOptions; z++) {
			rnumTracker[z] = 0;
		}
		
		rnumCount = 0;
		
	#endif
	
	int tried[myStruc->length][numberOfOptions];
	
	int i;
	for (i = 0; i < myStruc->length; i++) {
		resetTriedOptions(tried[i], numberOfOptions, myStruc->seq, myStruc->skel, i);
		#ifdef DEBUG
			assert(optionsAtThisPosition(tried[i], numberOfOptions));
		#endif
	}
	
	int cannotInit = 0;
	int furthestBadPosition = 0;
	int failCount = 0;
	
	
	for (i = 0; i < myStruc->length || cannotInit; i++) {
		
		if (i < 0) {
			i = 0;
		}
		int found = 0;
		
		if (myStruc->tempPairing[i] > i) {
			int partner = myStruc->tempPairing[i];
			
			while (!found && optionsAtThisPosition(tried[i], numberOfOptions) && optionsAtThisPosition(tried[partner], numberOfOptions)) {
				int rnum = get_rand_int(numberOfOptions);
				
				#ifdef DEBUG_FULL
					rnumTracker[rnum]++;
					rnumCount++;
				#endif
				
				if (!tried[i][rnum] && !tried[partner][initMap[rnum]]) {
					tried[i][rnum] = 1;
					tried[partner][initMap[rnum]] = 1;
					if (nucArray[myStruc->skel[i]][initNucs[rnum]] && nucArray[myStruc->skel[partner]][initComps[rnum]]) {
						myStruc->seq[i] = initNucs[rnum];
						myStruc->seq[partner] = initNucs[initMap[rnum]];
						#ifdef DEBUG
						assert(initNucs[initMap[rnum]] == initComps[rnum]);
						#endif
						if (!violatesPattern(myStruc->seq, myStruc->skel,0)) {
							if (initMode == SSM_INIT) {
								if (!violatesSSM(myStruc, wordLength)) {
									found = 1;
								}
								else {
									found = 0;	
								}
								
							}
							else {
								found = 1;
							}
						}
					}
				}
			}
		}
		else if (myStruc->tempPairing[i] < 0) {
			while (!found && 
				(optionsAtThisPosition(tried[i], numberOfOptions))) 
			{
				int rnum = get_rand_int(numberOfOptions);
				#ifdef DEBUG_FULL
					rnumTracker[rnum]++;
					rnumCount++;
				#endif
				
				if (!tried[i][rnum]) {
					tried[i][rnum] = 1;
					if (nucArray[myStruc->skel[i]][initNucs[rnum]]) {
						myStruc->seq[i] = initNucs[rnum];
						
						if (!violatesPattern(myStruc->seq, myStruc->skel,0)) {
							if (initMode == SSM_INIT) {
								if (!violatesSSM(myStruc, wordLength)) {
									found = 1;
								}
								else {
									found = 0;	
								}
							}
							
							else {
								found = 1;
							}
						}
					}
				}
			}
		}
		else {
			#ifdef DEBUG
			assert(!violatesPattern(myStruc->seq, myStruc->skel,0));
			assert(!isAFamily(myStruc->seq[i]));
			
			#endif
			if (initMode == SSM_INIT) {
				
				#ifdef DEBUG
				assert(!violatesSSM(myStruc, wordLength));
				
				#endif
			}
			#ifdef DEBUG
			assert(myStruc->tempPairing[i] < i);
			#endif
			found = 1;
		}
		#ifdef DEBUG
		
		assert(found || (!found && !optionsAtThisPosition(tried[i],numberOfOptions)) 
			|| (myStruc->tempPairing[i] >= 0 && !found && !optionsAtThisPosition(tried[myStruc->tempPairing[i]],numberOfOptions)));
		#endif
		
		#ifdef DEBUG
				if (myStruc->tempPairing[i] >=0 ) {	
					assert(optionsAtThisPosition(tried[i], numberOfOptions) == optionsAtThisPosition(tried[myStruc->tempPairing[i]],numberOfOptions));
				}
		#endif
		
		if (!found && (!optionsAtThisPosition(tried[i],numberOfOptions) || 
			((myStruc->tempPairing[i] > i) && !optionsAtThisPosition(tried[myStruc->tempPairing[i]],numberOfOptions)))) 
		{
			if (i == 0) {
				cannotInit = 1;
				fprintf(stderr, "Cannot Initialize This Sequence (01)\n");  
				fprintf(stderr, "Inherent constraints conflict initialization process at position %d.\n", furthestBadPosition);
				exit(1);
			}
			else {

				failCount++;
				
				if (failCount > SSM_RETRY * myStruc->length) {
					cannotInit = 1;
					fprintf(stderr, "Cannot Initialize This Sequence (02)\n");  
					fprintf(stderr, "Inherent constraints conflict initialization process at position %d.\n", furthestBadPosition);
					exit(1);
				}
				else if (failCount > wordLength * SSM_RETRY) {
					int j;
					for (j = 0; j < myStruc->length; j++) {
						resetTriedOptions(tried[j], numberOfOptions, myStruc->seq, myStruc->skel, j);
					}
					i = -1;
				}
				else {
					
					if (i > 0 && myStruc->tempPairing[i - 1] >= 0 && myStruc->tempPairing[i - 1] < i - 1) {
						int target =   myStruc->tempPairing[i - 1];
						while (i >= target) {
							resetTriedOptions(tried[i], numberOfOptions, myStruc->seq, myStruc->skel, i);
							
							if (myStruc->tempPairing[i] >= i) {
								resetTriedOptions(tried[myStruc->tempPairing[i]], numberOfOptions, myStruc->seq, myStruc->skel, myStruc->tempPairing[i] );
							}
							else if (myStruc->tempPairing[i] >= 0 && myStruc->tempPairing[i] < i && myStruc->tempPairing[i] < target) {
								target = myStruc->tempPairing[i];	
							}
							
							i--;	
						}
					}
					else {

						resetTriedOptions(tried[i], numberOfOptions, myStruc->seq, myStruc->skel, i);
						
						if (myStruc->tempPairing[i] > i) {
							resetTriedOptions(tried[myStruc->tempPairing[i]], numberOfOptions, myStruc->seq, myStruc->skel, myStruc->tempPairing[i] );
						}
						#ifdef DEBUG
						else {
							
							assert(myStruc->tempPairing[i] < 0);	
						}
						#endif
						
						#ifdef DEBUG
						if (initMode == SSM_INIT) {
							assert(!violatesSSM(myStruc, wordLength));
							assert(optionsAtThisPosition(tried[i], numberOfOptions));
						}
						#endif

						i-=2;
					}

				}
			}
		}
		else {
			#ifdef DEBUG
			assert(found);
			#endif
			furthestBadPosition = i+1;
		}
	}
	
	#ifdef DEBUG
	assert(!violatesPattern(myStruc->seq, myStruc->skel,0));
	if (initMode == SSM_INIT) {
		assert(!violatesSSM(myStruc, wordLength));
	}
	assert(myStruc->seq[myStruc->length] == -1);
	#endif
	
	memcpy(myStruc->bestSeq, myStruc->seq, (myStruc->length + 1) * sizeof(int));
	
	#ifdef DEBUG_FULL
		printf("RANDOM DIST: %d\n", rnumCount);
		for (z = 0; z < numberOfOptions; z++) {
			printf("%d:%f ",z, (float)((float)rnumTracker[z] / (float)rnumCount));	
		}
		printf("\n");
	#endif
	
	#ifdef DEBUG
		printf("INIT SEQ: \n");
		printSeq(myStruc->seq);
		printSeq(myStruc->skel);
	#endif
}


int ssmWordLength(int length) {
	int wl = (int)ceil(( ((0.5 * log((double)length))) / log(2.0) ) + 2.0);
	#ifdef DEBUG
		assert(wl > 0);
	#endif
	return wl;
}


int seqAlreadyInit(secStruc *myStruc) {
	int i;
	for (i = 0; i < myStruc->length; i++) {
		if (isAFamily(myStruc->seq[i]) || myStruc->seq[i] == 'N') {
			return 0;
		}
	}
	return 1;
}


void initHardBases(secStruc *myStruc) {
	int i;
	for (i = 0; i < myStruc->length; i++) {
		if (!myStruc->skel[i] || !isAFamily(myStruc->skel[i])) {
			myStruc->seq[i] = myStruc->skel[i];		
		}
	}
	#ifdef DEBUG
		assert(!violatesPattern(myStruc->seq, myStruc->skel,0));
	#endif
}


void freeStepList(stepList *theList) {
	if (theList) {
		if (theList->next && theList->pos >= 0) {
			freeStepList((stepList *)theList->next);
		}
		free(theList);
	}
}


void initStepList(stepList *theList) {
	theList->pos = -1;
	theList->from = -1;
	theList->to = -1;
	theList->next = NULL;
}


void clearStepList(stepList *theList) {
	if (theList) {
		freeStepList(theList->next);
		initStepList(theList);
	}
}


void addToStepList(stepList *theList, int pos, int from, int to) {
	stepList *tempList = theList;
	while (tempList->pos != -1) {
		tempList = (stepList *)tempList->next;
	}
	tempList->pos = pos;
	tempList->from = from;
	tempList->to = to;
	tempList->next = (stepList *)calloc(1, sizeof(stepList));
	initStepList(tempList->next);
}


int inStepList(stepList *theList, int pos, int to) {
	if (!theList || theList->pos == -1) {
		return 0;
	}
	stepList *tempList = theList;
	while (tempList->pos != -1) {
		if (tempList->pos == pos && (tempList->to == to)) {
			return 1;
		}
		tempList = (stepList *)tempList->next;
	}
	return 0;
}


void resetSeqPosToBest(secStruc *myStruc, int pos) {
	myStruc->seq[pos] = myStruc->bestSeq[pos];
}


void computeMFE(secStruc *myStruc) {
	struct timeval start_timeval;
	struct timeval end_timeval;
	struct timeval elapsed_timeval;
		
	gettimeofday(&start_timeval, 0);
	createExpandedSeq(myStruc->expandedSeq, myStruc->seq, myStruc->length, myStruc->numStrands, myStruc->strandBreaks);

	mfeFull(myStruc->expandedSeq, myStruc->expandedLength, myStruc->mfePairs, 3, DNARNACOUNT, DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_SALT);
	
	myStruc->optVal = (DBL_TYPE)getMFEDifference(myStruc);
	myStruc->mfeDifference = (int)myStruc->optVal;
	
	myStruc->timesOptimized++;
	gettimeofday(&end_timeval, 0);
  timeDifference (&elapsed_timeval, &end_timeval, &start_timeval);
	myStruc->mfeTimeSpent = (DBL_TYPE)elapsed_timeval.tv_sec + 1e-6 * (DBL_TYPE)elapsed_timeval.tv_usec;
}


void computeP(secStruc *myStruc) {
  struct timeval start_timeval;
	struct timeval end_timeval;
	struct timeval elapsed_timeval;
		
	gettimeofday(&start_timeval, 0);
	createExpandedSeq(myStruc->expandedSeq, myStruc->seq, myStruc->length, myStruc->numStrands, myStruc->strandBreaks);

	myStruc->optVal = psStarPairsOrParensFull( myStruc->tempPairing, NULL, myStruc->expandedSeq, 3, DNARNACOUNT, DANGLETYPE, TEMP_K - ZERO_C_IN_KELVIN, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_SALT);
	myStruc->pVal = myStruc->optVal;
	
	myStruc->timesOptimized++;
	gettimeofday(&end_timeval, 0);
  timeDifference (&elapsed_timeval, &end_timeval, &start_timeval);
	myStruc->pTimeSpent = (DBL_TYPE)elapsed_timeval.tv_sec + 1e-6 * (DBL_TYPE)elapsed_timeval.tv_usec;
}


#ifdef OUTPUT_TREE_SNAPSHOT
void outputTreeSnapshot(secStruc *myStruc, FILE *out, int k, int l) {

	if (!myStruc->parent) {
		assert(k == 1);
		assert(l == 1);
		assert(out == NULL);
		char *filename = calloc(1000, sizeof(char));
		sprintf(filename, "%s_%d.dot", myStruc->foldFile, myStruc->timesOptimized);
		out = fopen(filename,"w");
		fprintf(out, "graph \"\" {\n");
		
	}
	
	assert(out != NULL);
	if (hasChildren(myStruc)) {
		DBL_TYPE nativeSum = getNativeN(myStruc->leftChild) + getNativeN(myStruc->rightChild);
		fprintf(out, "node%d_%d [label=\"%Lf : %Lf : %d\"];\n", k,l,myStruc->optVal, nativeSum, myStruc->length);
	}
	else {
		fprintf(out, "node%d_%d [label=\"%Lf : %d\"];\n", k,l,myStruc->optVal, myStruc->length);
	}	
	
	
	
	if (myStruc->leftChild) {
		outputTreeSnapshot(myStruc->leftChild, out,  k+1, 2*l - 1);
		outputTreeSnapshot(myStruc->rightChild, out, k+1, 2*l);
		fprintf(out,"node%d_%d -- node%d_%d;\n", k,l,k+1, 2*l - 1);
		fprintf(out,"node%d_%d -- node%d_%d;\n", k,l,k+1, 2*l);

	}
		
	if (!myStruc->parent) {
		fprintf(out, "}");
	}
	
	
}
#endif


void computeNSstar(secStruc *myStruc) {
  struct timeval start_timeval;
	struct timeval end_timeval;
	struct timeval elapsed_timeval;
		
	gettimeofday(&start_timeval, 0);
  
  
	#ifdef DEBUG
		assert(myStruc->bestOptVal >= 0);
		assert(!violatesPattern(myStruc->seq, myStruc->skel,0));
		assert(!violatesSkeleton(myStruc->length, myStruc->seq, myStruc->skel));
		if (myStruc->numStrands > 1) {
        	assert(myStruc->expandedLength == myStruc->length + myStruc->numStrands - 1); 
    }
		assert(myStruc->numStrands >= 1);
	#endif
	
	pairPr = myStruc->myPairPr; 
	
	myStruc->nsCalls++;
	
	createExpandedSeq(myStruc->expandedSeq, myStruc->seq, myStruc->length, myStruc->numStrands, myStruc->strandBreaks);
	#ifdef DEBUG
		assert(getSeqLength(myStruc->seq) == myStruc->length);
		assert(getSeqLength(myStruc->expandedSeq) == getSeqLength(myStruc->seq) + myStruc->numStrands - 1);
		assert(getSeqLength(myStruc->seq) == getSeqLength(myStruc->bestSeq));
    assert(myStruc->expandedLength == getSeqLength(myStruc->expandedSeq));

	#endif

	if (myStruc->parent) {
	  myStruc->optVal = nsStarPairsOrParensFull( myStruc->length,  myStruc->expandedSeq, myStruc->tempPairing, NULL, 3, DNARNACOUNT, DANGLETYPE,  TEMP_K - ZERO_C_IN_KELVIN, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_SALT);
	}
	else {
		myStruc->optVal = nsStarPairsOrParensCorrected( myStruc->length, myStruc->numStrands, myStruc->expandedSeq, myStruc->tempPairing, NULL, myStruc->fakeBases, &myStruc->pfVal, 3, DNARNACOUNT, DANGLETYPE, 
		       TEMP_K - ZERO_C_IN_KELVIN, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_SALT);
	}
	
	myStruc->nVal = myStruc->optVal;
	
	myStruc->timesOptimized++;
	#ifdef OUTPUT_TREE_SNAPSHOT
		if (!myStruc->parent) {
			outputTreeSnapshot(myStruc, NULL, 1, 1);
		}
	#endif
	
	
	#ifdef DEBUG
	if (hasChildren(myStruc)) {
			printSeq(myStruc->expandedSeq);
			printf("VAL: %Lf:%d\n", myStruc->optVal, myStruc->relevantLength);
			printf("SUM: %Lf,:%Lf\n", myStruc->leftChild->optVal, myStruc->rightChild->optVal);
	}
	#endif
	gettimeofday(&end_timeval, 0);
  timeDifference (&elapsed_timeval, &end_timeval, &start_timeval);
	myStruc->pfTimeSpent = (DBL_TYPE)elapsed_timeval.tv_sec + 1e-6 * (DBL_TYPE)elapsed_timeval.tv_usec;
	myStruc->checked = 1;
  
  #ifdef EARLY_TERMINATE
  if (!myStruc->parent) {
    FILE *timerOut = fopen(timerFilename,"w");
    struct timeval timer_timeval;
    timeDifference (&timer_timeval, &end_timeval, global_timeval);
    DBL_TYPE timerSpent = (DBL_TYPE)elapsed_timeval.tv_sec + 1e-6 * (DBL_TYPE)elapsed_timeval.tv_usec;
    fprintf(timerOut, "%Lf\n", timerSpent);
    fclose(timerOut);
  }
  #endif
  
  
}


int getMFEDifference(secStruc *myStruc) {
	myStruc->mfeDifference = 0;
	int j;		
	for (j = 0; j < myStruc->length; j++) {
		myStruc->myProbs[j] = MFE_EPS;
		if (myStruc->tempPairing[j] != myStruc->mfePairs[j]) {
			myStruc->mfeDifference++;
			myStruc->myProbs[j] += 1.0;
		}
	}
	return myStruc->mfeDifference;
}


void calculateObjective(secStruc *myStruc) {
  


	if (designMode == N_OPTIMIZATION) {
		computeNSstar(myStruc);
	}
	else if (designMode == MFE_OPTIMIZATION) {
		computeMFE(myStruc);
	}
	else if (designMode == P_OPTIMIZATION) {
		computeP(myStruc);
	}

}


int mutationProhibited(secStruc *myStruc, mutation *candidateMutation, int *currentPatternViolations, int acceptablePatternViolations) {
	if (!nucArray[myStruc->skel[candidateMutation->posA]][candidateMutation->nucA]) {
		#ifdef DEBUG
			if (candidateMutation->posB >= 0) {
				assert(!nucArray[myStruc->skel[candidateMutation->posB]][candidateMutation->nucB]);	
			}
		#endif
		return 1;
	}
	#ifdef DEBUG
	else {
		if (candidateMutation->posB >= 0) {
			assert(nucArray[myStruc->skel[candidateMutation->posB]][candidateMutation->nucB]);	
		}	
	}
	#endif

	int *tempSeq = (int *)calloc(myStruc->length + 1, sizeof(int));
	memcpy(tempSeq, myStruc->seq, sizeof(int) * (myStruc->length + 1));
	
	tempSeq[candidateMutation->posA] = candidateMutation->nucA;
	if (candidateMutation->posB >= 0) {
		tempSeq[candidateMutation->posB] = candidateMutation->nucB;
	}
	#ifdef DEBUG
		assert(intArrayLen(tempSeq) == myStruc->length);
	#endif
	
	*currentPatternViolations = violatesPattern(tempSeq, myStruc->skel,acceptablePatternViolations);
	free(tempSeq);
	if (*currentPatternViolations == 0) {
		return 0;	
	}
	else if (*currentPatternViolations >= acceptablePatternViolations) {
		return 1;	
	}
	
	return 0;
}


void applyMutation(secStruc *myStruc, mutation *candidateMutation) {	

	myStruc->seq[candidateMutation->posA] = candidateMutation->nucA;

	if (candidateMutation->posB >= 0) {	
		myStruc->seq[candidateMutation->posB] = candidateMutation->nucB;
	}
}


int mutationUnfavorable(mutation *candidateMutation, stepList *list) {
	int retVal = 0;
	if (inStepList(list, candidateMutation->posA, candidateMutation->nucA)) {
			retVal = 1;
	}
	
	#ifdef DEBUG
		if (candidateMutation->posB >= 0) {
			assert(candidateMutation->nucB >= 0);
			int pairStep = inStepList(list, candidateMutation->posB, candidateMutation->nucB);
			if (retVal) {
				assert(pairStep);
			}
			else {
				assert(!pairStep);	
			}
		}
		else {
			assert(candidateMutation->posB == -1);
			assert(candidateMutation->nucB == -1);
		}
	#endif
	
	return retVal;
}


int favorableMove(DBL_TYPE newVal, DBL_TYPE oldVal, int optType) {
  if (optType != P_OPTIMIZATION) {
   return (newVal < oldVal);    
  }
  else {
    return (newVal > oldVal);
  }
}


#ifdef FREEZEBASES
void optimizeLeaf(secStruc *myStruc, int freeze) {
#else
void optimizeLeaf(secStruc *myStruc) {
#endif
  
	#ifdef DEBUG
		assert(intArrayLen(myStruc->seq) == myStruc->length);
		assert(!violatesPattern(myStruc->seq, myStruc->skel,0));
		printf("Optimizing %s\n",myStruc->tempFolding);
	#endif
		
	int successes[7];
	int tries[7];
	int lastone;
	
	DBL_TYPE pACC = 0.0;
	
	if (designMode == MFE_OPTIMIZATION) {
		pACC = B_ACCEPT;
	}
	
	int length = myStruc->length;
	
	myStruc->bestOptVal = (DBL_TYPE)myStruc->length;

	memcpy(myStruc->bestSeq, myStruc->seq, sizeof(int)*myStruc->length);
	
	int z;
	for (z = 0; z < 7; z++) {
		successes[z] = 0;
		tries[z] = 0;
	}
	lastone = 0;
	
		
	DBL_TYPE n_m_unfavorable = M_UNFAVORABLE * (DBL_TYPE)length;
	
	stepList theList;
	initStepList(&theList);
	
	calculateObjective(myStruc);

	calculateProbs(myStruc);
	setNewBestValue(myStruc);
	
	#ifdef DEBUG
	printf("Starting val: %Lf\n",myStruc->bestOptVal);
	#endif

	int m_unfavor = 0;
	
	int m_eval = 0;

	DBL_TYPE rnum = 100.0;
	while (!satisfiedLeafObjective(myStruc->bestOptVal,myStruc->targetObjective,designMode) && m_unfavor < n_m_unfavorable && (DBL_TYPE)m_eval < M_EVAL) { 

		
		mutation candidateMutation;
		
		if (bypassGuidance) {
			mutateRandomBaseOrPair(myStruc, &candidateMutation);
			lastone = 2;
		}
		else {
			mutateRandomConflict(myStruc, &candidateMutation); 
			lastone = 1;
		}

		int currentViolations = 0;

    #ifdef FREEZEBASES
    if (freeze && (!myStruc->purelyNativeBases[candidateMutation.posA] || !myStruc->purelyNativeBases[candidateMutation.posB])) {
      // do nothing
    }
 		else if ((mutationUnfavorable(&candidateMutation, &theList) || mutationProhibited(myStruc, &candidateMutation, &currentViolations, 0))) {
    #else
		if ((mutationUnfavorable(&candidateMutation, &theList) || mutationProhibited(myStruc, &candidateMutation, &currentViolations, 0))) {
    #endif
			m_unfavor++;	
		}
		else {
			applyMutation(myStruc, &candidateMutation);	
			calculateObjective(myStruc);
			tries[lastone]++;
			rnum = (DBL_TYPE)genrand_real1();
			
			if (rnum < pACC || favorableMove(myStruc->optVal,myStruc->bestOptVal,designMode)) {
				m_unfavor = 0;
				
				clearStepList(&theList);
				successes[lastone]++;
				calculateProbs(myStruc);
				setNewBestValue(myStruc);

				#ifdef DEBUG
					printf("SUCCESS: %d| %Lf| %d %d %d %d | %d %d %d %d | %d \n", myStruc->length, myStruc->optVal, successes[1], successes[2], successes[5], successes[6], tries[1], tries[2],  tries[5], tries[6],  m_unfavor);
				#endif

			}
			else {
        #ifdef DEBUG
        if (designMode == P_OPTIMIZATION) {
          printf("FAILURE: %.14Le\n", myStruc->optVal);
        }
        #endif
        int posA = candidateMutation.posA;
        int posB = candidateMutation.posB;
        addToStepList(&theList, posA, myStruc->seq[posA], candidateMutation.nucA);
        if (posB >= 0) {
          addToStepList(&theList, posB, myStruc->seq[posB], candidateMutation.nucB);
        }

				m_unfavor++;
			}

			revertToOldValue(myStruc);


		}
		rnum = (DBL_TYPE)genrand_real1(); // keeping this for random seed backwards compatibility
		if (designMode == MFE_OPTIMIZATION) {
			m_eval++;
		}

	}
	revertToOldValue(myStruc);
	
	#ifdef DEBUG
		assertValueInSync(myStruc);
		assertInSync(myStruc);
	#endif
	
	freeStepList(theList.next);
}


void mutationAtPos(secStruc *myStruc, mutation *candidateMutation, int pos) {
	candidateMutation->posA = pos;
#ifdef DEBUG
	assert(pos < myStruc->length);
	assert(myStruc->seq[pos] >= 0);	
#endif
	candidateMutation->nucA = getNucNot(myStruc->seq[pos]);
	

	if (myStruc->tempPairing[pos] >= 0) {
		candidateMutation->posB = myStruc->tempPairing[pos];
		candidateMutation->nucB = newNucsComps[candidateMutation->nucA];
	}
	else {
		candidateMutation->posB = -1;
		candidateMutation->nucB = -1;
	}
}

void calculateProbs(secStruc *myStruc) {

	DBL_TYPE *incorrectProbs = myStruc->myProbs;

	int length = myStruc->length;
	int i;
	DBL_TYPE probSum = 0.0;
	for (i = 0; i < length; i++) {
		
		if (designMode == N_OPTIMIZATION) {
			DBL_TYPE probUnpaired = myStruc->myPairPr[i*(myStruc->length+1) + myStruc->length];
			if (myStruc->tempPairing[i] >= 0) {
				incorrectProbs[i] = 1 - myStruc->myPairPr[i*(myStruc->length+1)+myStruc->tempPairing[i]];
			}
			else if (myStruc->tempPairing[i] < 0) {
				incorrectProbs[i] = 1 - probUnpaired;
			}
		}
		#ifdef DEBUG
    if (designMode == N_OPTIMIZATION) {
      assert(incorrectProbs[i] <= 1.01);
    }
    else if (designMode == MFE_OPTIMIZATION) {
      assert(incorrectProbs[i] <= 1.11);
 
    }
		#endif
		probSum += incorrectProbs[i];

	}
	for (i = 0; i < length; i++) {
		incorrectProbs[i] = incorrectProbs[i] / probSum;
	}
}


void mutateRandomConflict(secStruc *myStruc, mutation *candidateMutation) {
	int length = myStruc->length;
	
	
	DBL_TYPE *incorrectProbs = myStruc->myProbs;
	
	
	DBL_TYPE randNum = genrand_real2();
	
	int theOne;

	int i;
	DBL_TYPE probSum = 0.0;
	int cont = 1;
	for (i = 0; i < length && cont; i++) {
		if (randNum >= probSum && randNum < probSum + incorrectProbs[i]) {
			theOne = i;
			cont = 0;
		}
		probSum += incorrectProbs[i];
		#ifdef DEBUG
			assert(probSum <= 1.01);
		#endif
		
	}
	mutationAtPos(myStruc, candidateMutation, theOne);
}


void mutateRandomBaseOrPair(secStruc *myStruc, mutation *candidateMutation) {
	int pos = get_rand_int(myStruc->length);
	mutationAtPos(myStruc, candidateMutation, pos);
}


void initConflictList(conflictList *theList) {
	theList->pos = -1;
	theList->base = -1;
	theList->next = NULL;
}


void freeConflictList(conflictList *theList) {
	if (theList) {
		if (theList->next != NULL) {
			freeConflictList(theList->next);
		}
		free(theList);
	}
}


void addToConflictList(conflictList *theList, int pos, int base) {
	conflictList *tempList = theList;
	while (tempList->pos != -1) {
		tempList = (conflictList *)tempList->next;
	}
	
	tempList->pos = pos;
	tempList->base = base;
	tempList->next = (conflictList *)calloc(1, sizeof(conflictList));
	initConflictList((conflictList *)tempList->next);
}


int inConflictList(conflictList *theList, int pos, int base) {
	int retVal = 0;
	if (!theList) {
		return 0;
	}
	conflictList *tempList = theList;
	while (tempList->pos != -1 && !retVal) {
		if (tempList->pos == pos && tempList->base == base) {
			
			retVal = 1;
		}
		tempList = (conflictList *)tempList->next;
	}
	
	return retVal;
}


int conflictListSize(conflictList *theList) {
	int retVal = 0;
	
	conflictList *tempList = theList;
	while (tempList->pos != -1) {
		retVal++;
		tempList = (conflictList *)tempList->next;
	}
	return retVal;
}


int inEitherList(stepList *theList, conflictList *conflicts, int pos, int base, int *calcsaves) {
	if (inStepList(theList, pos, base)) {
		*calcsaves = *calcsaves + 1;
		return 1;
	}
	else return inConflictList(conflicts, pos, base);
}


void convertConflictArrayToList(secStruc *myStruc, conflictList *list) {
	initConflictList(list);
	int z;
	for (z = 0; z < myStruc->length; z++) {
		if (myStruc->conflictArray[z] == 1) {
			addToConflictList(list, z,myStruc->seq[z]);
		}
	}
}


void mapStrandBreaksDown(secStruc *parent, secStruc *child) {
	int *mapping;
	child->strandBreaks = (int *)calloc(parent->numStrands, sizeof(int));
	if (child->amILeft) {
		mapping = parent->leftChildMapping;

	}
	else {
		mapping = parent->rightChildMapping;
	}
	child->numStrands = 1;
	int i;
	int breakCounter = 0;
	for (i = 0; i < parent->numStrands - 1; i++) {
		int potentialMap = mapping[parent->strandBreaks[i]];
		
		#ifdef DEBUG
			assert((potentialMap < child->length && potentialMap >= 0) || potentialMap == -99);
		#endif
		if (potentialMap >= 0) {
			child->strandBreaks[breakCounter] = potentialMap;
			child->numStrands++;
			breakCounter++;
		}
	}
	
	if (child->amILeft) {
		int j;
		int breakPoint = child->loopStart;
		int cont = 1;
		for (j = 0; cont && j < child->numStrands - 1; j++) {
			if (child->strandBreaks[j] > breakPoint) {
				j--;
				cont = 0;
			}
		}
		child->numStrands++;
		int lastVal = breakPoint;
		for (; j < child->numStrands - 1; j++) {
			int currentVal = child->strandBreaks[j];
			child->strandBreaks[j] = lastVal;
			lastVal = currentVal;
		}
	}
	
  	child->expandedLength = child->length + child->numStrands - 1;
	#ifdef DEBUG
		assert(child->numStrands > 0);
	#endif
}

#ifdef DEBUG
void assertWatsonCrick(secStruc *myStruc) {
	if (myStruc->amILeft) {
		printf("LEFT:\n");	
	}
	printSeq(myStruc->seq);
	int i = 0;
	for (i = 0; i < myStruc->length; i++) {
		if (myStruc->tempPairing[i] >= 0) {
			
			assert(myStruc->seq[i] == newNucsComps[myStruc->seq[myStruc->tempPairing[i]]]);	
		}
	}
	
}

#endif

void mapSequenceDown(secStruc *parent, secStruc *child, int skel) {
	int *parentSequence;
	int *childSequence;
	if (skel) {
		parentSequence = parent->skel;
		childSequence = child->skel;
	}
	else {
		parentSequence = parent->seq;
		childSequence = child->seq;
	}

	int *mapping;
	if (child->amILeft) {
		mapping = parent->leftChildMapping;
	}   
	else {
		mapping = parent->rightChildMapping;
	}
	
	int i;
	for (i = 0; i < parent->length; i++) {
		if (mapping[i] >= 0) {
			childSequence[mapping[i]] = parentSequence[i];
		}
	}
	#ifdef DEBUG
		assertWatsonCrick(parent);
		assertWatsonCrick(child);
	#endif
	
}


void mapSequenceUp(int *parentSeq, int *childSequence, secStruc *child, int left) {
	int *mapping = child->parentMapping;
	int i;
	
	int offset = hMin;
	
	for (i = 0; i < child->length; i++) {
		if (left && (i < child->loopStart - offset || i >= child->loopStart + offset)) {
			#ifdef DEBUG
				assert(mapping[i] >= 0);
			#endif
			parentSeq[mapping[i]] = childSequence[i];
		}
		else if (!left && i > (offset - 1) && i < (child->length - offset)) {
			parentSeq[mapping[i]] = childSequence[i];
		}
	}
	
	#ifdef DEBUG_FULL
		printf("CHILD: %d\n", left);
		printSeq(childSequence);
		printf("%s\n",child->tempFolding);
		printf("PARENT: \n");
		printSeq(parentSeq);
		printf("%s\n",child->parent->tempFolding);
	#endif

}


void decompose(secStruc *parent) {
	#ifdef DEBUG
		assert(hasChildren(parent));
	#endif


  
	mapSequenceDown(parent, parent->leftChild, 0);
	mapSequenceDown(parent, parent->rightChild, 0);
	
	fixPatternViolations(parent->leftChild);
	fixPatternViolations(parent->rightChild);
	
  #ifdef LEAFCORRECTION
  overwriteChildrenLeafArrays(parent);
  #endif
  
}


void merge(secStruc *left, secStruc *right) {
	#ifdef DEBUG
		assert(left->parent == right->parent);
		assert(left->amILeft);
		assert(!right->amILeft);
		
		
	#endif	
	
	secStruc *parent = left->parent;
	mapSequenceUp(parent->seq, left->seq, left, 1);
	mapSequenceUp(parent->seq, right->seq, right, 0);
	
	#ifdef DEBUG
		printf("MERGED: %s\n", parent->tempFolding);
		printSeq(parent->seq);
		printf("%s\n",left->tempFolding);
		printf("%Lf: ", left->optVal);
		printSeq(left->seq);
		printf("%s\n", right->tempFolding);
		printf("%Lf: ", right->optVal);
		printSeq(right->seq);


	#endif
		
	parent->checked = 0;
	/* did emergent pattern violations occur */
	fixPatternViolations(parent);	
	calculateObjective(parent);

}


#ifdef DEBUG
void assertValueInSync(secStruc *myStruc) {
	if (mReopt >= 0) {
		assertInSync(myStruc);
	}
}
#endif


#ifdef DEBUG
void printTree(secStruc *myStruc, int level, FILE *treeOut) {
	
    
	
  int i;
  for (i = 0; i < level; i++) {
    printf("\t");
  }
  
  fprintf(treeOut, "node_%d_%d [label = \"%d\"];", myStruc->k, myStruc->l, myStruc->length);
  printf("%s\n",myStruc->tempFolding);

  if (hasChildren(myStruc)) {
    printTree(myStruc->leftChild, level+1, treeOut);
    fprintf(treeOut, "node_%d_%d -> node_%d_%d;", myStruc->k, myStruc->l, myStruc->leftChild->k, myStruc->leftChild->l); 
    printTree(myStruc->rightChild, level+1, treeOut);
    fprintf(treeOut, "node_%d_%d -> node_%d_%d;", myStruc->k, myStruc->l, myStruc->rightChild->k, myStruc->rightChild->l); 
  }
}
#endif


DBL_TYPE calculateLeafSum(secStruc *myStruc) {
  
  DBL_TYPE *nValues = (DBL_TYPE*)calloc(myStruc->length,sizeof(DBL_TYPE));

  getLeafNValues(myStruc, nValues);
  int i;
  DBL_TYPE sum = 0.0;
  for (i = 0; i < myStruc->length; i++) {
     sum += nValues[i];
  }
  
  return sum;
}


void getLeafNValues(secStruc *myStruc, DBL_TYPE *nValues) {
  if (hasChildren(myStruc)) {
    if (mReopt != 0) {
      decompose(myStruc);	
    }
    
    
    DBL_TYPE *leftValues = (DBL_TYPE*)calloc(myStruc->leftChild->length,sizeof(DBL_TYPE));
    DBL_TYPE *rightValues = (DBL_TYPE*)calloc(myStruc->rightChild->length,sizeof(DBL_TYPE));
    
    getLeafNValues(myStruc->leftChild, leftValues); 
    getLeafNValues(myStruc->rightChild, rightValues); 
    
    
#ifdef DEBUG
    int count = 0;
#endif
    int i;
    for (i = 0; i < myStruc->leftChild->length; i++) {
      if (myStruc->leftChild->parentMapping[i] >= 0) {
        nValues[myStruc->leftChild->parentMapping[i]] = leftValues[i];
#ifdef DEBUG
        count++;
#endif
      }
      #ifdef DEBUG
      else {
        assert(myStruc->leftChild->purelyNativeBases[i] == 0); 
      }
      #endif
    }
    for (i = 0; i < myStruc->rightChild->length; i++) {
      if (myStruc->rightChild->parentMapping[i] >= 0) {
        nValues[myStruc->rightChild->parentMapping[i]] = rightValues[i];
#ifdef DEBUG
        count++;
#endif
      }
      #ifdef DEBUG
      else {
        assert(myStruc->rightChild->purelyNativeBases[i] == 0); 
      }
      #endif
    }
#ifdef DEBUG
    assert(count == myStruc->length);
#endif
    
    return;
    
  }
  else {
    if (mReopt != 0) {
      calculateObjective(myStruc);
    }
    getNValues(myStruc, nValues);
    return;
  }
}


void optimizeNode(secStruc *myStruc, DBL_TYPE *nValues) {
	myStruc->checked = 0;
	#ifdef DEBUG
		assert(hasChildren(myStruc));
		assertInSync(myStruc);
	#endif

	secStruc *leftChild = myStruc->leftChild;
	secStruc *rightChild = myStruc->rightChild;
	
	decompose(myStruc);
	
	int i;
	int *leftMap = myStruc->leftChildMapping;
	int *rightMap = myStruc->rightChildMapping;
	DBL_TYPE *leftValues = (DBL_TYPE*)calloc(leftChild->length, sizeof(DBL_TYPE));
	DBL_TYPE *rightValues = (DBL_TYPE*)calloc(rightChild->length, sizeof(DBL_TYPE));


	DBL_TYPE rightTotal = 0.0;
	DBL_TYPE leftTotal = 0.0;
	for (i = 0; i < myStruc->length; i++) {	
		if (rightMap[i] >= 0 && rightChild->parentMapping[rightMap[i]] >= 0) {
			rightValues[rightMap[i]] = nValues[i];
			rightTotal += nValues[i];
		}
		if (leftMap[i] >= 0 && leftChild->parentMapping[leftMap[i]] >= 0) {
			leftValues[leftMap[i]] = nValues[i];
			leftTotal += nValues[i];
		}
	}

	DBL_TYPE rnum = (DBL_TYPE)genrand_real1();
	#ifdef DEBUG
	printf("LEFT:RIGHT: %Lf\n", leftTotal / (leftTotal + rightTotal));
	#endif
	
	if (rnum < leftTotal / (leftTotal + rightTotal)) {
		designSeq(leftChild, leftValues);
	}
	else {
		designSeq(rightChild, rightValues);
	}

	merge(leftChild, rightChild);

  
	#ifdef DEBUG
		assertInSync(myStruc);

	#endif
	free(leftValues);
	free(rightValues);
}


int hasChildren(secStruc *parent) {
	if (parent->leftChild) {
		#ifdef DEBUG
			assert(parent->rightChild);
			assert(parent->leftChild->parent == parent);
			assert(parent->rightChild->parent == parent);
			assert(parent->leftChild->amILeft);
			assert(!(parent->rightChild->amILeft));
		#endif
		return 1;
	}
	else {
		#ifdef DEBUG
			assert(!parent->rightChild);
		#endif
		return 0;	
	}
}

DBL_TYPE getNativeN(secStruc *myStruc) {
  #ifdef DEBUG 
  	assert(floatsEqual(myStruc->optVal,myStruc->bestOptVal));
  #endif
  int i;
  DBL_TYPE value = 0.0;
  int pair;
  int dummyCount = 0;

  int seqlength = myStruc->length;
  int *structPairs = myStruc->tempPairing;
  DBL_TYPE *pr = myStruc->myPairPr;
  
  value = 0;
  for( i = 0; i< seqlength; i++) {
	  if (myStruc->parentMapping[i] >= 0) {
		  pair = structPairs[i];
		  if( pair < 0) {
			  value += pr[ i*(seqlength+1) + seqlength];
		  }
		  else {
			  value += pr[ i*(seqlength+1) + pair];
		  }
	  }
	  else {
	    dummyCount++;
	  }
  }

  #ifdef DEBUG
  	printf("%d:D: %d:%d\n", myStruc->amILeft, dummyCount, hMin * 2);
	printf("LCVAL: %Lf\n", value);
  #endif
  
  return (DBL_TYPE)myStruc->length - (value + (DBL_TYPE)dummyCount);
	
}


DBL_TYPE getChildContributingN(secStruc *myStruc) {

  
  
	//DBL_TYPE *incorrectProbs = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
  DBL_TYPE incorrectSum = 0.0;
	int length = myStruc->length;
	int i;

  #ifdef DEBUG
  int count = 0;
  #endif
	for (i = 0; i < length; i++) {
		
		if (myStruc->parentMapping[i] >= 0) {
      #ifdef DEBUG
      count++;
      #endif
      if (designMode == N_OPTIMIZATION) {
        DBL_TYPE probUnpaired = myStruc->myPairPr[i*(myStruc->length+1) + myStruc->length];
        if (myStruc->tempPairing[i] >= 0) {
          incorrectSum += 1.0 - myStruc->myPairPr[i*(myStruc->length+1)+myStruc->tempPairing[i]];
        }
        else if (myStruc->tempPairing[i] < 0) {
          incorrectSum += 1.0 - probUnpaired;
        }
      }
      else {
        if (myStruc->tempPairing[i] >= 0 && myStruc->mfePairs[i] != myStruc->tempPairing[i]) {
          incorrectSum += 1.0;
        }
        else if (myStruc->tempPairing[i] < 0 && myStruc->mfePairs[i] >= 0) {
          incorrectSum += 1.0;
        }
        
      }
		}
	}
  #ifdef DEBUG
  assert(count == length - 2 * hMin);
  #endif
  return incorrectSum;
}

int satisfiedLeafObjective(DBL_TYPE bestValue, DBL_TYPE targetObjective, int optType) {
	
#ifdef DEBUG
  printf("TARGET: %Lf\n",targetObjective);
#endif
  if (optType != P_OPTIMIZATION && bestValue > targetObjective) {
    return 0;	
  }
  else if (optType == P_OPTIMIZATION && bestValue < targetObjective) {
    return 0; 
  }
  else if (optType == P_OPTIMIZATION && bestValue >= targetObjective) {
    return 1; 
  }

	return 1;
	
}


int satisfiedObjective(secStruc *myStruc) {
	
	if (designMode != P_OPTIMIZATION) {// && myStruc->optVal <= myStruc->targetObjective) {
		//return 1;	
	}
  else if (designMode == P_OPTIMIZATION && myStruc->optVal >= myStruc->targetObjective) {
    return 1; 
  }
  else if (designMode == P_OPTIMIZATION && myStruc->optVal < myStruc->targetObjective) {
    return 0; 
  }

	secStruc *left = myStruc->leftChild;
	//DBL_TYPE maxLeft;// = getChildContributingN(left);

  secStruc *right = myStruc->rightChild;
	//DBL_TYPE maxRight;// = getChildContributingN(right);

	
	DBL_TYPE maxLeft = myStruc->leftChild->leafNative[myStruc->leftChild->id];//getChildContributingN(left);
  if (maxLeft < left->relevantLength * nRatio) {
   maxLeft = left->relevantLength * nRatio;  
  }
  
  
  assertNativeChecksum(myStruc->leftChild);
	DBL_TYPE maxRight = myStruc->rightChild->leafNative[myStruc->rightChild->id];//getChildContributingN(right);
  if (maxRight < right->relevantLength * nRatio) {
    maxRight = right->relevantLength * nRatio;  
  }
  assertNativeChecksum(myStruc->rightChild);

	#ifdef DEBUG
		assert(floatsEqual(left->optVal,left->bestOptVal));
		assert(floatsEqual(right->optVal,right->bestOptVal));
	#endif
	
	DBL_TYPE childSum = maxLeft + maxRight;
	if (myStruc->optVal <= (childSum)) {
		return 1;
	}
	return 0;
	
}

void getNValues(secStruc *myStruc, DBL_TYPE *nValues) {
	
	int j;
	if (designMode == MFE_OPTIMIZATION) {
		for (j = 0; j < myStruc->length; j++) {
			nValues[j] = MFE_EPS;
			if (myStruc->tempPairing[j] != myStruc->mfePairs[j]) {
				nValues[j] += 1.0;
			}
		}
	}
	else {
		for (j = 0; j < myStruc->length; j++) {
			nValues[j] = 0.0;
	
			if (myStruc->tempPairing[j] >= 0) {
				nValues[j] = 1.0 - myStruc->myPairPr[j*(myStruc->length+1)+myStruc->tempPairing[j]];
			}
			else if (myStruc->tempPairing[j] < 0) {
				nValues[j] = 1.0 - myStruc->myPairPr[j*(myStruc->length+1) + myStruc->length];
			}
		
		}
	}	
}

void designSeq(secStruc *myStruc, DBL_TYPE *nValues) {

	/* have we initialized the best seq */
	if (myStruc->bestSeq[0] == BASE_N) {
		setNewBestValue(myStruc);
	}
	
  #ifdef LEAFCORRECTION
  myStruc->modified = 1;
  myStruc->leafNativeCheck[myStruc->id] = 0;
  #endif
  
	#ifdef DEBUG
		assert(myStruc->bestOptVal >= 0);
	#endif
		
  #ifdef DEBUG
  if (bypassGuidance && nValues != NULL) {
    int i = 0;
    for (i = 0; i < myStruc->length; i++) {
      assert(floatsEqual(nValues[i],1.0) || floatsEqual(nValues[i],0.0)); 
    }
  }
  #endif
  
	secStruc *leftChild = myStruc->leftChild;
	secStruc *rightChild = myStruc->rightChild;
	
	if (hasChildren(myStruc)) {

		if (nValues == NULL) {
			decompose(myStruc);
			designSeq(leftChild, NULL);
			designSeq(rightChild, NULL);
   		merge(leftChild, rightChild);
  
		}
		else {
			optimizeNode(myStruc, nValues);
		}

		int m = 0;
		

    
		setNewBestValue(myStruc);
		revertToOldValue(myStruc);
    decompose(myStruc);
    #ifdef LEAFCORRECTION
    assertNativeChecksum(myStruc);
    #endif
    
		while (m < mReopt && !satisfiedObjective(myStruc)) {
			
			DBL_TYPE *myValues = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));
      if (bypassGuidance) {
        int i;
        for (i = 0; i < myStruc->length; i++) {
          myValues[i] = 1.0;  
        }
      }
      else {
        getNValues(myStruc, myValues);
      }
      optimizeNode(myStruc, myValues);

      if (favorableMove(myStruc->optVal,myStruc->bestOptVal,designMode)) {
        setNewBestValue(myStruc);
      }
      revertToOldValue(myStruc);
      decompose(myStruc);
			m++;

			free(myValues);
		}

	}
	else {
		/* we're at the bottom of the tree */	
		#ifdef DEBUG
			assert(!hasChildren(myStruc));
		#endif

    #ifdef FREEZEBASES
    int freeze = 0;
    #endif
		if (nValues != NULL) {
			initSeq(myStruc);
      #ifdef FREEZEBASES
      freeze = 1;
      #endif
		}
		
    #ifdef FREEZEBASES
		optimizeLeaf(myStruc, freeze);
    
    #else
		optimizeLeaf(myStruc);
    #endif
    
		DBL_TYPE bestValue;	
		
		int *bestSeq = (int *)calloc(myStruc->length, sizeof(int));
		DBL_TYPE *bestPr = (DBL_TYPE *)calloc((myStruc->length + 1)*(myStruc->length+1), sizeof(DBL_TYPE));
		DBL_TYPE *bestProbs = (DBL_TYPE *)calloc(myStruc->length, sizeof(DBL_TYPE));

		extractBestValues(myStruc, &bestValue, bestSeq, bestPr, bestProbs) ;
		
		int m = 0;
		while (m < mLeafopt && !satisfiedLeafObjective(bestValue,myStruc->targetObjective,designMode)) {
			initSeq(myStruc);
      #ifdef FREEZEBASES
      optimizeLeaf(myStruc, freeze);
      #else
      optimizeLeaf(myStruc);
      #endif
			if (favorableMove(myStruc->bestOptVal,bestValue,designMode)) {
				extractBestValues(myStruc, &bestValue, bestSeq, bestPr, bestProbs); 
			}
			m++;
		}

		setAndRevertToOldValue(myStruc, bestValue, bestSeq, bestPr, bestProbs); 
    assertNativeChecksum(myStruc);
		#ifdef DEBUG
			assertInSync(myStruc);
		#endif
		
		free(bestSeq);
		free(bestPr);
		free(bestProbs);


	}
	#ifdef DEBUG
		assertInSync(myStruc);
	#endif

}
