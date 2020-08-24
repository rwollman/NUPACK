#ifndef DESIGN_H
#define DESIGN_H


#include <shared.h>
#include <thermo.h>
#include "design_utils.h"

#ifdef DEBUG
#include "design_test.h"
#endif

//int naType;
//int dangles;
int initMethod;
int leafsum_flag;
int quickFlag;
int mReopt;
int mLeafopt;
DBL_TYPE nRatio;
DBL_TYPE ppairsCutoff;

int fakeSize;

DBL_TYPE mpC;
DBL_TYPE mpCprime;

int bypassDesign;
int bypassHierarchy;
int bypassGuidance;
int initMode;
int designMode; 
int loadInit;
int hMin;

FILE *trailOut;

#ifdef EARLY_TERMINATE
    struct timeval *global_timeval;
    char *timerFilename;
#endif

time_t tstart;
char *trailFile;

typedef struct intList {
	int value;
	struct intList *next;
	
} intList;

typedef struct mutation {
	int posA;
	int nucA;
	int posB;
	int nucB;	
	
} mutation;


typedef struct secStruc {
  
  int id;
  
  int length;
  int relevantLength;
  
  int k;
  int l;
  int indyLength;
  int expandedLength;
  
  //int parentLow;
  //int parentHigh;
  
  int checked;
  
  int level;
  
  DBL_TYPE targetObjective;
  
  int *fakeBases;
  int *purelyNativeBases;
  
  int numStrands;
  int *strandBreaks;
  
  int *tempPairing;
  char *tempFolding;
  char *expFolding;

  int loopStart;
		
  struct secStruc *leftChild;
  struct secStruc *rightChild;
  struct secStruc *parent;
  
  int amILeft;
  int amIChild;
  
  int *parentMapping;
  int *leftChildMapping;
  int *rightChildMapping;
  
  char *foldFile;
  DBL_TYPE pfVal;
  
  int *seq;
  int *skel;
  int *expandedSeq;
  
  DBL_TYPE bestNVal;
  DBL_TYPE nVal;
  int *bestSeq;
  DBL_TYPE *bestPr;
  DBL_TYPE *bestProbs;
  
  int *bestConflictArray;
  int *conflictArray;
  
  int *bestPairs;
  int *mfePairs;
  int bestMFEDifference;
  int mfeDifference;
  
  DBL_TYPE optVal; 
  DBL_TYPE bestOptVal;
  
  DBL_TYPE pVal;
  DBL_TYPE bestPVal;
  
  int bestNDifference;
  int nDifference;
  
  DBL_TYPE *myPairPr;
  DBL_TYPE *myProbs;
  
  int timesOptimized;
  int totalTimeSpent;
  int nsCalls;
  
  DBL_TYPE pfTimeSpent;
  DBL_TYPE mfeTimeSpent;
  DBL_TYPE pTimeSpent;

  
  int maximumUnfavorableEvaluations;
  int maximumMutationSearch;
  
  int *initNucs;
  int *initComps;
  int *initMap;
  int numInitBases;
  
  int ssmWordLength;
  
  #ifdef LEAFCORRECTION
  DBL_TYPE *leafArray;
  int *leafCheck;
  DBL_TYPE *leafNative;
  int *leafNativeCheck;
  int modified;
  #endif
  
  int failures;
  
} secStruc;


typedef struct stepList {
	int pos;
	int from;
	int to;
	struct stepList *next;
} stepList;


typedef struct conflictList {
	int pos;
	int base;
	struct conflictList *next;
} conflictList;


void freeEngine(void);
int* nucToIntTable(void);
void getNValues(secStruc *myStruc, DBL_TYPE *nValues);
void initConflictList(conflictList *theList);
void freeConflictList(conflictList *theList);
void freeStruc(secStruc *myStruc);
void designSeq(secStruc *myStruc, DBL_TYPE *nValues);
int conflictListSize(conflictList *theList);
int inEitherList(stepList *theList, conflictList *conflicts, int pos, int base, int *calcsaves);
void addToConflictList(conflictList *theList, int pos, int base);
int changeUnpaired(int pos, stepList *theList, conflictList *conflicts, secStruc *myStruc, int *calcsaves, int *currentViolations, int acceptablePatternViolations);
void createTree(secStruc *parentStruc);
void parseStructure(secStruc *myStruc, char *foldFile);
void initSeq(secStruc *myStruc);
void computeP(secStruc *myStruc);
void computeNSstar(secStruc *myStruc);
void computeMFE(secStruc *myStruc);
void initStructures(secStruc *wholeStruc, int justEval, int alreadyInit);
int getSSMWordLength(int length);
int violatesSSM(secStruc *myStruc, int wordLength);
long int q2d(int *entry, int wordLength);
int mapFakeBases(secStruc *parent, secStruc *child);
int performMutation(secStruc *myStruc, int pos, stepList *theList, conflictList *conflicts, int *calcsaves, int *currentViolations, int acceptablePatternViolations);
void initStepList(stepList *theList);
int getMFEDifference(secStruc *myStruc);
int getNDifference(secStruc *myStruc);
void setAllBestValues(secStruc *myStruc);
int findSplitPoint(int *pairing, int iLength, int numStrands, int *breaks);
void initializeMapping(int *mapping, int length);
void setNewBestP(secStruc *myStruc);
void setNewBestN(secStruc *myStruc);
void setNewBestMFE(secStruc *myStruc);
int ssmWordLength(int length);
int violatesSkeleton(int length, int *skel, int *seq);
void setupInitNucs(secStruc *myStruc);
void resetOptVal(secStruc *myStruc);
void initEngine(char *psFile);
void mapSequenceDown(secStruc *parent, secStruc *child, int struc);
void mapSequenceUp(int *parentSeq, int *childSequence, secStruc *child, int left);
void resetNodeVals(secStruc *myStruc, DBL_TYPE optVal);
void setNewBestValue(secStruc *myStruc);
void initStruc(secStruc *myStruc, int length, int expandedLength);
void clearHash(long int *hash, int wordLength);
long int inHash(long int *hash,int *entry, int wordLength) ;
void removeFromHash(long int *hash, int *entry, int wordLength);
int invalidComplementInHash(long int *hash, int *word, int pos, int *pairing, int wordLength);
void addToHash(long int *hash, int *entry, int pos, int wordLength);
void resetTriedOptions(int *tried, int numberOfOptions, int *seq, int *skel, int pos);
int optionsAtThisPosition(int *tried, int numberOfPossibleOptions);
void pureRandomInit(secStruc *myStruc);
int seqAlreadyInit(secStruc *myStruc);
void initHardBases(secStruc *myStruc);
void freeStepList(stepList *theList);
void clearStepList(stepList *theList) ;
void addToStepList(stepList *theList, int pos, int from, int to);
int inStepList(stepList *theList, int pos, int to);
void resetSeqPosToBest(secStruc *myStruc, int pos);
void mutateInheritedDefects(secStruc *myStruc, conflictList *conflicts);
void calculateProbs(secStruc *myStruc);
int inConflictList(conflictList *theList, int pos, int base);
void mapStrandBreaksDown(secStruc *parent, secStruc *child);
int hasChildren(secStruc *parent);
void decompose(secStruc *parent);
void merge(secStruc *left, secStruc *right);
void optimizeNode(secStruc *myStruc, DBL_TYPE *nValues);

#ifdef FREEZEBASES
void optimizeLeaf(secStruc *myStruc, int freeze);
#else
void optimizeLeaf(secStruc *myStruc);
#endif
void convertConflictArrayToList(secStruc *myStruc, conflictList *list);
int nearStrandBreak(int pos, int length, int *pairing, int numStrands, int *breaks, int hmin);
void fixPatternViolations(secStruc *myStruc);
void mutateRandomBaseOrPair(secStruc *myStruc, mutation *candidateMutation);
void mutateRandomConflict(secStruc *myStruc, mutation *candidateMutation);
int mutationProhibited(secStruc *myStruc, mutation *candidateMutation, int *currentPatternViolations, int acceptablePatternViolations);
void applyMutation(secStruc *myStruc, mutation *candidateMutation);
void calculateObjective(secStruc *myStruc);
int mutationUnfavorable(mutation *candidateMutation, stepList *list);
void mutationAtPos(secStruc *myStruc, mutation *candidateMutation, int pos);
void revertToOldValue(secStruc *myStruc);
DBL_TYPE calculateLeafSum(secStruc *myStruc);
void getLeafNValues(secStruc *myStruc, DBL_TYPE *nValues);
void designHeader(int argc, char ** argv, char * filename, struct tm * start_time,char * structure);
void outputConstantParameters(char *filename);
int helixExtension(int pos, int *seq, int *pairing);
void printTree(secStruc *myStruc, int level, FILE*);
DBL_TYPE getNativeN(secStruc *myStruc);
int satisfiedObjective(secStruc *myStruc);
void setAndRevertToOldValue(secStruc *myStruc, DBL_TYPE bestVal, int *bestSeq, DBL_TYPE *bestPr, DBL_TYPE *bestProbs);
void extractBestValues(secStruc *myStruc, DBL_TYPE *bestVal, int *bestSeq, DBL_TYPE *bestPr, DBL_TYPE *bestProbs); 
DBL_TYPE getChildContributingN(secStruc *myStruc);
#ifdef DEBUG
void assertWatsonCrick(secStruc *myStruc);
void assertValueInSync(secStruc *myStruc);
void assertInSync(secStruc *myStruc);
void assertValidStrandBreaks(secStruc *myStruc);


#endif
int satisfiedLeafObjective(DBL_TYPE bestValue, DBL_TYPE targetObjective, int optType);
int favorableMove(DBL_TYPE newVal, DBL_TYPE oldVal, int optType);

#ifdef LEAFCORRECTION
void assertChecksum(secStruc *myStruc);
void copyLeafArrayFromChild(secStruc *parent, secStruc *child);
void overwriteChildrenLeafArrays(secStruc *myStruc);
DBL_TYPE sumLeafIDs(secStruc *myStruc);
DBL_TYPE sumNativeIDs(secStruc *myStruc);
void assertNativeChecksum(secStruc *myStruc);
#endif

#endif
