#include "design_constants.h"

#ifdef DEBUG

#include "design_test.h"

void unitTester() {
	testParensConverter(); 
	testConstraints();
	testPatternPrevention();
	//testNSTheta();
}

void testNSTheta() {
	
	
	char *seq1 = "CCCCCAGACGGGGGG";
	int pairing1[] = {14,13,12,11,10,-1,-1,-1,-1,-1,4,3,2,1,0};
	int fakeBases1[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int numStrands1 = 1;
	int *breaks1 = NULL;

	int length1 = strlen(seq1);
	DBL_TYPE *pairPr1 = (DBL_TYPE *)calloc((length1 + 1) * (length1 + 1), sizeof(DBL_TYPE));
	
	
	char *seq2 = "ACGUACUGACGCGCGCGCGCCAGACUGCGCGCGCGC";
	int pairing2[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,34,33,32,31,30,29,28,27,26,25, -1,-1,-1,-1,-1, 19, 18, 17,16,15,14,13,12,11,10};
	int fakeBases2[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int numStrands2 = 2;
	int breaks2[] = {15};

	int length2 = strlen(seq2);
	DBL_TYPE *pairPr2 = (DBL_TYPE *)calloc((length2 + 1) * (length2 + 1), sizeof(DBL_TYPE));
	
	
	char *seq3 = "ACGAUCGAUUAGCUAGCUAACGAUCGAUCGAUUCAGACGAUCGAUUAGCUAGCUAACGAUCGAUCGAUUCAGACGAUCGAUUAGCUAGCUAACGAUCGAUCGAUUCAG";
	int pairing3[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	int fakeBases3[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int numStrands3 = 1;
	int *breaks3 = NULL;

	int length3 = strlen(seq3);
	DBL_TYPE *pairPr3 = (DBL_TYPE *)calloc((length3 + 1) * (length3 + 1), sizeof(DBL_TYPE));
	
	
	DBL_TYPE val1 = createNSThetaCase(seq1, pairing1, fakeBases1, numStrands1, breaks1, pairPr1);

	DBL_TYPE val2 = createNSThetaCase(seq2, pairing2, fakeBases2, numStrands2, breaks2, pairPr2);

	DBL_TYPE val3 = createNSThetaCase(seq3, pairing3, fakeBases3, numStrands3, breaks3, pairPr3);

	
	printf("VAL: %LE:%LE\n", val1, val2);

	int i;
	for (i = 0; i < 10; i++) {
		
		DBL_TYPE tempVal1a = createNSThetaCase(seq1, pairing1, fakeBases1, numStrands1, breaks1, pairPr1);
		DBL_TYPE tempVal1b = val1;
		
		
		assert(floatsEqual(tempVal1b, tempVal1a));
		
		DBL_TYPE tempVal2a = createNSThetaCase(seq2, pairing2, fakeBases2, numStrands2, breaks2, pairPr2);
		DBL_TYPE tempVal2b = val2;
		
		
		assert(floatsEqual(tempVal2b, tempVal2a));
		
		
		DBL_TYPE tempVal3a = createNSThetaCase(seq3, pairing3, fakeBases3, numStrands3, breaks3, pairPr3);
		DBL_TYPE tempVal3b = val3;
		
		
		assert(floatsEqual(tempVal3b, tempVal3a));
		
	}
	
	


	
}


DBL_TYPE createNSThetaCase(char *seq, int *pairing, int *fakeBases, int numStrands, int *breaks, DBL_TYPE *pm) {
	
	int length = strlen(seq);
	
	int *intSeq = (int *)calloc(length + 1, sizeof(int));
	int *intExpandedSeq = (int *)calloc(length + numStrands,sizeof(int));
		
	pairPr = pm;
	
	
	convertBasesToInts(seq, intSeq);
	

	createExpandedSeq(intExpandedSeq, intSeq, length, numStrands, breaks);
		
	DBL_TYPE pfVal;
	
	
	DBL_TYPE optVal = nsStarPairsOrParensCorrected( length, numStrands, intExpandedSeq, pairing, NULL, fakeBases, &pfVal, 3, DNARNACOUNT, DANGLETYPE, 
							TEMP_K - ZERO_C_IN_KELVIN, SODIUM_CONC, MAGNESIUM_CONC, USE_LONG_SALT);

	return optVal;
	
}

void testParensConverter() {

	printf("TESTING PARENS CASE 1\n");
	char *parensWithPlus = "(((((+)))))";
	char *parensWithoutPlus = "((((()))))";
	int breaks[] = {5};
	int numBreaks = 1;
	createParensConverterCase(parensWithPlus, parensWithoutPlus, breaks, numBreaks);	
	printf("PASSED\n");
	
	printf("TESTING PARENS CASE 2\n");
	parensWithPlus = "(((((+)))+))...";
	parensWithoutPlus = "((((()))))...";
	int breaks2[] = {5, 8};
	numBreaks = 2;
	createParensConverterCase(parensWithPlus, parensWithoutPlus, breaks2, numBreaks);	
	printf("PASSED\n");

	
	printf("TESTING PARENS CASE 3\n");
	parensWithPlus = "(((((...)))))";
	parensWithoutPlus = "(((((...)))))";
	int *breaks3 = NULL;
	numBreaks = 0;
	createParensConverterCase(parensWithPlus, parensWithoutPlus, breaks3, numBreaks);
	printf("PASSED\n");
	
	printf("PASSED ALL PARENS CONVERSION TESTS!\n");

	
}

void createParensConverterCase(char *parensWithPlus, char *parensWithoutPlus, int *breaks, int numBreaks) {
	char *pWithP = (char *)calloc(strlen(parensWithPlus) + 1, sizeof(char));
	char *pWithoutP = (char *)calloc(strlen(parensWithoutPlus) + 1, sizeof(char));
	int *checkBreaks = (int *)calloc(strlen(parensWithPlus), sizeof(int));
	strcpy(pWithP, parensWithPlus);
	strcpy(pWithoutP, parensWithoutPlus);
	int totalBreaks = removeStrandBreaks(pWithP, checkBreaks, strlen(pWithP));
	assert(totalBreaks == numBreaks);
	assert(!strcmp(pWithP,pWithoutP));
	int i;
	for (i = 0; i < totalBreaks; i++) {
		printf("%d:%d\n",checkBreaks[i], breaks[i]);
		assert(checkBreaks[i] == breaks[i]);		
	}
	
	free(pWithP);
	free(pWithoutP);
	free(checkBreaks);
	
}


void testConstraints() {
	
	char *seq = "CCCCCCCCCCC";
	char *skel ="NNNNNNNNNNN";
	createSkeletonCase(seq, skel, 0);
	
	seq = "AAAAA";
	skel ="CCCCC";
	createSkeletonCase(seq, skel, 1);

	seq = "AACAA";
	skel ="WWSWW";
	createSkeletonCase(seq, skel, 0);
	
		
	printf("PASSED ALL CONSTRAINT PREVENTION UNIT TESTS!\n");
	
	
}

void createSkeletonCase(char *seq, char *skel, int violates) {
	assert(strlen(seq) == strlen(skel));
	
	int *intSeq = (int *)calloc(strlen(seq)+1, sizeof(int));
	int *intSkel = (int *)calloc(strlen(seq)+1,sizeof(int));
		
	convertBasesToInts(seq, intSeq);
	convertBasesToInts(skel,intSkel);
	
	char *testSeq = (char *)calloc(strlen(seq)+1, sizeof(char));
	convertIntsToBases(intSeq, testSeq, strlen(seq));
	assert(strcmp(testSeq, seq) == 0);
	
	assert(violatesSkeleton(strlen(seq), intSeq,intSkel) == violates);

	free(intSeq);
	free(intSkel);
	free(testSeq);	
	
}

void testPatternPrevention() {
	
	if (maxBadStrings < 1) {
		return;	
	}
	
	printf("THIS UNIT TEST REQUIRES A SPECIFIC PREVENTED STRINGS FILE\n");
	
	char *seq = "CCCCCCCCCCC";
	char *skel ="NNNNNNNNNNN";
	createPPCase(seq, skel, 0, 1);
	createPPCase(seq, skel, 2, 3);

	
	seq = "CACACACACA";
	skel ="NNNNNNNNNN";
	createPPCase(seq, skel, 0, 0);
	
	seq = "CCCACACAACCC";
	skel ="CCCNNNNNNCCC";
	createPPCase(seq, skel, 0, 0);
	
	seq = "CCCACACACCCC";
	skel ="CCCNNNNNNCCC";
	createPPCase(seq, skel, 1, 1);

	seq = "GGGGAAAACCCC";
	skel ="NNNNNNNNNNNN";
	createPPCase(seq, skel, 1,2);
	createPPCase(seq, skel, 2,3);
	createPPCase(seq, skel, 3,3);
	
	seq = "CCCCCCCCC";
	skel ="SSSSSSSSS";
	createPPCase(seq, skel, 0, 0);
	
	seq = "CCCCCCCCC";
	skel ="SSSSNSSSS";
	createPPCase(seq, skel, 0, 1);
	createPPCase(seq, skel, 1, 2);
	createPPCase(seq, skel, 2, 3);
	createPPCase(seq, skel, 3, 4);
	createPPCase(seq, skel, 4, 5);
	createPPCase(seq, skel, 5, 6);
	createPPCase(seq, skel, 6, 7);
	createPPCase(seq, skel, 7, 8);
	createPPCase(seq, skel, 8, 8);

	
		
	printf("PASSED ALL PATTERN PREVENTION UNIT TESTS!\n");
	
}


void createPPCase(char *seq, char *skel, int violates, int actualViolations) {
	assert(strlen(seq) == strlen(skel));
	
	int *intSeq = (int *)calloc(strlen(seq)+1, sizeof(int));
	int *intSkel = (int *)calloc(strlen(seq)+1,sizeof(int));
	
	convertBasesToInts(seq, intSeq);
	convertBasesToInts(skel,intSkel);
	
	char *testSeq = (char *)calloc(strlen(seq)+1, sizeof(char));
	convertIntsToBases(intSeq, testSeq, strlen(seq));
	assert(strcmp(testSeq, seq) == 0);

	assert(violatesPattern(intSeq,intSkel,violates) == actualViolations);
	
	free(intSeq);
	free(intSkel);
	free(testSeq);
}

#endif
