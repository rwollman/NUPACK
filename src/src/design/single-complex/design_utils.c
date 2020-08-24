#include "design_utils.h"


int** loadedBadStrings;

void loadNucNumbers(void) {
	
	newNucsComps[BASE_N] = BASE_N;
	newNucsComps[BASE_A] = BASE_U;
	newNucsComps[BASE_C] = BASE_G;
	newNucsComps[BASE_G] = BASE_C;
	newNucsComps[BASE_U] = BASE_A;
	newNucsComps[BASE_R] = BASE_Y;
	newNucsComps[BASE_Y] = BASE_R;
	newNucsComps[BASE_M] = BASE_K;
	newNucsComps[BASE_K] = BASE_M;
	newNucsComps[BASE_S] = BASE_S;
	newNucsComps[BASE_W] = BASE_W;
	newNucsComps[BASE_V] = BASE_B;
	newNucsComps[BASE_B] = BASE_V;	
	newNucsComps[BASE_H] = BASE_D;
	newNucsComps[BASE_D] = BASE_H;

	nucArray[BASE_A][BASE_A] = 1;
	nucArray[BASE_C][BASE_C] = 1;
	nucArray[BASE_G][BASE_G] = 1;
	nucArray[BASE_U][BASE_U] = 1;
	nucArray[BASE_R][BASE_A] = 1;
	nucArray[BASE_R][BASE_G] = 1;
	nucArray[BASE_Y][BASE_C] = 1;
	nucArray[BASE_Y][BASE_U] = 1;
	nucArray[BASE_M][BASE_A] = 1;
	nucArray[BASE_M][BASE_C] = 1;
	nucArray[BASE_K][BASE_G] = 1;
	nucArray[BASE_K][BASE_U] = 1;
	nucArray[BASE_S][BASE_G] = 1;
	nucArray[BASE_S][BASE_C] = 1;
	nucArray[BASE_W][BASE_A] = 1;
	nucArray[BASE_W][BASE_U] = 1;
	nucArray[BASE_V][BASE_A] = 1;
	nucArray[BASE_V][BASE_C] = 1;
	nucArray[BASE_V][BASE_G] = 1;
	nucArray[BASE_H][BASE_A] = 1;
	nucArray[BASE_H][BASE_U] = 1;
	nucArray[BASE_H][BASE_C] = 1;
	nucArray[BASE_B][BASE_G] = 1;
	nucArray[BASE_B][BASE_U] = 1;
	nucArray[BASE_B][BASE_C] = 1;
	nucArray[BASE_D][BASE_G] = 1;
	nucArray[BASE_D][BASE_A] = 1;
	nucArray[BASE_D][BASE_U] = 1;
	nucArray[BASE_N][BASE_A] = 1;
	nucArray[BASE_N][BASE_C] = 1;
	nucArray[BASE_N][BASE_G] = 1;
	nucArray[BASE_N][BASE_U] = 1;
	
	convertBaseToInt[(int)'N'] = BASE_N;	
	convertBaseToInt[(int)'A'] = BASE_A;
	convertBaseToInt[(int)'C'] = BASE_C;
	convertBaseToInt[(int)'G'] = BASE_G;
	convertBaseToInt[(int)'U'] = BASE_U;
	convertBaseToInt[(int)'R'] = BASE_R;
	convertBaseToInt[(int)'Y'] = BASE_Y;
	convertBaseToInt[(int)'M'] = BASE_M;
	convertBaseToInt[(int)'K'] = BASE_K;
	convertBaseToInt[(int)'S'] = BASE_S;
	convertBaseToInt[(int)'W'] = BASE_W;
	convertBaseToInt[(int)'V'] = BASE_V;
	convertBaseToInt[(int)'H'] = BASE_H;
	convertBaseToInt[(int)'B'] = BASE_B;
	convertBaseToInt[(int)'D'] = BASE_D;
	convertBaseToInt[(int)'+'] = STRAND_PLUS;
	
	convertIntToBase[BASE_N] = 'N';
	convertIntToBase[BASE_A] = 'A';
	convertIntToBase[BASE_C] = 'C';
	convertIntToBase[BASE_G] = 'G';
        if(DNARNACOUNT == DNA) {
          convertIntToBase[BASE_U] = 'T';
        } else {
	  convertIntToBase[BASE_U] = 'U';
        }

	convertIntToBase[BASE_R] = 'R';
	convertIntToBase[BASE_Y] = 'Y';
	convertIntToBase[BASE_M] = 'M';
	convertIntToBase[BASE_K] = 'K';
	convertIntToBase[BASE_S] = 'S';
	convertIntToBase[BASE_W] = 'W';
	convertIntToBase[BASE_V] = 'V';
	convertIntToBase[BASE_H] = 'H';
	convertIntToBase[BASE_B] = 'B';
	convertIntToBase[BASE_D] = 'D';
	convertIntToBase[STRAND_PLUS] = '+';
}


void freeLoadedStrings(void) {
	int i;
	for (i = 0; i < maxBadStrings; i++) {
		free(loadedBadStrings[i]);	
	}
}


void convertIntsToBases(int *ints, char* bases, int length) {
	int i;
	for (i = 0; i < length; i++) {
		bases[i] = convertIntToBase[ints[i]];	
	}
	assert(ints[length] == -1);
}


void convertBasesToInts(char *bases, int *ints) {
	int length = strlen(bases);
	int i;
	for (i = 0; i < length; i++) {
		ints[i] = convertBaseToInt[(int)bases[i]];	
	}
	ints[length] = -1;
	
}


void loadBadStrings(char *filename) {
	
	loadedBadStrings = (int **)calloc(1000, sizeof(int *));
	maxBadStringLength = 0;
	
	FILE *in = fopen(filename, "r");
	int i = 0;
	while (!feof(in)) {
		char *tmp = (char *)calloc(100, sizeof(char));
		fscanf(in, "%s\n", tmp);
		//free(tmp);
		if (strcmp(tmp, "") != 0) {
			int len = strlen(tmp);
			int *currArray = (int *)calloc(len + 1, sizeof(int));
			int j;
			for (j = 0; j < len; j++) {
				currArray[j] = convertBaseToInt[(int)tmp[j]];
			}
			currArray[j] = -1;
			loadedBadStrings[i] = currArray;
			if (len > maxBadStringLength) {
				maxBadStringLength = len;	
			}
			i++;
		}
		free(tmp);
	}
	
	maxBadStrings = i;
}


int seqCmp(int *seqA, int *seqB, int length) {

	int i;
	for (i = 0; i < length; i++) {
		if (seqA[i] != seqB[i]) {
			return 1;	
		}
	}
	return 0;
}


int intArrayLen(int *arr) {
	int len = 0;
	while(arr[len] != -1) {
		len++;	
	}
	return len;
}


int violatesPattern(int *seq, int *skel, int maxViolations) {
	
	if (maxBadStrings < 1) {
		return 0;
	}

	int badCount = 0;
	int length = intArrayLen(seq);
	#ifdef DEBUG
		assert(length == intArrayLen(skel));
	#endif
	int i, j;
	for (i = 0; i < length; i++) {
		for (j = 0; j < maxBadStrings; j++) {
			int *badWord = loadedBadStrings[j];
			int badWordLength = intArrayLen(badWord);
			int z = i;
			
			int counter = 0;
			int conflicts = 0;
			int conflictsInSkel = 0;
			while (z < length && counter < badWordLength) {
				if (nucArray[badWord[counter]][seq[z]]) {
					conflicts++;
					if (skel[z] != BASE_N) {
						conflictsInSkel++;	
					}
				}
				z++;
				counter++;
			}
			#ifdef DEBUG
				assert(conflicts <= badWordLength);
				assert(conflictsInSkel <= badWordLength);
			#endif
			
			if (conflicts == badWordLength && conflictsInSkel < badWordLength) {
				badCount++;	
			}		
			if (badCount > maxViolations) {

				#ifdef DEBUG
					assert(badCount == maxViolations + 1);
				#endif
				return badCount;	
			}
		}
	}

	return badCount;	
}


void initSeqToN(int *seq, int length) {
	int i;
	for (i = 0; i < length; i++) {
		seq[i] = BASE_N;
	}
}


void loadStringFromInitFile(FILE *initFile, int length, int *seq) {
	
	char *tmpString = (char *)calloc(2 * (length+1), sizeof(char));
	fscanf(initFile, "%s", tmpString);
	
	int tmpLength;
	tmpLength = strlen( tmpString);

	if (tmpLength != length) {
		fprintf(stderr, "ERROR IN INIT FILE:\n%s\n",tmpString);	
	}

	convertBasesToInts(tmpString, seq);
	free(tmpString);
	return;
	
}


int isAFamily(int base) {
	if (base == 0 || base > 4) {
		return 1;	
	}
	else {
		return 0;
	}
}


int isBaseSubset(int skel, int nuc) {
	return nucArray[skel][nuc];
}


void changeTtoU(char *seq) {
	
	int length = strlen(seq);
	int i;
	for (i = 0; i < length; i++) {
		if (seq[i] == 'T') {
			seq[i] = 'U';
		}
		
	}
	
}


void printPairs(int *pairs, int length) {
	int i;
	for (i = 0; i < length; i++) {
		printf("%d ",pairs[i]);	
	}
	printf("\n");
}


int getIndependentLength(int *pairing, int length) {
	int k;
	int indyLength = 0;
	for (k = 0; k < length; k++) {
		if (pairing[k] > k || pairing[k] < 0) {
				indyLength++;
		}
		
	}	
	return indyLength;
}


int removeStrandBreaks(char *parens, int *breaks, int length) {
	char *tempParens = (char *)calloc(length + 1, sizeof(char));
	int breakCounter = 0;
	int i;
	int counter = 0;
	for (i = 0; i < length; i++) {
		
		if (parens[i] != '+') {
			parens[counter] = parens[i];
			counter++;
			
		}
		else {
			breaks[breakCounter] = i - breakCounter;
			breakCounter++;		
		}
	}

	parens[counter] = 0;
	free(tempParens);

	return breakCounter;	
}


#ifdef DEBUG
void assertExpandedMatches(int *expanded, int *seq) {

	int i = 0;
	int counter = 0;
	while (expanded[i] != -1) {
		if (expanded[i] != STRAND_PLUS) {	
			if (expanded[i] != seq[counter]) {
				printf("%d,%d: %d,%d\n", i, counter, expanded[i],seq[counter]);	
			}
			assert(expanded[i] == seq[counter]);
			counter++;
		}
		i++;
	}
	
	
}
#endif

int getSeqLength(int *seq) {
	int length = 0;
	int i = 0;
	while (seq[i] != -1) {
		#ifdef DEBUG
			assert(seq[i] < 100);
		#endif
		length++;
		i++;
		
	}
	return length;
}

void createExpandedSeq(int *expanded, int *seq, int length, int numStrands, int *breaks) {
	#ifdef DEBUG
		assert(numStrands > 0);
		assert(seq[length] == -1);
	#endif
	
	if (numStrands == 1) {
		memcpy(expanded, seq, (length + 1) * sizeof(int));
		#ifdef DEBUG
			assert(expanded[length] == -1);
		#endif
		return;
	}
	
	int numBreaks = numStrands - 1;
	int i;

	expanded[length + numBreaks] = -1;

	
	#ifdef DEBUG
		assert(seq[length] == -1);
	#endif

	int newExpStart = 0;
	int newStart = 0;
	for (i = 0; i < numBreaks; i++) {
		memcpy(&expanded[newExpStart],&seq[newStart],(breaks[i]-newStart)*sizeof(int));
		expanded[breaks[i] + i] = STRAND_PLUS;

		newExpStart = breaks[i] + i + 1;
		newStart = breaks[i];
	}
	 
	for (i = newExpStart; i < length + numBreaks; i++) {
		expanded[i] = seq[i - numBreaks];	
		
	}
	
	#ifdef DEBUG
		assert(expanded[length + numBreaks] == -1);
		assert(getSeqLength(expanded) == length + numBreaks);
		assertExpandedMatches(expanded, seq);
	#endif 
}


int convertParensToPairs(char *parens, int *pairs, int length, int *breaks) {
	
	int numBreaks = removeStrandBreaks(parens, breaks, length);
	int numStrands = numBreaks + 1; 
	length = length - numBreaks;
	
	int i;
	for (i = 0; i < length; i++) {
		if (parens[i] == '(') {
				int potential = 1;
				int j = i+1;
				while (potential && j < length) {
					if (parens[j] == '(') {
							potential++;
					}
					else if (parens[j] == ')') {
					potential--;
					}
					j++;
				}
				assert(potential == 0);
				pairs[i] = j - 1;
				pairs[j-1] = i;
		}
		else if (parens[i] == '.') {
			pairs[i] = -1;	
		}
		
	}
	#ifdef DEBUG
		assert(numStrands > 0);
	#endif
	return numStrands;
}


void printSeq(int seqnum[]){
  int i = 0;
  while (seqnum[i] != -1) {
	  #ifdef DEBUG
	  	assert(seqnum[i] < 100 && seqnum[i] >= 0);
	  #endif
    printf("%c",convertIntToBase[seqnum[i]]);
    i++;
  }
  printf("\n");
}


void printSeqWithBreaks(char *seq, int numStrands, int *breaks) {

	int numBreaks = numStrands - 1;
	int i;
	int j = 0;
	for (i = 0; i < numBreaks && seq[j] != 0; i++) {
		while (j < breaks[i] && seq[j] != 0) {
			printf("%c",seq[j]);
			j++;
		}
		printf("+");	
		
	}
	while (seq[j] != 0) {
		printf("%c",seq[j]);
		j++;
	}
	
	printf("\n");
}


int fixSkeleton(int *skel, int *pairs, int length) {

	int i;
	int counter = 0;
	for (i = 0; i < length; i++) {
		skel[counter] = skel[i];

		if (skel[i] != STRAND_PLUS) {
			counter++;
		}
		
		
	}
	skel[counter] = -1;
	for (i = 0; i < counter; i++) {
		if (pairs[i] >= 0 && skel[i] != newNucsComps[skel[pairs[i]]] && skel[i] != BASE_N && skel[pairs[i]] != BASE_N) {
			fprintf(stderr, "ERROR IN CONSTRAINTS:\n %d at POS %d NOT COMPATIBLE WITH PAIR %d at POS %d\n", skel[i], i, skel[pairs[i]], pairs[i]);
			exit(1);

		}
		else {
			if (pairs[i] >= 0 && skel[i] != newNucsComps[skel[pairs[i]]]) {
				if (skel[i] == BASE_N && skel[pairs[i]] != BASE_N ) {
					skel[i] = newNucsComps[skel[pairs[i]]]; 	
				}
				else if (skel[i] != BASE_N && skel[pairs[i]] == BASE_N) {
					skel[pairs[i]] = newNucsComps[skel[i]];	
				}
				else if (skel[i] != newNucsComps[skel[pairs[i]]]) {
					fprintf(stderr, "2 ERROR IN CONSTRAINTS:\n %d at POS %d NOT COMPATIBLE WITH PAIR %d at POS %d\n", skel[i], i, skel[pairs[i]], pairs[i]);
					exit(1);	
				}
			}
		}

	}
	return 1;
	
}

void getWordAtPosition(char *seq, char *word, int wordLength, int pos) {
	assert(validSSMWordPosition(strlen(seq), wordLength, pos));
	strncpy(word, &seq[pos - wordLength + 1], wordLength);	
}


int validSSMWordPosition(int length, int wordLength, int pos) {
	if (pos >= (wordLength - 1) && pos < length) return 1;
	return 0;
}


int wordFinished(char *word) {
	int i;
	int length = strlen(word);
	for (i = 0; i < length; i++) {
		if (word[i] == 'N') return 0;	
	}
	return 1;
}


void calculateBaseDist(DBL_TYPE *base_dist, int *seq, int length) {

	int i;
	int *counter = (int*)calloc(NUM_BASES+1, sizeof(int));
	
	for (i = 0; i < length; i++) {
		counter[seq[i]]++;
	}
	
	for (i = 1; i < NUM_BASES+1; i++) {
		base_dist[i] = (DBL_TYPE)( (DBL_TYPE)counter[i] / (DBL_TYPE)length);
	}
	
	free(counter);
}


int get_rand_int(int max) {
	int rand_num = (int)(genrand_real2() * (DBL_TYPE)max);
	return rand_num;
}


int getPot(char *parens) {
	int length = strlen(parens);
	int i;
	int left = 0;
	int right = 0;
	for (i = 0; i < length; i++) {
		if (parens[i] == '(') left++;
			else if (parens[i] == ')') right++;
			
	}
	if (left != right) {
	}
	
	return 0;
}


int getNuc(void) {
	return get_rand_int(NUM_BASES) + 1;
}


int getNucNot(int nucIn) {
	int nucOut = nucIn;
	while (nucOut == nucIn) {
		nucOut = getNuc();
	}
	#ifdef DEBUG
		assert(nucOut >0);
		assert(nucOut != nucIn);
	#endif
	return nucOut;
}


void getComplement(int *str, int *newStr, int length) {
	int i;
	
	for (i = 0; i < length; i++) {
		newStr[i] = newNucsComps[str[i]];
	}

	newStr[length] = -1;	
}




int * reverseStr(int *str, int length) {
	int i;
	int *tempStr;
	
	tempStr = (int *)calloc(length+1, sizeof(int));
	memcpy(tempStr, str, length * sizeof(int));
	
	tempStr[length] = -1;
	for (i = 0; i < length; i++){
		tempStr[i] = str[length - 1 - i];
	}
	return tempStr;
	
}

int verifyStructure(char *parens, char *targetParens) {
	
	int retVal;

#ifdef DEBUG
	int length = strlen(parens);
	
	assert(length == strlen(targetParens));
#endif
	
	if (strcmp(parens, targetParens) != 0) {
		#ifdef DEBUG
			printf(":( ** MFE CONFLICTS EXIST\n");
				
			int k;
			printf("TARGET: ");
			for (k = 0; k < length; k++) {
				if (parens[k] == targetParens[k]) {
					printf("%c",targetParens[k]);
				}
				else {
					printf("!%c",targetParens[k]);
				}
				
				
			}
			printf("\n");
			printf("SENT: %s\n",parens);
		#endif
		
		retVal = 0;
	}
	else {
		#ifdef DEBUG
			printf("*** MFE VERIFIED\n");
		#endif
		retVal = 1;
	}
	
	return retVal;
	
}


int violatesSkeleton(int length, int *seq, int *skel) {
	int i;
	for (i = 0; i < length; i++) {
		if (!nucArray[skel[i]][seq[i]]) {
			printf("VIOLATION: %d: %d:%d\n", i, skel[i], seq[i]);
			printSeq(skel);
			printSeq(seq);
			return 1;	
		}
		
	}
	return 0;
}


int getIndyLength (int length, int *pairing) {
	
	int i;
	int counter = 0;
	for (i = 0; i < length; i++) {
		if (pairing[i] < i) {
			counter++;
		}
		
		
	}
	return counter;
}



int timeDifference (struct timeval *difference, struct timeval *end, struct timeval *start)
{
	if (end->tv_usec < start->tv_usec) {
		int nsec = (start->tv_usec - end->tv_usec) / 1000000 + 1;
		start->tv_usec -= 1000000 * nsec;
		start->tv_sec += nsec;
	}
	if (end->tv_usec - start->tv_usec > 1000000) {
		int nsec = (end->tv_usec - start->tv_usec) / 1000000;
		start->tv_usec += 1000000 * nsec;
		start->tv_sec -= nsec;
	}

	difference->tv_sec = end->tv_sec - start->tv_sec;
	difference->tv_usec = end->tv_usec - start->tv_usec;
	
	return end->tv_sec < start->tv_sec;
}

int floatsEqual(DBL_TYPE val1, DBL_TYPE val2) {
	if (fabsl(val1 - val2) < FLT_EPSILON) {
		return 1;	
	}
	return 0;
	
}
