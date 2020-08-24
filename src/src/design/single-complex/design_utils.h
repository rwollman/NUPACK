
#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>

#include <shared.h>
#include "design_constants.h"
#include "read_command_line.h"


#define NUM_BASES 4

#define baseR 6 //AG
#define baseY 11 //CU
#define baseM 7 //AC
#define baseK 10 //GU
#define baseS 8 //CG
#define baseW 9 //AU

#define baseV 12 //ACG
#define baseH 13 //ACU
#define baseD 14 //AGU
#define baseB 15 //CGU

#define baseN 0

int newNucsComps[20];
int nucArray[20][20];
int convertBaseToInt[100];
char convertIntToBase[20];


int maxBadStrings;
int maxBadStringLength;

#define BADSTRINGSNO 0
#define BADSTRINGSLENGTH 0


void freeLoadedStrings(void);
int isAFamily(int base);
int inNucArray(char family, char member);
DBL_TYPE uniformDist(void);
int get_rand_int(int max);
char chooseNotNuc(int nuc);
int getNuc(void);
int getNucNot(int nucIn);
char complementOf(char nuc1);
int isComplement(char nuc1, char nuc2);
void getComplement(int *str, int *newStr, int length);
void changeNuc(char *search);
int choose_nucleotide(void);
int * reverseStr(int *str, int length);
int verifyStructure(char *parens, char *filename);
int getPot(char *parens);
void calculateBaseDist(DBL_TYPE *baseDist, int *seq, int length);
int getNucNumber(char nuc);
int getIndyLength (int length, int *pairingPoint);
DBL_TYPE getCtargetDist(char *tempseq, DBL_TYPE targetGC, int *pairing, int paired);
void changeTtoU(char *seq);
void loadStringFromInitFile(FILE *initFile, int length, int *seq);
void loadBadStrings(char *filename);
int inBadStringsArray(char *seq);
int timeDifference(struct timeval *result, struct timeval *end, struct timeval *start);
int strFamCmp(char *seqA, char *badString);
void loadNucNumbers(void);
int convertParensToPairs(char *parens, int *pairs, int length, int *breaks);
void printPairs(int *pairs, int length);
int getIndependentLength(int *pairing, int length);
int fixSkeleton(int *skel, int *pairs, int length);
int isBaseSubset(int skel, int nuc);
int violatesPattern(int *seq, int *skel, int maxViolations);
void getWordAtPosition(char *seq, char *word, int wordLength, int pos);
int wordFinished(char *word);
int validSSMWordPosition(int length, int wordlength, int pos);
void initSeqToN(int *seq, int length);
void printSeq(int seqnum[]);
void convertBasesToInts(char *bases, int *ints);
void convertIntsToBases(int *ints, char* bases, int length);
int seqCmp(int *seqA, int *seqB, int length);
int removeStrandBreaks(char *parens, int *breaks, int length);
void createExpandedSeq(int *expanded, int *seq, int length, int numStrands, int *breaks);
void printSeqWithBreaks(char *seq, int numStrands, int *breaks);
int intArrayLen(int *arr);
int getSeqLength(int *seq);
int violatesSkeleton(int length, int *skel, int *seq);
void assertExpandedMatches(int *expanded, int *seq);
int floatsEqual(DBL_TYPE val1, DBL_TYPE val2);
#endif
