#include <stdio.h>

int ReadCommandLine(int, char**);

int ReadInputFileComplexes(char *filePrefix, int *nStrands,
                           char ***seqs, int **seqlength,
                           int *maxLength, int *maxComplexSize);

void printHeader(int nStrands, char **seqs, int maxComplexSize,
                 int totalOrders, int nNewPerms, int nSets, int nNewComplexes,
                 FILE *F_cx, int nargs, char **argv, int isPairs);

void print_deprecation_info(FILE *out);
