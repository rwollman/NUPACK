#ifndef NUPACK_THERMO_DISTRIBUTIONS_READCOMMANDLINE_H__
#define NUPACK_THERMO_DISTRIBUTIONS_READCOMMANDLINE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>

void ReadCommandLine(int nargs, char **args, char *cxFile, char *countFile, 
         char *logFile, char *distFile, char *LambdaFile, 
         int *SortOutput, int *WriteLambda, double *MaxSizeLambda,
         double *kT, int *quiet, int *WriteLogFile, int *Toverride,
         int *NoPermID, int * NUPACK_VALIDATE, int *v3);

void DisplayDistributionsHelp(int DummyArgument);

void print_deprecation_info(FILE *out);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_DISTRIBUTIONS_READCOMMANDLINE_H__ */
