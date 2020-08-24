#ifndef NUPACK_THERMO_CONCENTRATIONS_READCOMMANDLINE_H__
#define NUPACK_THERMO_CONCENTRATIONS_READCOMMANDLINE_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>

void ReadCommandLine(int nargs, char **args, char *cxFile, char *conFile, 
         char *logFile, char *eqFile, char *pairsFile, char *fpairsFile,
         int *SortOutput, int *MaxIters, double *tol, double *kT,
         int *MaxNoStep, int *MaxTrial, double *PerturbScale, int *quiet,
         int *WriteLogFile, int *Toverride, int *NoPermID, 
         int *DoBPfracs, unsigned long *seed, double *cutoff, int * NUPACK_VALIDATE,
         int *v3);

void DisplayHelpConc(void);

void print_deprecation_info(FILE *out);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_READCOMMANDLINE_H__ */
