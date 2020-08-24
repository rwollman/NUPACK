#ifndef NUPACK_THERMO_CONCENTRATIONS_INPUTFILEREADER_H__
#define NUPACK_THERMO_CONCENTRATIONS_INPUTFILEREADER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
       int **numPermsArray, char *cxFile, char *conFile, int quiet);

double ReadInputFiles(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
          double **x0, int *numSS, int *numSS0, int *numTotal, 
          int *numPermsArray, char *cxFile, char *conFile, double *kT, 
          int Toverride, char  *logFile, char  *eqFile, 
          char *fpairsFile, int quiet, int WriteLogFile, int DoBPfrac,
          int NoPermID);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_INPUTFILEREADER_H__ */
