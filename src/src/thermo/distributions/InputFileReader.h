#ifndef NUPACK_THERMO_DISTRIBUTIONS_INPUTFILEREADER_H__
#define NUPACK_THERMO_DISTRIBUTIONS_INPUTFILEREADER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
       int **numPermsArray, int *nComments, char *cxFile, char *countFile, 
       int quiet);

void ReadInputFiles(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
        int **m0, double *M, int *numSS, int *numSS0, int *numTotal,
        int *numPermsArray, char *cxFile, char *countFile, double *kT,
        int Toverride, char *logFile, char *eqFile, int quiet, 
        int WriteLogFile);

void ReadInputFilesPerm(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
      int **m0, double *M, int *numSS, int *numSS0, 
      int *newnTotal, int nTotal, char *cxFile, char *countFile, 
      double *kT, int Toverride, char *logFile, char  *eqFile, 
      int quiet, int WriteLogFile);

int InputCompare(const void *p1, const void *p2);

int InputComparePerm(const void *p1, const void *p2);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_DISTRIBUTIONS_INPUTFILEREADER_H__ */
