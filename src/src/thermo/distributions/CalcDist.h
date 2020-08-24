#ifndef NUPACK_THERMO_DISTRIBUTIONS_CALCDIST_H__
#define NUPACK_THERMO_DISTRIBUTIONS_CALCDIST_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void CalcDist(double *mEq, double **Pmn, int **A, double *G, int *m0, double M, 
        int numSS, int numTotal, double MaxSizeLambda, char *LambdaFile, 
        double kT, int WriteLambda, int *CompIDArray, int *PermIDArray,
        int quiet, char *logFile, int WriteLogFile, int NoPermID);

int next(int *mComplex, int *m0, int **A, int numSS, int numComplex, int numTotal,
   int *mMax, int LastInc);

void UpdateLambda(int ***Lambda, int k, int *mComplex, int *m0, int **A, int numSS,
      int numTotal);

int NegCheck(int *mComplex, int *m0,int **A,int numSS,int numTotal);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_DISTRIBUTIONS_CALCDIST_H__ */
