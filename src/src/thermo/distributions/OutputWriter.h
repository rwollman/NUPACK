#ifndef NUPACK_THERMO_DISTRIBUTIONS_OUTPUTWRITER_H__
#define NUPACK_THERMO_DISTRIBUTIONS_OUTPUTWRITER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void WriteOutput(double *mEq, double **Pmn, double *G, int **A, int *CompIDArray, 
     int *PermIDArray, int LargestCompID, int numSS, int numTotal, 
     int nTotal, int nComments, int maxm0, double kT, char *cxFile, 
     int SortOutput, char *distFile, int quiet, int NoPermID,int NUPACK_VALIDATE);

int Compare21(const void *p1, const void *p2);

int Compare22(const void *p1, const void *p2);

int Compare23(const void *p1, const void *p2);

int Compare24(const void *p1, const void *p2);

int Compare25(const void *p1, const void *p2);

int Compare26(const void *p1, const void *p2);

int Compare27(const void *p1, const void *p2);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_DISTRIBUTIONS_OUTPUTWRITER_H__ */
