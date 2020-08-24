#ifndef NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__
#define NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void WriteOutput(double *x, double *G, int *CompIDArray, int LargestCompID, 
     int numSS, int numTotal, int nTotal, double kT, char *cxFile, 
     int SortOutput, char  *eqFile, double MolesWaterPerLiter, int quiet, 
     int NoPermID,int NUPACK_VALIDATE);

int Compare11(const void *p1, const void *p2);  // Comparison function for sorting

int Compare12(const void *p1, const void *p2);  // Comparison function for sorting

int Compare13(const void *p1, const void *p2);  // Comparison function for sorting

int Compare14(const void *p1, const void *p2);  // Comparison function for sorting

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_OUTPUTWRITER_H__ */
