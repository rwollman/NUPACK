#ifndef NUPACK_THERMO_CONCENTRATIONS_FRACPAIR_H__
#define NUPACK_THERMO_CONCENTRATIONS_FRACPAIR_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void FracPair(int numSS, int nTotal, int quiet, int NoPermID, int LargestCompID, 
        int *numPermsArray, char *eqFile, char *conFile, char *pairsFile, 
        char *fpairsFile, double cutoff, int NUPACK_VALIDATE);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_FRACPAIR_H__ */
