#ifndef MIN_H
#define MIN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*These functions are used in minimum free energy calculations,
  and closely mimic their partition function counterparts */

//Returns hairpin energy, unless nicked (returns NAD_INFINITY)
DBL_TYPE MinHairpin(int i, int j, int seq[], int seqlength, int **etaN);

//finds the minimum energy multiloop closed by i,j. (complexity = 3)
DBL_TYPE MinMultiloops(int i, int j, int seq[],
                       DBL_TYPE *Fms, DBL_TYPE *Fm, int seqlength,
                       int **etaN);

//finds minimum energy exterior loop
DBL_TYPE MinExteriorLoop(int i,int j, int seq[], int seqlength,
                         DBL_TYPE *F, int *nicks, int **etaN);

//finds the minimum interior or multiloop (complexity > 3)
DBL_TYPE MinInterior_Multi(int i, int j, int seq[], int seqlength,
                           DBL_TYPE *Fm, DBL_TYPE *Fb,
                           int *nicks, int **etaN);

//These functions find minimum energy interior loop (complexity = 3)
void MinFastILoops( int i, int j, int L, int seqlength, int seq[],
                    int **etaN, DBL_TYPE *Fb, DBL_TYPE *Fx, DBL_TYPE *Fx_2,
                    DBL_TYPE *minILoopEnergyBySize);

void makeNewFx( int i, int j, int seq[], int seqlength,
                int **etaN, DBL_TYPE Fb[], DBL_TYPE Fx[]);

void extendOldFx( int i, int j, int seqlength, DBL_TYPE Fx[], DBL_TYPE Fx_2[]);
DBL_TYPE MinInextensibleIL( int i, int j, int seq[], int seqlength,
                            DBL_TYPE Fb[], int **etaN, DBL_TYPE *minILoopEnergyBySize);

//Finds the minimum values for Fs, Fms, F, Fm respectively (complexity = 3)
void MakeFs_Fms( int i, int j, int seq[], int seqlength,
                 DBL_TYPE *Fs, DBL_TYPE *Fms, DBL_TYPE *Fb,
                 int *nicks, int **etaN);

void MakeF_Fm_N3( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *F, DBL_TYPE *Fs,
                  DBL_TYPE *Fms, DBL_TYPE *Fm,
                  int *nicks, int **etaN);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
