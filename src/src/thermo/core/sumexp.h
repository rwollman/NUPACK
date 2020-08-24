#ifndef __SUMEXP_H__
#define __SUMEXP_H__

#include "ene.h"
#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

//Hairpin energy (exp)
DBL_TYPE ExplHairpin( int i, int j, int seq[], int seqlength, int **etaN);

//Calculates the contribution to the partition function of multiloops (non-nicked)
DBL_TYPE SumExpMultiloops( int i, int j, int seq[],
                           DBL_TYPE *Qms, DBL_TYPE *Qm, int seqlength,
                           int **etaN);
//Calculates the contribution of exterior loops
DBL_TYPE SumExpExteriorLoop( int i,int j, int seq[], int seqlength,
                             DBL_TYPE *Q,
                             int *nicks, int **etaN);

//Computes Qs, Qms  (pairs in exterior loops, multi loops)
void MakeQs_Qms( int i, int j, int seq[], int seqlength,
                 DBL_TYPE *Qs, DBL_TYPE *Qms, DBL_TYPE *Qb,
                 int *nicks, int **etaN);

//Computes Q, Qm for complexity = 3 algorithm
void MakeQ_Qm_N3( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qs,
                  DBL_TYPE *Qms, DBL_TYPE *Qm,
                  int *nicks, int **etaN);

//void MakeQ_Qm_N4( int i, int j, int seq[], int seqlength,
//                  DBL_TYPE *Q, DBL_TYPE *Qm, DBL_TYPE *Qb );

//Calculates contribution of interior loops and multiloops
//for complexity >= 4 methods.  Less memory usage, more time required
//than complexity = 3 method.  Same result.
DBL_TYPE SumExpInterior_Multi( int i, int j, int seq[], int seqlength,
                               DBL_TYPE *Qm, DBL_TYPE *Qb);

//Efficiently calculates the contribution of large interior loops
void fastILoops( int i, int j, int L, int seqlength, int seq[],
                 int **etaN,
                 DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2,
                 DBL_TYPE *Qb_bonus);


//makeNewQx creates new "extensible" base cases for the interval i,j.
void makeNewQx( int i, int j, int seq[], int seqlength,
                int **etaN, DBL_TYPE Qb[], DBL_TYPE Qx[]);
//extendOldQx extends Qx for the i-1, j+1 case
void extendOldQx( int i, int j, int seqlength,
                  DBL_TYPE Qx[], DBL_TYPE Qx_2[]);

//Directly calculates the contribution of small interior loops
DBL_TYPE SumExpInextensibleIL( int i, int j, int seq[], int seqlength,
                               DBL_TYPE Qb[],  int **etaN);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
