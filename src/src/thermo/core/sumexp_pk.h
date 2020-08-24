#ifndef __SUMEXP_PK_H__
#define __SUMEXP_PK_H__

#include "ene.h"
#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

//Calculates Qb for complexity = 5.
DBL_TYPE SumExpQb_Pk( int i, int j, int seq[], int seqlength,
                      DBL_TYPE *Qp, DBL_TYPE *Qm );

//Calculates Qp for complexity = 5.  Has extra cases not in JCC papers
//that allow for single pairs to form pseudoknot "gap-spanning" regions.
DBL_TYPE SumExpQp_N5( int i, int j, int seq[], int seqlength,
                      DBL_TYPE *Qgl, DBL_TYPE *Qgr, DBL_TYPE *Qg, DBL_TYPE *Qz);

//The following 6 functions calculate the respective Q values in their names.
//Diagrams for these are in our JCC papers
void MakeQg_N5( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
                DBL_TYPE *Qm, DBL_TYPE *Qgls, DBL_TYPE *Qgrs, DBL_TYPE *QgIx,
                DBL_TYPE *QgIx_2, short *possiblePairs);
void MakeQgls( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm,  DBL_TYPE *Qgls);
void MakeQgrs( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm, DBL_TYPE *Qgrs);
void MakeQgl( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Qg, DBL_TYPE *Qgl, DBL_TYPE *Qz);
void MakeQgr( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Qgr, DBL_TYPE *Qgl, DBL_TYPE *Qz);
void MakeQ_Qm_Qz( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qm, DBL_TYPE *Qz,
                  DBL_TYPE *Qb, DBL_TYPE *Qp);

//A fast interior loops method for gap matrices.  Analogous to the fastIloops routine
void fastIloop_Qg(int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Qg, DBL_TYPE *QgIx, DBL_TYPE *QgIx_2,
                  short *possiblePairs);
DBL_TYPE SumexplInextensibleIL_Qg( int i, int j, int d, int e,
                                   int seq[], int seqlength,
                                   DBL_TYPE *Qg);
void makeNewQgIx( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Qg, DBL_TYPE *QgIx);
void extendOldQgIx( int i, int j, int d, int e, int seq[], int seqlength,
                    DBL_TYPE *Qg, DBL_TYPE *QgIx, DBL_TYPE *QgIx_2);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
