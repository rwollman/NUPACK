#ifndef __PAIRSPR_H__
#define __PAIRSPR_H__

#include "pfuncUtils.h"


#ifdef __cplusplus
extern "C" {
#endif

//These files deal with calculating pair probabilities after calculating the partition function

//complexity==3 bactrack to calculate P matrices
void calculatePairsN3( DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE *Qm, DBL_TYPE *Qms,
                       DBL_TYPE *Qs, /*DBL_TYPE *Qn, DBL_TYPE *Qsn,*/ DBL_TYPE **Qx,
                       DBL_TYPE **Qx_1, DBL_TYPE **Qx_2, DBL_TYPE * Qb_bonus,
                       DBL_TYPE *P, DBL_TYPE *Pb, DBL_TYPE *Pm, DBL_TYPE *Pms,
                       DBL_TYPE *Ps, int seqlength,
                       int seq[], int *nicks, int** etaN);

void MakeP_Pm_N3( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qs,
                  DBL_TYPE *Qms, DBL_TYPE *Qm,
                  DBL_TYPE *P,  DBL_TYPE *Ps,
                  DBL_TYPE *Pms, DBL_TYPE *Pm,
                  int **etaN);
void MakePs_Pms( int i, int j, int seq[], int seqlength,
                 DBL_TYPE *Qs, DBL_TYPE *Qms, DBL_TYPE *Qb,
                 DBL_TYPE *Ps, DBL_TYPE *Pms, DBL_TYPE *Pb,
                 int *nicks, int **etaN);

//Consider Exterior loops
void prExterior_N3( int i,int j, int seq[], int seqlength,
                    DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE *Qb_bonus,
                    DBL_TYPE *P, DBL_TYPE *Pb,
                    int *nicks, int **etaN);

//Consider multiloops
void prMultiBp_N3( int i, int j, int seq[], int seqlength,
                   DBL_TYPE *Qb, DBL_TYPE *Qms, DBL_TYPE *Qm, DBL_TYPE *Qb_bonus,
                   DBL_TYPE *Pb, DBL_TYPE *Pms, DBL_TYPE *Pm, int **etaN);


//calculate contribution of interior loops to Pb
void prFastILoops( int i, int j, int L, int seqlength, int seq[],
                   DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2, 
                   DBL_TYPE *Qb_bonus,
                   DBL_TYPE *Pb, DBL_TYPE *Px, DBL_TYPE *Px_2,
                   int *nicks, int **etaN, float *preX, float *preX_2);

//calculate Pb contributions of small interior loops
void smallInteriorLoop( int pf_ij, int seq[], int seqlength, int i, int j,
                        int d, int e, int leftNick, int rightNick, DBL_TYPE *Qb,
                        DBL_TYPE * Qb_bonus,
                        DBL_TYPE *Pb, int error);

//manage Qx arrays during pairPr calculations
void prManageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1,
                 DBL_TYPE **Qx_2, DBL_TYPE **Px, DBL_TYPE **Px_1,
                 DBL_TYPE **Px_2, float **preX, float **preX_1,
                 float **preX_2, int len, int seqlength);

//complexity == 5 backtrack / P matrices calculations
void  calculatePairsN5( DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE *Qm,
                        DBL_TYPE *Qp, DBL_TYPE *Qz, DBL_TYPE *Qg,
                        DBL_TYPE *Qgl, DBL_TYPE *Qgr, DBL_TYPE *Qgls,
                        DBL_TYPE *Qgrs,
                        DBL_TYPE **QgIx, DBL_TYPE **QgIx_1, DBL_TYPE **QgIx_2,
                        DBL_TYPE *P,
                        DBL_TYPE *Pb, DBL_TYPE *Pp, DBL_TYPE *Pz,
                        DBL_TYPE *Pg, DBL_TYPE *Pbg,
                        DBL_TYPE *Pm, DBL_TYPE *Pgl, DBL_TYPE *Pgr,
                        DBL_TYPE *Pgls, DBL_TYPE *Pgrs,
                        int seqlength,
                        int seq[]);

//contributions from pseudoknots
void PseudoknotLoopN5( int i, int j, int pf_ij,
                       DBL_TYPE *Qp, DBL_TYPE *Qgl, DBL_TYPE *Qgr, DBL_TYPE *Qg,
                       DBL_TYPE *Qz,
                       DBL_TYPE *Pp, DBL_TYPE *Pgl, DBL_TYPE *Pgr, DBL_TYPE *Pg,
                       DBL_TYPE *Pz, DBL_TYPE *Pbg,
                       int seq[], int seqlength);

/*The functions with Make in their title, closely mimic the analogous MakeQ functions,
  but are used to calculate the P matrices */

void MakeP_Pm_Pz( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Q, DBL_TYPE *Qm, DBL_TYPE *Qz,
                  DBL_TYPE *Qb, DBL_TYPE *Qp, DBL_TYPE *P,
                  DBL_TYPE *Pm, DBL_TYPE *Pz, DBL_TYPE *Pb,
                  DBL_TYPE *Pp);

void MakePgr( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Qgr, DBL_TYPE *Qgl, DBL_TYPE *Qz,
              DBL_TYPE *Pgr, DBL_TYPE *Pgl, DBL_TYPE *Pz);

void MakePgl( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Qg, DBL_TYPE *Qgl, DBL_TYPE *Qz,
              DBL_TYPE *Pg, DBL_TYPE *Pgl, DBL_TYPE *Pz,
              DBL_TYPE *Pbg);

void MakePgrs( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm, DBL_TYPE *Qgrs, DBL_TYPE *Pg, DBL_TYPE *Pm,
               DBL_TYPE *Pgrs);

void MakePgls( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
               DBL_TYPE *Qm,  DBL_TYPE *Qgls, DBL_TYPE *Pg, DBL_TYPE *Pm,
               DBL_TYPE *Pgls);

void MakePg_N5( int i, int j, int seq[], int seqlength, DBL_TYPE *Qg,
                DBL_TYPE *Qm, DBL_TYPE *Qgls, DBL_TYPE *Qgrs,
                DBL_TYPE *QgIx, DBL_TYPE *QgIx_2,
                DBL_TYPE *Pg, DBL_TYPE *Pm, DBL_TYPE *Pgls, DBL_TYPE *Pgrs,
                DBL_TYPE *PgIx, DBL_TYPE *PgIx_2,
                float *preX, float *preX_2);

void MakePb_N5( int i, int j, int seq[], int seqlength,
                DBL_TYPE *Qm, DBL_TYPE *Qb, DBL_TYPE *Qp,
                DBL_TYPE *Pm, DBL_TYPE *Pb, DBL_TYPE *Pp);



//fast interior loop treatment of gap matrices
void prFastILoopsN5( int i, int j, int seq[], int seqlength,
                     DBL_TYPE *Qg, DBL_TYPE *QgIx, DBL_TYPE *QgIx_2,
                     DBL_TYPE *Pg,
                     DBL_TYPE *PgIx, DBL_TYPE *PgIx_2,
                     float *preX, float *preX_2);

void prManageQgIx( DBL_TYPE **QgIx, DBL_TYPE **QgIx_1,
                   DBL_TYPE **QgIx_2, DBL_TYPE **PgIx, DBL_TYPE **PgIx_1,
                   DBL_TYPE **PgIx_2,
                   float **preX, float **preX_1, float **preX_2,
                   int d, int seqlength);

void MakePg_Inextensible( int i, int j,int d, int e, int seq[], int seqlength,
                          DBL_TYPE *Qg, DBL_TYPE *Pg);


/* For pairPr calculations, the fast interior loop routines involve subtracting
   flaoting point numbers.  This can introduce precision errors, so the following two
   functions help keep track of this. */

//subtractLongDouble sets a = a - b, and returns the bits of precision that are lost
float subtractLongDouble( DBL_TYPE *a, DBL_TYPE b);

//The following two functions will recalculate Qx, QgIx from scratch, whenever
//the amount of error accumulated from subtractions exceeds a threshold (MAXPRECERR)
void recalculateQx( int i, int j, int size, int fbix, int seq[],
                    int seqlength, DBL_TYPE *Qx, DBL_TYPE *Qb,
                    int *nicks, int **etaN,
                    int side);
void recalculateQgIx( int i, int j, int d, int e, int size, int qgix, int seq[],
                      int seqlength, DBL_TYPE *QgIx, DBL_TYPE *Qb,
                      int side);


//Calculate interior loop probabilities the slow, N^4 way.
void prInteriorLoopsN4MS(int i, int j, int seq[], int seqlength,
                         DBL_TYPE *Qb, DBL_TYPE *Pb, int *nicks);
//old code used with prInteriorLoopN4MS.  There is no need to use this
void findNicks( int *nicks, int *leftNickIndex, int *rightNickIndex,
                int *nNicks, int leftEdge, int rightEdge);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
