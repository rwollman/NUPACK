#ifndef PFUNC_PKNOTS_H__
#define PFUNC_PKNOTS_H__

#include <shared.h>
#include "ene.h"
#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

//Finds minimum energy structures for pseudknot matrices

//find minimum energy, rightmost pseudoknot contained within closing pair i,j
DBL_TYPE MinFb_Pk( int i, int j, int seq[], int seqlength, DBL_TYPE *Fp,
                   DBL_TYPE *Fm);

//Fill out minimum values for all Fg matrices with outer pair i,j
void MakeFg_N5( int i, int j, int seq[], int seqlength, DBL_TYPE *Fg,
                DBL_TYPE *Fm, DBL_TYPE *Fgls, DBL_TYPE *Fgrs, DBL_TYPE *FgIx,
                DBL_TYPE *FgIx_2, short *possiblePairs);

//These four functions compute interior loops in Fg matrices (complexity = 5)
void fastIloop_Fg(int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Fg, DBL_TYPE *FgIx, DBL_TYPE *FgIx_2,
                  short *possiblePairs);
DBL_TYPE MinInextensibleIL_Fg( int i, int j, int d, int e,
                               int seq[], int seqlength,
                               DBL_TYPE *Fg);
void makeNewFgIx( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *Fg, DBL_TYPE *FgIx);

void extendOldFgIx( int i, int j, int d, int e, int seq[], int seqlength,
                    DBL_TYPE *Fg, DBL_TYPE *FgIx, DBL_TYPE *FgIx_2);


//Find min values for Fgls, Fgrs, Fgl, Fgr, F, Fm, Fz, respectively
void MakeFgls( int i, int j, int seq[], int seqlength, DBL_TYPE *Fg,
               DBL_TYPE *Fm,  DBL_TYPE *Fgls);
void MakeFgrs( int i, int j, int seq[], int seqlength, DBL_TYPE *Fg,
               DBL_TYPE *Fm, DBL_TYPE *Fgrs);
void MakeFgl( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Fg, DBL_TYPE *Fgl, DBL_TYPE *Fz);
void MakeFgr( int i, int j, int seq[], int seqlength,
              DBL_TYPE *Fgr, DBL_TYPE *Fgl, DBL_TYPE *Fz);
void MakeF_Fm_Fz( int i, int j, int seq[], int seqlength,
                  DBL_TYPE *F, DBL_TYPE *Fm, DBL_TYPE *Fz,
                  DBL_TYPE *Fb, DBL_TYPE *Fp);

//Find the min energy pseudoknot with outer boundaries i and j
DBL_TYPE MinFp_N5( int i, int j, int seq[], int seqlength,
                   DBL_TYPE *Fgl, DBL_TYPE *Fgr, DBL_TYPE *Fg,
                   DBL_TYPE *Fz);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
