/** \file backtrack.h
 * backtrack functions are used to compute the mfe structures once the
 * mfe energies are known.  Due to the possibility of symmetry in
 * multi-stranded problems, sequences with a periodic repeat of
 * sequences will need to enumerate structures within epsilon of the
 * algorithmic mfe, in order to compute the true mfe (see SIAM review
 * paper for description.)  The algorithm used to enumerate is similar
 * in flavor to the Vienna Package's RNAsubopt.  As a consequence, it
 * is possible for these backtracking algorithms to take an
 * exponentially large amount of time, although in practice they run
 * very quickly.
 */

#ifndef NUPACK_THERMO_CORE_BACKTRACK_H__
#define NUPACK_THERMO_CORE_BACKTRACK_H__

#include "ene.h"
#include "min.h"
#include "pfuncUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

//The thepairs[i] = j, thepairs[j] = i.  If either is already set, return an error
void SetPair( int i, int j, int *thepairs);

//Call SetPair for all structures in the list held in *ds
void SetPairs( int i, int j, dnaStructures *ds);

//complexity = 3 backtrack
//These functions are all analogous to their partition function counterparts.
void bktrF_Fm_N3( int i, int j, int seq[], int seqlength,
                  const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                  const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks, int **etaN,
                  dnaStructures *dnaStr,
                  const char *type,
                  const int *maxILoopSize,
                  const DBL_TYPE mfeEpsilon,
                  const int onlyOne);
void bktrFs_Fms( int i, int j, int seq[], int seqlength,
                 const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                 const DBL_TYPE *Fs, const DBL_TYPE *Fms,
                 const int *nicks, int **etaN, dnaStructures *dnaStr, const char *type,
                 const int *maxILoopSize,
                 const DBL_TYPE mfeEpsilon,
                 const int onlyOne);

void bktrFb_N3( int i, int j, int seq[], int seqlength, const DBL_TYPE *F, const DBL_TYPE *Fb,
                const DBL_TYPE *Fm, const DBL_TYPE *Fs, const DBL_TYPE *Fms,
                const int *nicks, int **etaN, dnaStructures *dnaStr,
                const int *maxILoopSize,
                const DBL_TYPE mfeEpsilon,
                const int onlyOne);
int bktrMinMultiloops( int i, int j, int seq[], int seqlength,
                       const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                       const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks,
                       int **etaN, dnaStructures *dnaStr,
                       const int *maxILoopSize,
                       const DBL_TYPE mfeEpsilon,
                       const int onlyOne);
int bktrMinExteriorLoop( int i, int j, int seq[], int seqlength,
                         const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                         const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks,
                         int **etaN, dnaStructures *dnaStr,
                         const int *maxILoopSize,
                         const DBL_TYPE mfeEpsilon,
                         const int onlyOne);
int bktrMinInteriorLoop( int i, int j, int seq[], int seqlength,
                         const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                         const DBL_TYPE *Fs, const DBL_TYPE *Fms, const int *nicks,
                         int **etaN, dnaStructures *dnaStr,
                         const int *maxILoopSize,
                         const DBL_TYPE mfeEpsilon,
                         const int onlyOne);


//complexity = 5
void bktrF_Fm_FzN5( int i, int j, int seq[], int seqlength,
                    const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
                    const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
                    const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
                    const DBL_TYPE *Fgr, dnaStructures *dnaStr, const int *nicks,
                    int **etaN, DBL_TYPE mfeEpsilon, const char *type);

void bktrFbN5( int i, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr, dnaStructures *dnaStr,
               const int *nicks, int **etaN,
               const DBL_TYPE mfeEpsilon);

void bktrFgN5( int i, int d, int e, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr,  dnaStructures *dnaStr,
               const int *nicks, int **etaN,
               const DBL_TYPE mfeEpsilon);

void bktrFpN5( int i, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr, dnaStructures *dnaStr, const int *nicks,
               int **etaN, DBL_TYPE mfeEpsilon);

void bktrFgls( int i, int d, int e, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr, dnaStructures *dnaStr,
               const int *nicks, int **etaN,
               const DBL_TYPE mfeEpsilon);

void bktrFgrs( int i, int d, int e, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr, dnaStructures *dnaStr,
               const int *nicks, int **etaN,
               const DBL_TYPE mfeEpsilon);

void bktrFgl(  int i, int e, int f, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr, dnaStructures *dnaStr,
               const int *nicks, int **etaN,
               const DBL_TYPE mfeEpsilon);

void bktrFgr(  int i, int d, int e, int j, int seq[], int seqlength,
               const DBL_TYPE *F, const DBL_TYPE *Fb, const DBL_TYPE *Fm,
               const DBL_TYPE *Fp, const DBL_TYPE *Fz, const DBL_TYPE *Fg,
               const DBL_TYPE *Fgls, const DBL_TYPE *Fgrs, const DBL_TYPE *Fgl,
               const DBL_TYPE *Fgr, dnaStructures *dnaStr,
               const int *nicks, int **etaN,
               const DBL_TYPE mfeEpsilon);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
