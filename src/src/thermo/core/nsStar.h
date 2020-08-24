#ifndef PFUNC_NSSTAR_H__
#define PFUNC_NSSTAR_H__

#include "shared.h"
#include "pf.h"

#ifdef __cplusplus
extern "C" {
#endif

// Computer n(s*) from a structure specified in pairs or parens format
// Either pairs or parents should be NULL for the function to work properly
DBL_TYPE nsStarPairsOrParensFull( int seqlength, int seq[], int *pairs,
                                  char *parens, int complexity, int naType, int dangles,
                                  DBL_TYPE temperature, DBL_TYPE sodiumconc,
                                  DBL_TYPE magnesiumconc, int uselongsalt);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
