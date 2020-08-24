#ifndef NUPACK_COMPLEXES_COMPLEXESUTILS_H__
#define NUPACK_COMPLEXES_COMPLEXESUTILS_H__

#include "complexesStructs.h"
#include <shared.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int getPermutation(int, int, int*);
void resetNicks(int, int*);
void nextMultiset(int, int*, int*, int*, int); //generate a multiset
int isCyclicP(int, int*, int*); //check if a cyclic permutation
void symmetryCheck(multiset*, int, permutation*); //check if symmetry
void printPerms(FILE*, int, int, multiset*);
void printMfesToFile(const dnaStructures *ds, FILE *fp,
                       const int *nicks);
int compareMultisets(const void*, const void*);
int comparePermutations(const void *, const void *);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_COMPLEXES_COMPLEXESUTILS_H__ */
