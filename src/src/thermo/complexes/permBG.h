#ifndef NUPACK_COMPLEXES_PERMBG_H__
#define NUPACK_COMPLEXES_PERMBG_H__

#include "complexesStructs.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void neg(int t, int n, int k);
void gen(int t, int n, int k);
void PrintIt(void);
void BigGen(int);
void swap(int, int, int);
void initializeMP(int, int*);
void freeMP(void);
int nextPerm(void);
void setPerm(int*);
// Count the number of sets in a list of permutations
int CountSets(permutation * perms, int nPerms, int nStrands);
// Get the maximum complex size
int GetMaxComplexSize(multiset * allSets, int totalSets);
// Fill in the set information using the permutations
int FillSets(multiset * allSets, permutation * allPermutations,
             int totalSets, int totalPerms,
             int nStrands, int * seqlength);
// Generate fixed content necklaces.
// There is a more efficient algorithm
// (see Sawada 2003) but we don't need this for the commonly used (and actually
// documented) use cases of the executable. The only time this is used is when a
// C line stands alone and we need to generate all necklaces for that complex.
int makeFCPermutations(permutation * loc, int * composition, int length, int nStrands);
// fill out circular permutations of length "length" and
// alphabet size "nStrands"
int makePermutations(permutation * loc, int length, int nStrands);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_COMPLEXES_PERMBG_H__ */
