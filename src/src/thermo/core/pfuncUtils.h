#ifndef __PFUNCUTILS_H__
#define __PFUNCUTILS_H__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <ctype.h>

#include "shared.h"
#include "mfeUtils.h"

#ifdef __cplusplus
extern "C" {
#endif


//pf_index calculates the array index for a Q-type array
int pf_index_old( int i, int j, int N);

#define pf_index(i,j,N) ((j)==(i)-1?(N)*((N)+1)/2 + (i) : ((i)*(N)+(j)-(i)*(1+(i))/2))

#define pf_index_same(i,N) ((i)*(N)-(i)*((i)-1)/2)

#define EtaNIndex(i,j,N) pf_index((int)(i),(int)(j),N)

#define EtaNIndex_same(i,N) pf_index_same((int)(i),N)
//((int)(j)==(int)(i)-1?(int)((N)*((N)+1)/2 + (int)(i)) : ((int)(i)*(N)+(j)-(i)*(1+(i))/2))

//converts a pair to an index (For energy calculations)
int GetMismatchShift( int base1, int base2);

//Returns the pair type (similar to GetMismatchShift), except assumes
//a Watson Crick (non-Wobble) pair
int GetPairType( int b);

//returns TRUE if two bases can form a pair.
int CanPairOld( int i, int j);

#ifdef NOGU
#define CanPair(i,j) ((i) + (j) == 5 ? TRUE : FALSE) 
#else
#define CanPair(i,j) ((i) + (j) == 5 || (i) + (j) == 7 ? TRUE : FALSE)
#endif
#define CanWCPair(i, j) ((i) + (j) == 5 ? TRUE : FALSE)

//gap_index calculates the array index of a "gap" matrix.
int gap_index( int h, int r, int m, int s, int seqlength);

//fbixIndex computes the array index for a Qx/Fx array (fast i loops)
int fbixIndexOld( int d, int i, int size, int N );
#ifndef DEBUG
#define fbixIndex(d, i, size, N ) ((i)*((d)-1) + (size))
#else
#define fbixIndex(d, i, size, N ) fbixIndexOld(d, i, size, N )
#endif

//QgIxIndex computes the array index for a QgIx/FgIx array (N^5 fast i loops)
int QgIxIndex( int d, int i, int size, int h1, int m1, int N);

//WithinEps checks if x-tol1 < y < x+tol2
int WithinEps( DBL_TYPE x, DBL_TYPE y, DBL_TYPE tol1, DBL_TYPE tol2);

//converts A->1, C->2, G->3, T/U->4
int Base2int( char base);
char Int2base( int base);
// Converts sequence of letters to sequence of
// returns -1 on success, index of faulty base letter on error
int convertSeq(char *seqchar, int* seqnum, int seqlength);
int printSeqNum(int seqnum[]);

//A speed enhancement for the N^5 Algorithms.  Checks if there exists
//a possible pair between (i-b,j+b) where is a non-negative integer.
//Sets possiblePairs[ pf_index(i,j,seqlength)] == TRUE if possible
void CheckPossiblePairs( short **possiblePairs, int seqlength, int seq[]);

//LoadFold loads a structure of type fold from a data file.
void LoadFold( fold *thefold, char filename[]);

/* expectedCorrectBases uses the calculated pair probabilities (stored in pairPr)
   to determine, on average, how many bases are in the same state as the
   secondary structure described by the integer array pairs.
   base i paired to base j <=> pairs[i] = j && pairs[j] = i.
   Base i is unpaired <=> pairs[i];
*/
DBL_TYPE expectedCorrectBases( int *pairs, int seqlength);

/* checkSymmetry() determines whether a specified secondary structure
   (specifed by thepairs), has rotational symmetries.  Only symmetries that
   evenly divide possibleSymmetry are checked. */
int checkSymmetry( const int *thepairs, int seqlength,
                   const int *nicks, int possibleSymmetry,
                   int nStrands);

void findUniqueMins( dnaStructures*, const int*, int, int, DBL_TYPE);

int comparePairs( const int *p1, const int *p2,
                  const int seqlength, const int unitLength,
                  const int possibleSymmetry);

/* getStructureFromParens() computes an integer array of the pairs from the dot-parens array (line).
   Parens of the same type () or [] or {} or <> are assumed to be nested. Different types need not be
   nested with one another. */
void getStructureFromParens( char *line, int *pairs, int seqlength);


/* ******************************************************************************** */
//   functions for maintaining a list of mfe structures
//   (or structures within mfeEpsilon of the algorithmic mfe).  See above for a
//   description of the dnaStructures struct */
/* ******************************************************************************** */
//delete all sequences stored in dnaStructures
void clearDnaStructures( dnaStructures*);

//make the first argument a duplicate of the second (copy values)
void copyDnaStructures( dnaStructures*, const dnaStructures*);

//append the structures in the second argument to the first
//all the new structures deviate from the mfe by newErr, so only copy those
//whose current error + newErr are still less than the allowable mfeEpsilon
void addDnaStructures( dnaStructures*, const dnaStructures*, DBL_TYPE newErr,
                       DBL_TYPE mfeEpsilon, int only_min);

//Print a single secondary structure, as described by thepairs (see nsStar),
//including strand breaks as '+'.  Different symbols for pairs are cycled through
//as pseudoknots are introduced.
void PrintStructure( char *thefold, const int *thepairs, int **etaN,
                     int seqlength, char *filename);

//Print all structures saved in *ds, using PrintStructure
void PrintDnaStructures( const dnaStructures *ds, int **etaN, const int *nicks,
                         int symmetry, char *filename);

//A dumbed down version of PrintDnaStructures, but only uses . ( ), ignoring multistrands and
//pseudoknots.  Used only for debugging
void PrintS( const dnaStructures *ds);

//check if a file exists
int fileExists( char*);

// Print a matrix with the 
void print_pf_mat(DBL_TYPE *mat, int len, char *name, FILE *fname);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
