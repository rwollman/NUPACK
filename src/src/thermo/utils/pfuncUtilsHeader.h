/*
  pfuncUtilsHeader.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 7/2006, Justin Bois 1/2007

  Note that DBL_TYPE is defined in pfuncUtilsConstants.h, and
  indicates what type of floating point variables should be used
  (float, double, long double)

  See below for descriptions of each function
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif



#include "pfuncUtilsConstants.h"
#include "runtime_constants.h"
#include "physical_constants.h"
#include "../../shared/utilsHeader.h"

/* ************************ */
// Macros

// C compiler needs this! DO NOT REMOVE!  JZ
#ifndef MIN
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#endif
//This macro returns the min of x,y, regardless of type
#ifndef MAX
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#endif

/* ****************************** */

//fold struct describes a secondary structure and sequence
typedef struct{
  int *seq; //the sequence
  int seqlength; //sequence length
  int *pairs; //array indicating what is paired with what (see nsStarPairsOrParens)
  int *pknots; //similar to pairs, but indicates ends of pseudoknots
  int *fixedBases; //this is used in design code only, and restricts the identity of a postion
  int  *isNicked; //indicates if a position is right before the end of a strand
  int nStrands; //number of strands in a multi-stranded complex
} fold;

//paramter sets.  DNA, RNA are mfold 2.3 parameters that allow for temperature dependence
//RNA37 is the mfold3.0 parameter set which is only good at 37C.  COUNT sets energies to zero,
//so that the "partition function" is simply a count of the number of possible structures
enum { DNA, RNA, RNA37, USE_SPECIFIED_PARAMETERS_FILE, COUNT};
enum { FALSE, TRUE};

//oneDnaStruct and dnaStructures are used for enumerating sequences
typedef struct {
  int *theStruct; //describes what is paired to what
  DBL_TYPE error; //accumulated error (from the mfe) for a structure
  DBL_TYPE correctedEnergy; //actual energy of a structure
  int slength;
  //(accounting for symmetry).
} oneDnaStruct;

typedef struct {
  oneDnaStruct *validStructs;
  int nStructs; //# of structures stored
  int nAlloc; //# of structures allocated
  int seqlength;
  DBL_TYPE minError; //minimum deviation from mfe for all seqs
  //in validStructs

} dnaStructures;


// From utils.c in shared directory:
  /*
     Calculates the number of moles of water per liter at temperature
     (T in degrees C).
  */
double WaterDensity(double T);


/* ******************************************************************************** */
/*                          BEGIN FUNCTIONS FROM PF.C                               */
/* ******************************************************************************** */
/*
   pfuncFull: Calculates the partition function.
   Arguments:
   InputSeq is the sequence, including '+' characters to demarcate strand breaks.

   Complexity = 3 is the multi-stranded, pseudoknot-free method  [ O(N^3)]
   Complexity = 5 is the single-stranded, pseudoknot algorithm   [ O(N^5)]

   naType is chosen from the enum list given above, and chooses the parameter set

   dangles = 0, no dangle energies
   dangles = 1, mfold treatment of dangles
   dangles = 2, dangle energies are always summed, regardless of nearby structures
   (same as the Vienna package -D2 option)
   temperature = the temperature in celsius
   calcPairs = 1 indicates that pair probabilities are calculated and stored in global pairPr
   calcPairs = 0 skips pair probability calculations

   Ignores the possibility of symmetry
*/

DBL_TYPE pfuncFull( int inputSeq[], int complexity, int naType, int dangles,
                    DBL_TYPE temperature, int calcPairs,
                    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//pfuncFullWithSym is the Same as pfuncFull, but divides
//the result by permSym to account for symmetries
DBL_TYPE pfuncFullWithSym( int inputSeq[], int complexity, int naType,
                           int dangles, DBL_TYPE temperature, int calcPairs, int permSymmetry,
                           DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

DBL_TYPE pfuncFullWithSymHelper( int inputSeq[], int seqlength, int nStrands,
                                 int complexity, int naType,
                                 int dangles, DBL_TYPE temperature, int calcPairs,
                                 int permSymmetry, DBL_TYPE sodiumconc,
                                 DBL_TYPE magnesiumconc, int uselongsalt);


/* pfunc
   Calls pfuncFull, and assumes complexity = 3, DNA parameters, T = 37, dangles = 1,
   calcPairs = 1, [Na+] = 1.0, [Mg++] = 0.0, and short helix model for salt correction
*/
DBL_TYPE pfunc( int seq[]);


/* ******************************************************************************** */
/*                            END FUNCTIONS FROM PF.C                               */
/* ******************************************************************************** */



/* ******************************************************************************** */
/*                     BEGIN FUNCTIONS FROM MFEUTILS.C                              */
/* ******************************************************************************** */

/* mfeFull() returns an integer array (user allocated) describing the
   minimum free energy structure, with
   base i paired to base j <=> pairs[i] = j && pairs[j] = i;
   base i is unpaired <=> pairs[i].

   The arguments are the same as nsStarFull */
DBL_TYPE mfeFull( int inputSeq[], int seqLen, int *thepairs, int complexity,
                  int naType, int dangles, DBL_TYPE temperature,
                  DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//mfe makes the same assumptions as nsStar
DBL_TYPE mfe( int seq[], int seqLength, int *thepairs);

/* mfeFullWithSym() is used when a strand permutation
   has a specified symmetry factor e.g. ABAB has a symmetry of 2.
   The algorithm first calculates the mfe structure assuming no symmetry (mfeFull).
   If the reported mfe has a symmetric secondary structure, then its corrected
   energy should is the reported value + kT log( symmetry).  Consequently, an
   enumeration is performed of all structures within at most kT log(symmetry) of the
   calculated mfe.  The interval may be less if a structure deviating by
   exactly one pair from the mfe has a smaller energy gap.  The true mfe, taking
   into account symmetry is then calculated and reported.
   The struct mfeStructures (allocated by the user)
   will contain the enumerated structures, with the best at the front of the list */
DBL_TYPE mfeFullWithSym( int inputSeq[], int seqLen,
              dnaStructures *mfeStructures, int complexity, int naType,
              int dangles, DBL_TYPE temperature, int symmetry, int onlyOne,
              DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//mfeFullWithSym_Subopt is similar to mfeFullWithSym, but enumerates all structures
//within range of the algorithmic mfe (by algorithmic, I mean ignoring symmetry)
DBL_TYPE mfeFullWithSym_SubOpt( int inputSeq[], int seqLen,
                               dnaStructures *mfeStructures, int complexity, int naType,
                               int dangles, DBL_TYPE temperature, int symmetry,
                                DBL_TYPE fixedSubOptRange,
                                int onlyOne, DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
                                int uselongsalt);


/* getStructure() will create a dot-parens notation (stored in structure, used allocated)
   for the structure described by pairs (see above).  Will be misleading for
   pseudoknotted structures */
void getStructure( int seqlength, const int *thepairs, char *structure);

//initialize a structure of type dnaStructures
void initMfeStructures( dnaStructures*, int);

//used by qsort to arrange enumerated structures by corrected energy.
int compareDnaStructs( const void *, const void *);

//used during validation runs to sort the output based on structure
//instead of based on floating point numbers.
int compareDnaStructsOutput(const void *, const void *);
/* ******************************************************************************** */
/*                       END FUNCTIONS FROM MFEUTILS.C                              */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                  BEGIN FUNCTIONS FROM CALCULATEENERGY.C                          */
/* ******************************************************************************** */

/* naEnergyPairsOrParensFull() Computes the energy of a structure
   given either an int array of pairs (thepairs), or a character array
   (parens) in dot-parens notation.  See nsStarPairsOrParensFull().
   Assumes no symmetry.  Include the bimolecular association penalty
*/
DBL_TYPE naEnergyPairsOrParensFull( int *thepairs, char *parens,
                                   int inputSeq[], int naType,
                                    int dangles, DBL_TYPE temperature,
                                    DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
                                    int uselongsalt);

//same as above, but allows a value for possible symmetry (v(pi) = symmetry of sequence permutation)
DBL_TYPE naEnergyPairsOrParensFullWithSym( int *thepairs, char *parens,
                                          int inputSeq[], int naType,
                                          int dangles, DBL_TYPE temperature,
                                           int possibleSymmetry,
                                           DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc,
                                           int uselongsalt);


//naEnergyFull() is the same as above, but obtains the secondary structure from prefix.fold
DBL_TYPE naEnergyFull( char prefix[], int inputSeq[], int naType,
                       int dangles, DBL_TYPE temperature,
                       DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//Same as naEnergyFull(), but allows possible symmetry (see above)
DBL_TYPE naEnergyFullWithSym( char prefix[], int inputSeq[], int naType,
                              int dangles, DBL_TYPE temperature, int possibleSymmetry,
                              DBL_TYPE sodiumconc, DBL_TYPE magnesiumconc, int uselongsalt);

//computes energy of secondary structure in file, using same assumptions as pfunc
DBL_TYPE naEnergy( char *prefix, int seq[]);
//same as above, but structure input in thepairs or parens
DBL_TYPE naEnergyPairsOrParens( int *thepairs, char *parens,
                                int inputSeq[]);

//Makefold() creates a structure of type fold from parens or thepairs.
void MakeFold( fold *thefold, int seqlength, int seq[], char *parens, int *thepairs);
/* ******************************************************************************** */
/*                    END FUNCTIONS FROM CALCULATEENERGY.C                          */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                         BEGIN FUNCTIONS FROM INIT.C                              */
/* ******************************************************************************** */
/* ReadSequence() reads in a nucleic acid sequence from a file
   and allocates *seq to store it */
void ReadSequence( int *seqlength, char **seq, char *filename );

/* getSequenceLength() counts the number of strands (stored as nStrands)  in a given sequence (seq)
   and returns the sequence length */
int getSequenceLength( char *seq, int *nStrands);

/* getSequenceLength() counts the number of strands (stored as nStrands)  in a given sequence (seq)
   and returns the sequence length */
int getSequenceLengthInt( int seq[], int *nStrands);

/* processMultiSeqence() copies input sequence (containing strandPLUS members) to seq, but without the
   strandPLUS members.  The location of the strand breaks, strandPLUS, are stored in the nicks array */
void processMultiSequence( int inputSeq[], int seqlength, int nStrands,
                           int seq[], int nicks[]);

//Allocates Q and sets the values to zero.
void InitLDoublesMatrix( DBL_TYPE **Q, int size, char name[]);

//Sets Q to all zero
void ClearLDoublesMatrix(DBL_TYPE **Q, int size, char name[]);

//Memory management for "fast" interior loops subroutine
void manageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1,
               DBL_TYPE **Qx_2, int len, int seqlength);
//manageFx is the mfe version
void manageFx( DBL_TYPE **Fx, DBL_TYPE **Fx_1,
               DBL_TYPE **Fx_2, int len, int seqlength);


//Memory management routines for QgIx ("fast" gap spanning array treatment)
void manageQgIx( DBL_TYPE **QgIx, DBL_TYPE **QgIx_1,
                 DBL_TYPE **QgIx_2, int d, int seqlength);
void manageFgIx( DBL_TYPE **FgIx, DBL_TYPE **FgIx_1,
                 DBL_TYPE **FgIx_2, int d, int seqlength);

// Salt correction
  DBL_TYPE computeSaltCorrection(DBL_TYPE sodiumConc, DBL_TYPE magenesiumConc,
                                 int useLongHelix);

//Load energy parameters.  Global variable DNARNACOUNT determines parameter set
void LoadEnergies( void);
void setParametersToZero( void);

//Set Q[ pf_index(i, i-1, seqlength)] = 1;
void nonZeroInit( DBL_TYPE Q[], int seq[], int seqlength);

/*InitEtaN() initializes the etaN array.  etaN[ EtaNIndex(i,j,seqlength)][0] is the
  number of nicks between i and  j (i,j in [0.5, 1.5, 2.5, ...]).
  etaN[ EtaNIndex(i,j,seqlength)][1] is the index of the leftmost nick
  in the region between i and j, i.e. nicks[ EtaNIndex...] is the position
  of the leftmost nick between i and j.
*/
void InitEtaN( int **etaN, const int *nicks, int seqlength);
int EtaNIndex_old( float i, float j, int seqlength);

/* These functions set the size of Qg and Fg matrices is correct
   between successive calls to pfunc or mfe.  They also call PrecomputeValuesN5(f)  */
void initPF( int seqlength);
void initMfe( int seqlength);

//precomputes the sizeTerm array, used in extending interior loops from
//size i to i+2
void PrecomputeValuesN5( int seqlength);
//PrecomputeValuesN5f is for the mfe algorithm
void PrecomputeValuesN5f( int seqlength);
/* ******************************************************************************** */
/*                           END FUNCTIONS FROM INIT.C                              */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                     BEGIN FUNCTIONS FROM PFUNCUTILS.C                            */
/* ******************************************************************************** */
//pf_index calculates the array index for a Q-type array
int pf_index_old( int i, int j, int N);
#define pf_index(i,j,N) ((j)==(i)-1?(N)*((N)+1)/2 + (i) : ((i)*(N)+(j)-(i)*(1+(i))/2))
#define pf_index_same(i,N) ((i)*(N)-(i)*((i)-1)/2)
#define EtaNIndex(i,j,N) pf_index((int)(i),(int)(j),N)
#define EtaNIndex_same(i,N) pf_index_same((int)(i),N)
//((int)(j)==(int)(i)-1?(int)((N)*((N)+1)/2 + (int)(i)) : ((int)(i)*(N)+(j)-(i)*(1+(i))/2))
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

//converts a pair to an index (For energy calculations)
int GetMismatchShift( int base1, int base2);

//Returns the pair type (similar to GetMismatchShift), except assumes
//a Watson Crick (non-Wobble) pair
int GetPairType( int b);

//returns TRUE if two bases can form a pair.
int CanPairOld( int i, int j);
#ifdef NOGU
#define CanPair(i,j) ((i)+(j)==5 ?TRUE : FALSE)
#else
#define CanPair(i,j) ((i)+(j)==5 || (i)+(j)==7?TRUE : FALSE)
#endif

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
/* ******************************************************************************** */
/*                       END FUNCTIONS FROM PFUNCUTILS.C                            */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                       BEGIN FUNCTIONS FROM ENERGY.C                              */
/* ******************************************************************************** */
//Nearest neighbor energy of two consecutive base pairs
DBL_TYPE HelixEnergy( int i, int j, int h, int m);

//interior mismatch energy
DBL_TYPE InteriorMM( char a, char b, char x, char y);

//hairpin energy
DBL_TYPE HairpinEnergy( int i, int j, int seq[] );

//interior loop energy
DBL_TYPE InteriorEnergy(  int i, int j, int h, int m, int seq[]);

//interior loop energy, but allows for the exclusion of the i-j mismatch
//contribution (for use with fast i loop calculations)
DBL_TYPE InteriorEnergyFull( int i, int j, int h, int m, int seq[], int);

//Calculate dangle energies, assuming the entire structure is known
//(and can hence accurately calculate wobble pair dangles)
DBL_TYPE DangleEnergyWithPairs( int i, int j, fold* thefold);

//Calculates the dangle energy of a subsequence (i,j), assuming
//i-1, j+1 are paired (unless near a strand break).
//DangleEnergy miscalculates the dangles of nearby wobble pairs
DBL_TYPE DangleEnergy( int i, int j, int seq[], int seqlength);

//Calculates exp(-(dangle energy)/RT) and
//exp( -(interior loop energy)/RT), respectively
extern unsigned int seqHash;
DBL_TYPE ExplDangle( int i, int j, int seq[], int seqlength);
DBL_TYPE ExplInternal( int i, int j, int h, int m, int seq[]);

//NickDangle calculates the dangle energy, taking into account the effects
//of strand breaks (nicks).  If hairpin == TRUE, then this region is a nicked hairpin
//and may be closed by a wobble pair
DBL_TYPE NickDangle(  int i, int j, const int *nicks, int **etaN, int hairpin,
                      int seq[], int seqlength);

/* Computes the energy of an exterior loop with no secondary structure,
   and returns either exp( -energy/RT) or simply energy*/
DBL_TYPE NickedEmptyQ( int i, int j, const int nicks[], int seq[],
                       int seqlength, int **etaN);
DBL_TYPE NickedEmptyF( int i, int j, const int nicks[], int seq[],
                       int seqlength, int **etaN);

// Lookup table for the function 1.75*kB*TEMP_K*LOG_FUNC( size/30.0)
DBL_TYPE sizeLog(int size);

// Lookup table for contraction part of prFastILoops
DBL_TYPE sizeEnergyLog(int size);

//Computes the contribution of asymmetry and size to a large interior loop
DBL_TYPE asymmetryEfn( int L1, int L2, int size);
/* ******************************************************************************** */
/*                         END FUNCTIONS FROM ENERGY.C                              */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                         BEGIN FUNCTIONS FROM SUMEXP.C                            */
/* ******************************************************************************** */
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
                 DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2);


//makeNewQx creates new "extensible" base cases for the interval i,j.
void makeNewQx( int i, int j, int seq[], int seqlength,
                int **etaN, DBL_TYPE Qb[], DBL_TYPE Qx[]);
//extendOldQx extends Qx for the i-1, j+1 case
void extendOldQx( int i, int j, int seqlength,
                  DBL_TYPE Qx[], DBL_TYPE Qx_2[]);

//Directly calculates the contribution of small interior loops
DBL_TYPE SumExpInextensibleIL( int i, int j, int seq[], int seqlength,
                               DBL_TYPE Qb[],  int **etaN);
/* ******************************************************************************** */
/*                           END FUNCTIONS FROM SUMEXP.C                            */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                        BEGIN FUNCTIONS FROM SUMEXP_PK.C                          */
/* ******************************************************************************** */
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
/* ******************************************************************************** */
/*                          END FUNCTIONS FROM SUMEXP_PK.C                          */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                          BEGIN FUNCTIONS FROM PAIRSPR.C                          */
/* ******************************************************************************** */
//These files deal with calculating pair probabilities after calculating the partition function

//complexity==3 bactrack to calculate P matrices
void calculatePairsN3( DBL_TYPE *Q, DBL_TYPE *Qb, DBL_TYPE *Qm, DBL_TYPE *Qms,
                       DBL_TYPE *Qs, /*DBL_TYPE *Qn, DBL_TYPE *Qsn,*/ DBL_TYPE **Qx,
                       DBL_TYPE **Qx_1, DBL_TYPE **Qx_2,
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
                    DBL_TYPE *Q, DBL_TYPE *Qb,
                    DBL_TYPE *P, DBL_TYPE *Pb,
                    int *nicks, int **etaN);

//Consider multiloops
void prMultiBp_N3( int i, int j, int seq[], int seqlength,
                   DBL_TYPE *Qb, DBL_TYPE *Qms, DBL_TYPE *Qm,
                   DBL_TYPE *Pb, DBL_TYPE *Pms, DBL_TYPE *Pm, int **etaN);


//calculate contribution of interior loops to Pb
void prFastILoops( int i, int j, int L, int seqlength, int seq[],
                   DBL_TYPE *Qb, DBL_TYPE *Qx, DBL_TYPE *Qx_2,
                   DBL_TYPE *Pb, DBL_TYPE *Px, DBL_TYPE *Px_2,
                   int *nicks, int **etaN, float *preX, float *preX_2);

//calculate Pb contributions of small interior loops
void smallInteriorLoop( int pf_ij, int seq[], int seqlength, int i, int j,
                        int d, int e, int leftNick, int rightNick, DBL_TYPE *Qb,
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
/* ******************************************************************************** */
/*                            END FUNCTIONS FROM PAIRSPR.C                          */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                            BEGIN FUNCTIONS FROM MIN.C                            */
/* ******************************************************************************** */
/*These functions are used in minimum free energy calculations,
  and closely mimic their partition function counterparts */

//Returns hairpin energy, unless nicked (returns NAD_INFINITY)
DBL_TYPE MinHairpin( int i, int j, int seq[], int seqlength, int **etaN);

//finds the minimum energy multiloop closed by i,j. (complexity = 3)
DBL_TYPE MinMultiloops( int i, int j, int seq[],
                        DBL_TYPE *Fms, DBL_TYPE *Fm, int seqlength,
                        int **etaN);

//finds minimum energy exterior loop
DBL_TYPE MinExteriorLoop( int i,int j, int seq[], int seqlength,
                          DBL_TYPE *F, int *nicks, int **etaN);

//finds the minimum interior or multiloop (complexity > 3)
DBL_TYPE MinInterior_Multi( int i, int j, int seq[], int seqlength,
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
/* ******************************************************************************** */
/*                              END FUNCTIONS FROM MIN.C                            */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                           BEGIN FUNCTIONS FROM PKNOTS.C                          */
/* ******************************************************************************** */
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
/* ******************************************************************************** */
/*                             END FUNCTIONS FROM PKNOTS.C                          */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                        BEGIN FUNCTIONS FROM BACKTRACK.C                          */
/* ******************************************************************************** */
/*
   backtrack functions are used to compute the mfe structures once the
   mfe energies are known.  Due to the possibility of symmetry in
   multi-stranded problems, sequences with a periodic repeat of
   sequences will need to enumerate structures within epsilon of the
   algorithmic mfe, in order to compute the true mfe (see SIAM review
   paper for description.)  The algorithm used to enumerate is similar
   in flavor to the Vienna Package's RNAsubopt.  As a consequence, it
   is possible for these backtracking algorithms to take an
   exponentially large amount of time, although in practice they run
   very quickly.
*/

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
/* ******************************************************************************** */
/*                          END FUNCTIONS FROM BACKTRACK.C                          */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                        BEGIN FUNCTIONS FROM GETENERGY.C                          */
/* ******************************************************************************** */
//These functions are used to calculate the energy of a structure (linear time algorithm)

//compute the energy of the structure described in *thefold
DBL_TYPE GetEnergy( fold *thefold);

//Calculate the energy of an exterior loop substructure
DBL_TYPE EnergyF( int start, int stop, fold *thefold);
//Energy of Paired substructure
DBL_TYPE EnergyFb( int start, int stop, fold *thefold);
//Energy of pseudoknotted substructure
DBL_TYPE EnergyPk( int i, int j, fold *thefold);
//Energy of gap spanning region
DBL_TYPE EnergyFg( int i, int d, int e, int j, fold *thefold);
//Energy of pseudoknot interior
DBL_TYPE EnergyFz( int start, int stop, fold *thefold);
/* ******************************************************************************** */
/*                          END FUNCTIONS FROM GETENERGY.C                          */
/* ******************************************************************************** */

/* ******************************************************************************** */
/*                          BEGIN FUNCTIONS FROM NSSTAR.C                           */
/* ******************************************************************************** */
// Computer n(s*) from a structure specified in pairs or parens format
// Either pairs or parents should be NULL for the function to work properly
DBL_TYPE nsStarPairsOrParensFull( int seqlength, int seq[], int *pairs,
                                  char *parens, int complexity, int naType, int dangles,
                                  DBL_TYPE temperature, DBL_TYPE sodiumconc,
                                  DBL_TYPE magnesiumconc, int uselongsalt);
/* ******************************************************************************** */
/*                           END FUNCTIONS FROM NSSTAR.C                            */
/* ******************************************************************************** */


/* ******************************************************************************** */
/*                 BEGIN FUNCTIONS FROM READCOMMANDLINENPK.C                        */
/* ******************************************************************************** */
//get all command line options, including a possible input file
int ReadCommandLineNPK(int nargs, char **args, char *inputfile);

//get help message
void DisplayHelp(void);
//print out the thermodynamic options
void PrintNupackThermoHelp(void);
//print out the standard options for NUPACK Utilities executables
void PrintNupackUtilitiesHelp(void);

//read in input file
int ReadInputFile( char*, char*, int*, float*, char*, int*);

//get input interactively
void getUserInput(char*, int*,  float*, char*);

//determine if a permutation has a cyclic symmetry
int calculateVPi( int *, int);

//print outputs to the screen
void printInputs( int, char**, const char *, int,  const float*, const char*, char*);

//Prints a header with copyright information
void header( int, char**, char*, char*);
/* ******************************************************************************** */
/*                   END FUNCTIONS FROM READCOMMANDLINENPK.C                        */
/* ******************************************************************************** */

#ifdef __cplusplus
}
#endif

#endif

