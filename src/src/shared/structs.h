#ifndef NUPACK_SHARED_STRUCTS_H
#define NUPACK_SHARED_STRUCTS_H

#include "constants.h"

#ifdef __cplusplus
extern "C" {
#endif

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
enum parameter_set { DNA, RNA, RNA37, USE_SPECIFIED_PARAMETERS_FILE, COUNT};
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


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // NUPACK_SHARED_STRUCTS_H
