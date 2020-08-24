/*
  OutputWriter.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 1/2007

  For use with Distributions.c.

  Writes the output to the appropriate output files.

  For input and output formatting, see associated manual.
*/


#include "OutputWriter.h"
  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"


// Structure for sorting output that includes permutations
struct PermSortStruct {
  int CompID; // The ID tag associated with a complex
  int PermID; // ID number for permutation
  int *Aj; // Array representing column j of A
  double FreeEnergy; // Free energy for permutation
  double mj; // The concentration of the permutation
  double mjc; // The complex concentration
  double *Pmn; // The probability distribution
  int numSS; // Number of entries in Aj (needed for sorting routine)
  char *AuxStr; // String containing auxillary information from input file
};

// Structure for sorting output that has only complexes
struct CompSortStruct {
  int CompID; // The ID tag associated with a complex
  double FreeEnergy; // Free energy for permutation
  int *Aj; // Array representing column j of A
  double mj; // The concentration of the complex
  double *Pmn; // The probability distribution
  int numSS; // Number of entries in Aj (needed for sorting routine)
  char *AuxStr; // String containing auxillary information from input file
};

/* ******************************************************************************** */
void WriteOutput(double *mEq, double **Pmn, double *G, int **A, int *CompIDArray, 
                 int *PermIDArray, int LargestCompID, int numSS, int numTotal, 
                 int nTotal, int nComments, int maxm0, double kT, char *cxFile, 
                 int SortOutput, char *distFile, int quiet, int NoPermID, int NUPACK_VALIDATE) {
 /*
   Writes the output for a distributions calculation.
     */
 
 // Go back through the input file and load in the permutation information.
 int i,j,k,n; // Counters
 char line[MAXLINE]; // A line from a file
 char *tok; // Token
 char tokseps[] = " \t\n"; // Token separators  
 int *CompLookup; // Lookup array matching complex ID's to indices in x, m, G, etc.
 int ComplexID; // Complex ID for complex read in from cx file
 double *mEqc; // Array of complex equlibrium counts
 char **AuxStrArray; // Array of auxillary strings
 struct CompSortStruct *CompOutStruct; // Output structure for complexes only
 struct PermSortStruct *PermOutStruct; // Output structure for sorting with perms
 FILE *fp; // File handle for the cx file
 FILE *fpdist; // File handles for the output files
 
 
 // Allocate memory for lookup table for complexes (to check if conc. is nonzero)
 CompLookup = (int *) malloc (LargestCompID * sizeof(int));
 // Initialize it to -1
 for (j = 0; j < LargestCompID; j++) {
   CompLookup[j] = -1;
 }
 // If included in problem change entry
 for (j = 0; j < numTotal; j++) {
   CompLookup[CompIDArray[j]-1] = j;
 }
 
 // Open the cx file
 if ((fp = fopen(cxFile,"r")) == NULL) {
   if (quiet == 0) {
     printf("Error in opening file %s!\n",cxFile);
     printf("\nExiting....\n\n");
   }
   exit(ERR_CX);
 }
 
 // Blow through comments
 for (j = 0; j < nComments; j++) {
   fgets(line,MAXLINE,fp);
 }
 
 if (NoPermID == 1) { // Only report complexes
   
   // Allocate memory for output structure
   CompOutStruct = (struct CompSortStruct *) 
     malloc(LargestCompID * sizeof(struct CompSortStruct));
   for (j = 0; j < LargestCompID; j++) {
     CompOutStruct[j].Aj = (int *) malloc(numSS * sizeof(struct CompSortStruct));
     CompOutStruct[j].Pmn = (double *) malloc((maxm0+1) * sizeof(double));
   }

   for (k = 0; k < nTotal; k++) {
     // Read in the line from the cx file
     while (fgets(line,MAXLINE,fp) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n'));
     // Get the complex ID
     tok = strtok(line,tokseps);
     ComplexID = atoi(tok) - 1;
     CompOutStruct[ComplexID].CompID = ComplexID + 1;
     
     // Pull out column of A corresponding to complex
     for (i = 0; i < numSS; i++) {
       tok = strtok(NULL,tokseps);
       CompOutStruct[ComplexID].Aj[i] = atoi(tok);
     }
     
     // Enter the free energy
     tok = strtok(NULL,tokseps);
     CompOutStruct[ComplexID].FreeEnergy = str2double(tok) / kT;
     
     // The rest of the auxillary information
     if ((tok = strtok(NULL,"")) != NULL) {
       CompOutStruct[k].AuxStr = (char *) malloc((strlen(tok)+1) * sizeof(char));
       strcpy(CompOutStruct[k].AuxStr,tok);
       // Delete null character if there is one
       if (CompOutStruct[k].AuxStr[strlen(CompOutStruct[k].AuxStr)-1] == '\n') {
         CompOutStruct[k].AuxStr[strlen(CompOutStruct[k].AuxStr)-1] = '\0';
       }
     }
     else {
       CompOutStruct[k].AuxStr = (char *) malloc(1 * sizeof(char));
       CompOutStruct[k].AuxStr[0] = '\0';
     }
     
     // Put numSS in just because we have to for qsort
     CompOutStruct[ComplexID].numSS = numSS;
     
     // mj and Pmn (this is done redundantly)
     if (CompLookup[ComplexID] == -1) {
       CompOutStruct[ComplexID].mj = 0;
       CompOutStruct[ComplexID].Pmn[0] = 1;
       for (n = 1; n <= maxm0; n++) {
         CompOutStruct[ComplexID].Pmn[n] = 0;
       }
     }
     else {
       CompOutStruct[ComplexID].mj = mEq[CompLookup[ComplexID]];
       for (n = 0; n <= maxm0; n++) {
         CompOutStruct[ComplexID].Pmn[n] = Pmn[CompLookup[ComplexID]][n];
       }
     }
   }
   
   // close the cx file
   fclose(fp);
   
   // Sort the results
   if (SortOutput == 1) {
     qsort(CompOutStruct,LargestCompID,sizeof(struct CompSortStruct),Compare25);
   }
   else if (SortOutput == 3) {
     qsort(CompOutStruct,LargestCompID,sizeof(struct CompSortStruct),Compare26);
   }
   else { // SortOutput == 4
     qsort(CompOutStruct,LargestCompID,sizeof(struct CompSortStruct),Compare27);
   }
   
   // Write out results
   if ((fpdist = fopen(distFile,"a")) == NULL) {
     if (quiet == 0) {
       printf("Error opening %s.\n\nExiting....\n",distFile);
     }
     exit(ERR_DIST);
   }
   for (j = 0; j < LargestCompID; j++) {
     // ComplexID and dummy permutation number (-1).
     fprintf(fpdist,"%d\t",CompOutStruct[j].CompID);
     if (NoPermID == 0) {
       fprintf(fpdist,"-1\t"); // Dummy perm ID
     }
     
     // The corresponding column in A
     for (i = 0; i < numSS; i++) {
       fprintf(fpdist,"%d\t",CompOutStruct[j].Aj[i]);
     }
     
     // Free energy in kcal/mol
     if(!NUPACK_VALIDATE) {
        fprintf(fpdist,"%8.6e\t",CompOutStruct[j].FreeEnergy * kT);
        // Expected count
        fprintf(fpdist,"%8.6e\t",CompOutStruct[j].mj);
        // Probability distribution for counts
        for (n = 0; n < maxm0+1; n++) {
          fprintf(fpdist,"%8.6e\t",CompOutStruct[j].Pmn[n]);
        }
     } else {
        fprintf(fpdist,"%.14e\t",CompOutStruct[j].FreeEnergy * kT);
        // Expected count
        fprintf(fpdist,"%.14e\t",CompOutStruct[j].mj);
        // Probability distribution for counts
        for (n = 0; n < maxm0+1; n++) {
          fprintf(fpdist,"%.14e\t",CompOutStruct[j].Pmn[n]);
        }
      }

     
     if (NoPermID == 1) { // Can still have auxillary information
       fprintf(fpdist,"%s\n",CompOutStruct[j].AuxStr);
     }
   }
   fclose(fpdist);
   
   
   // Free the allocated memory
   for (j = 0; j < LargestCompID; j++) {
     free(CompOutStruct[j].Aj);
     free(CompOutStruct[j].Pmn);
     if (NoPermID == 1) {
       free(CompOutStruct[j].AuxStr);
     }
   }
   free(CompOutStruct);
   
 }
 
 else { // Sorting done on permutations
   
   // Allocate memory for the output structure
   PermOutStruct = (struct PermSortStruct *) 
     malloc(nTotal * sizeof(struct PermSortStruct));
   AuxStrArray = (char **) malloc(nTotal * sizeof(char *));
   for (k = 0; k < nTotal; k++) {
     PermOutStruct[k].Aj = (int *) malloc(numSS * sizeof(struct PermSortStruct));
     PermOutStruct[k].Pmn = (double *) malloc((maxm0+1) * sizeof(double));
   }
   
   // Allocate memory for storage of mEqc and initialize
   mEqc = (double *) malloc(LargestCompID * sizeof(double));
   for (j = 0; j < LargestCompID; j++) {
     mEqc[j] = 0;
   }
   
   // Read in the data from the ocx file
   k = 0; // Counter in PermOutStruct
   j = 0; // Counter in AuxStrArray
   while (fgets(line,MAXLINE,fp) != NULL && line[0] != '\0' && line[0] != '\n') {
     
     // Get the complex ID
     tok = strtok(line,tokseps);
     ComplexID = atoi(tok) - 1;
     
     if (CompLookup[ComplexID] == -1) { // Not included in calculation
       // Permutation number
       tok = strtok(NULL,tokseps);
       PermOutStruct[k].PermID = atoi(tok);
       
       // Pull out column of A corresponding to complex
       for (i = 0; i < numSS; i++) {
         tok = strtok(NULL,tokseps);
         PermOutStruct[k].Aj[i] = atoi(tok);
       }
       
       // Enter the free energy
       tok = strtok(NULL,tokseps);
       PermOutStruct[k].FreeEnergy = str2double(tok) / kT;
       
       // The rest of the auxillary information
       if ((tok = strtok(NULL,"")) != NULL) {
         PermOutStruct[k].AuxStr = (char *) malloc((strlen(tok)+1) * sizeof(char));
         strcpy(PermOutStruct[k].AuxStr,tok);
         // Delete newline character if there is one
         if (PermOutStruct[k].AuxStr[strlen(PermOutStruct[k].AuxStr)-1] == '\n') {
           PermOutStruct[k].AuxStr[strlen(PermOutStruct[k].AuxStr)-1] = '\0';
         }
       }
       else {
         PermOutStruct[k].AuxStr = (char *) malloc(1 * sizeof(char));
         PermOutStruct[k].AuxStr[0] = '\0';
       }
       
       // Put numSS in just because we have to for qsort
       PermOutStruct[k].numSS = numSS;
       PermOutStruct[k].CompID = ComplexID + 1;
       PermOutStruct[k].mj = 0;
       PermOutStruct[k].Pmn[0] = 1;
       for (n = 1; n <= maxm0; n++) {
         PermOutStruct[k].Pmn[n] = 0;
       }
       // Advance counter
       k++;
     }
     else {  // Included in calculation
       // Permutation number
       tok = strtok(NULL,tokseps);
       
       // column of A corresponding to complex
       for (i = 0; i < numSS; i++) {
         tok = strtok(NULL,tokseps);
       }
       
       // Enter the free energy
       tok = strtok(NULL,tokseps);
       
       // The rest of the auxillary information
       if ((tok = strtok(NULL,"")) != NULL) {
         AuxStrArray[j] = (char *) malloc((strlen(tok)+1) * sizeof(char));
         strcpy(AuxStrArray[j],tok);
         // Delete newline character if there is one
         if (AuxStrArray[j][strlen(AuxStrArray[j])] == '\n') {
           AuxStrArray[j][strlen(AuxStrArray[j])] = '\0';
         }
       }
       else {
         AuxStrArray[j] = (char *) malloc(1 * sizeof(char));
         AuxStrArray[j][0] = '\0';
       }
       j++;
     }
   }
   
   // Enter values for entries included in calculation
   for (j = 0; j < nTotal-k; j++) {
     PermOutStruct[j+k].CompID = CompIDArray[j];
     PermOutStruct[j+k].PermID = PermIDArray[j];
     for (i = 0; i < numSS; i++) {
       PermOutStruct[j+k].Aj[i] = A[i][j];
     }
     PermOutStruct[j+k].FreeEnergy = G[j];
     PermOutStruct[j+k].numSS = numSS;
     PermOutStruct[j+k].mj = mEq[j];
     mEqc[CompIDArray[j]-1] += mEq[j];
     for (n = 0; n <= maxm0; n++) {
       PermOutStruct[j+k].Pmn[n] = Pmn[j][n];
     }
     PermOutStruct[j+k].AuxStr = 
       (char *) malloc((strlen(AuxStrArray[j])+1) * sizeof(char));
     strcpy(PermOutStruct[j+k].AuxStr,AuxStrArray[j]);
     free(AuxStrArray[j]);
   }
   
   // Free the AuxStrArray
   free(AuxStrArray);
   
   // Enter mjc entries
   for (k = 0; k < nTotal; k++) {
     PermOutStruct[k].mjc = mEqc[PermOutStruct[k].CompID-1];
   }
   
   // close the cx file
   fclose(fp);
   
   // Sort the results (no sorting necessary for SortOutput == 0)
   if (SortOutput == 1) {
     qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare21);
   }
   else if (SortOutput == 2) {
     qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare22);
   }
   else if (SortOutput == 3) {
     qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare23);
   }
   else if (SortOutput == 4) {
     qsort(PermOutStruct,nTotal,sizeof(struct PermSortStruct),Compare24);
   }
   
   
   // Write out results
   if ((fpdist = fopen(distFile,"a")) == NULL) {
     if (quiet == 0) {
       printf("Error opening %s.\n\nExiting....\n",distFile);
     }
     exit(ERR_DIST);
   }
   for (k = 0; k < nTotal; k++) {
     // ComplexID and permutation number.
     fprintf(fpdist,"%d\t%d\t",PermOutStruct[k].CompID,PermOutStruct[k].PermID);
     
     // The corresponding column in A
     for (i = 0; i < numSS; i++) {
       fprintf(fpdist,"%d\t",PermOutStruct[k].Aj[i]);
     }
     
     if(!NUPACK_VALIDATE) {
      // Free energy in kcal/mol
      fprintf(fpdist,"%8.6e\t",PermOutStruct[k].FreeEnergy * kT);
      
      // Expected value of count
      fprintf(fpdist,"%8.6e\t",PermOutStruct[k].mj);
      
      // Probability distribution for counts
      for (n = 0; n < maxm0 + 1; n++) {
        fprintf(fpdist,"%8.6e\t",PermOutStruct[k].Pmn[n]);
      }
     } else {
      fprintf(fpdist,"%.14e\t",PermOutStruct[k].FreeEnergy * kT);
      
      // Expected value of count
      fprintf(fpdist,"%.14e\t",PermOutStruct[k].mj);
      
      // Probability distribution for counts
      for (n = 0; n < maxm0 + 1; n++) {
        fprintf(fpdist,"%.14e\t",PermOutStruct[k].Pmn[n]);
      }
     }
     
     // Auxillary information
     fprintf(fpdist,"%s\n",PermOutStruct[k].AuxStr);
   }
   fclose(fpdist);
   
   // Free the allocated memory
   for (k = 0; k < nTotal; k++) {
     free(PermOutStruct[k].Aj);
     free(PermOutStruct[k].Pmn);
     free(PermOutStruct[k].AuxStr);
   }
   free(PermOutStruct);
   free(mEqc);
   
 }
 
 // Free CompLookup
 free(CompLookup);
 
}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare21(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Whichever has highest expected count (<mj>) returns -1.

     If the concentrations are the same (usually only the case if they're both zero),
     they are sorted by complex ID and then permutation ID.
  */

  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;

  if (ps1->mj < ps2->mj) {
    return 1;
  }
  else if (ps1->mj > ps2->mj) {
    return -1;
  }
  else { // Equal permutation concentration
    if (ps1->CompID < ps2->CompID) {
      return -1;
    }
    else if (ps1->CompID > ps2->CompID) {
      return 1;
    }
    else {  // same complex ID
      if (ps1->PermID < ps2->PermID) {
	return -1;
      }
      else if (ps1->PermID > ps2->PermID) {
	return 1;
      }
      else { // Shouldn't ever get here
	return 0;
      }
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare22(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Used for sorting first by complex concentration then permutation concentration
  */

  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;


  if (ps1->CompID == ps2->CompID) { // Same complex, sort by perm conc.
    if (ps1->mj < ps2->mj) {
      return 1;
    }
    else if (ps1->mj > ps2->mj) {
      return -1;
    }
    else { // Same permutation concentration (sort by Perm ID)
      if (ps1->PermID < ps2->PermID) {
	return -1;
      }
      else if (ps1-> PermID > ps2->PermID) {
	return 1;
      }
      else { // Same PermID
	return 0;
      }
    }
  }
  else { // Different complexes, sort by complex conc.
    if (ps1->mjc < ps2->mjc) {
      return 1;
    }
    else if (ps1->mjc > ps2->mjc) {
      return -1;
    }
    else { // Same complex concentration (sort by complex ID)
      if (ps1->CompID < ps2->CompID) {
	return -1;
      }
      else if (ps1-> CompID > ps2->CompID) {
	return 1;
      }
      else { // Shouldn't ever get here
	return 0;
      }
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare23(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     They are sorted by complex ID and then permutation ID.
  */

  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;

  if (ps1->CompID < ps2->CompID) {
    return -1;
  }
  else if (ps1->CompID > ps2->CompID) {
    return 1;
  }
  else {  // same complex ID
    if (ps1->PermID < ps2->PermID) {
      return -1;
    }
    else if (ps1->PermID > ps2->PermID) {
      return 1;
    }
    else { // Shouldn't ever get here
      return 0;
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare24(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Used for sorting first by number of strands in complex, and then
     alphabetically.  I.e., an order of complexes might be:
     A, B, AA, AB, BB

     Then it's sorted by permutation ID number.
  */

  int i; // Counter
  const struct PermSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermSortStruct *ps2 = p2;

  if (ps1->CompID == ps2->CompID) {  // Same complex
    if (ps1->PermID < ps2->PermID) {
      return -1;
    }
    if (ps1->PermID > ps2->PermID) {
      return 1;
    }
    else { // Same Perm ID number
      return 0;
    }
  }
  else { // Different complexes
    if (sumint(ps1->Aj,ps1->numSS) > sumint(ps2->Aj,ps2->numSS)) {
      return 1;
    }
    else if (sumint(ps1->Aj,ps2->numSS) < sumint(ps2->Aj,ps2->numSS)) {
      return -1;
    }
    else { // Have same number of strands
      for (i = 0; i < ps1->numSS; i++) {
	if (ps1->Aj[i] > ps2->Aj[i]) {
	  return -1;
	}
	else if (ps1->Aj[i] < ps2->Aj[i]) {
	  return 1;
	}
      }
      return 0; // Shouldn't ever get here
    }
  }

}
/* ******************************************************************************** */

/* ******************************************************************************** */
int Compare25(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Whichever has highest expected count (<mj>) returns -1.

     If the concentrations are the same (usually only the case if they're both zero),
     they are sorted by complex ID.
  */

  const struct CompSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct CompSortStruct *ps2 = p2;

  if (ps1->mj < ps2->mj) {
    return 1;
  }
  else if (ps1->mj > ps2->mj) {
    return -1;
  }
  else { // Equal complex concentration
    if (ps1->CompID < ps2->CompID) {
      return -1;
    }
    else if (ps1->CompID > ps2->CompID) {
      return 1;
    }
    else {  // Shouldn't ever get here
      return 0;
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare26(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Whichever has highest complex ID returns -1
  */

  const struct CompSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct CompSortStruct *ps2 = p2;

  if (ps1->CompID < ps2->CompID) {
    return -1;
  }
  else if (ps1->CompID > ps2->CompID) {
    return 1;
  }
  else { // Only get here if there's an error (already taken care of)
    return 0;
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int Compare27(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Used for sorting first by number of strands in complex, and then
     alphabetically.  I.e., an order of complexes might be:
     A, B, AA, AB, BB

     This is used when Method = 1.
  */

  int i; // Counter
  const struct CompSortStruct *ps1 = p1;  // Get the right type of pointer
  const struct CompSortStruct *ps2 = p2;

  if (sumint(ps1->Aj,ps1->numSS) > sumint(ps2->Aj,ps2->numSS)) {
    return 1;
  }
  else if (sumint(ps1->Aj,ps2->numSS) < sumint(ps2->Aj,ps2->numSS)) {
    return -1;
  }
  else { // Have same number of strands
    for (i = 0; i < ps1->numSS; i++) {
      if (ps1->Aj[i] > ps2->Aj[i]) {
	return -1;
      }
      else if (ps1->Aj[i] < ps2->Aj[i]) {
	return 1;
      }
    }
  }

  // Only get here if there's some sort of error, already taken care of in ReadInput
  return 0;

}
/* ******************************************************************************** */


