/*
  InputFileReader.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  For use with Distributions.c. 

  Reads the input from the .cx and .count files.  The input is sorted
  for use in CalcDist.

  After the comment lines, which must begin with a % character, each
  row of the cxFile contains the complex ID and its free energy in
  kcal/mol.  E.g., if permutation 2 of complex ID number 17 is ABBD
  and has a free energy -18.62 kcal/mol, and there are four monomer
  types, the corresponding row in the cxFile would be:
      17  2  1  2  0  1  -18.62

  If the input NoPermID == 1, i.e., the second column above doesn't
  exist, the entry in the file looks like:
      17  1  2  0  1  -18.62

  The countFile contains the initial counts of the single-strands
  (entered as integers).  The last line contains the volume of the
  solution (or "box") to be considered in the calculation.  The volume
  is entered in units of LITERS, and may be entered in scientific
  notation, e.g., 1.3e-21.

  For further info about input and output formatting, see associated
  manual.

  WARNING: This program does little format checking, so if there is an
  error in the input, it is likely to result in some strange error
  and/or seg fault.
*/

#include "InputFileReader.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"

// Structures for storing input and subsequent sorting
struct CompStruct { // Struct for complexes (used for output)
  int *Aj; // Array representing column j of A
  int numSS; // number of entries in Aj
  int CompID;
  double FreeEnergy; // Partition function for species
  char *AuxStr; // String containing auxillary information from input file
};

struct PermStruct { // Struct for complexes (used for output)
  int *Aj; // Array representing column j of A
  int numSS; // number of entries in Aj
  int CompID;
  int PermID;
  double FreeEnergy; // Partition function for species
  char *AuxStr; // String containing auxillary information from input file
};


/* ******************************************************************************** */
void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
             int *numPermsArray[], int *nComments, char *cxFile, char *countFile, 
             int quiet) {
 
 /*
   Finds the number of single-strands and the number of complexes
     in the input file using the information in cxFile and countFile.
     
     THE MEMORY FOR numPermsArray IS ALLOCATED IN THIS FUNCTION AND MUST
     BE FREED OUTSIDE OF IT.
     */
 
 int j; // Counter
 char line[MAXLINE]; // A junk buffer to dump the lines as we count them
 char *tok; // Token
 char tokseps[] = " \t\n"; // Token separators
 int CompID; // Complex ID number  
 FILE *fp; // The file we're reading from
 
 /* *************** Find the number of single species ********************* */
 // Open the count file
 if ((fp = fopen(countFile,"r")) == NULL) {
   if (quiet == 0) {
     printf("Error in opening file %s!\n",countFile);
     printf("\nExiting....\n\n");
   }
   exit(ERR_COUNT);
 }
 
 // Blow through comments and blank lines
 while (fgets(line,MAXLINE,fp) != NULL && 
        (line[0] == '%' || line[0] == '\0' || line[0] == '\n'));
 
 // Count the lines in the con file to see how many ss species there are
 *numSS = 1;  // We already read the first line of input
 while (fgets(line,MAXLINE,fp) != NULL && line[0] != '\0'
        && line[0] != '\n') {
          (*numSS)++;
        }
 fclose(fp);
 
 // last line of con file is volume, not another species count
 (*numSS)--;
 /* *********************************************************************** */
 
 
 /* *************** Find the maximum complex ID number. ******************* */
 // Open the cx file
 if ((fp = fopen(cxFile,"r")) == NULL) {
   if (quiet == 0) {
     printf("Error in opening file %s!\n",cxFile);
     printf("\nExiting....\n\n");
   }
   exit(ERR_CX);
 }
 
 // Blow through comments and blank lines
 *nComments = 0;
 while (fgets(line,MAXLINE,fp) != NULL && 
        (line[0] == '%' || line[0] == '\0' || line[0] == '\n')) {
          (*nComments)++;
        }
 
 *nTotal = 1; // Initialize count of total number of permutations counter
 *numTotal = 1;  // This is the total number of complexes (initialize to one --
 // assuming .cx file is nonempty).
 while (fgets(line,MAXLINE,fp) != NULL && line[0] != '\0'
        && line[0] != '\n') {
          if (line[0]=='%') continue;
          CompID = atoi(strtok(line,tokseps));
          if (CompID > *numTotal) {
            *numTotal = CompID;
          }
          (*nTotal)++; // Advance total number of permuations counter
        }
 fclose(fp);   
 *LargestCompID = *numTotal;
 /* *********************************************************************** */
 
 
 /* ******* Get the number of perumations for each complex **************** */
 // Allocate memory and initialize array of number of permutations
  (*numPermsArray) = (int *) malloc((*numTotal) * sizeof(int));
 for (j = 0; j < (*numTotal); j++) {
   (*numPermsArray)[j] = 0;
 }
 
 // Open the cx file
 if ((fp = fopen(cxFile,"r")) == NULL) {
   if (quiet == 0) {
     printf("Error in opening file %s!\n",cxFile);
     printf("\nExiting....\n\n");
   }
   exit(ERR_CX);
 }
 
 // Blow through comments and blank lines
 for (j = 0; j < *nComments; j++) {
   fgets(line,MAXLINE,fp);
 }

 while (fgets(line,MAXLINE,fp) != NULL && line[0] != '\0' && line[0] != '\n') {
   tok = strtok(line,tokseps);   // Complex ID
   if (tok[0]=='%') continue; // Skip comments
   CompID = atoi(tok) - 1;
   ((*numPermsArray)[CompID])++;
 }
 fclose(fp);
 
 // Check to make sure the complex ID's are sequential
 for (j = 0; j < *numTotal; j++) {
   if ((*numPermsArray)[j] == 0) {
     if (quiet == 0) {
       printf("Input cx file must contain all complex ID numbers between\n");
       printf("1 and the total number of complexes!\n\n");
       printf("Exiting....\n\n");
     }
     exit(ERR_NONSEQUENTIAL);
   }
 }
 /* *********************************************************************** */
 
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void ReadInputFiles(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
                    int **m0, double *M, int *numSS, int *numSS0, int *numTotal,
                    int *numPermsArray, char *cxFile, char *countFile, double *kT,
                    int Toverride, char *logFile, char *eqFile, int quiet, 
                    int WriteLogFile) {
  /*
    If one of the entries in the con file is zero, the problem is
      reformulated as if that strand does not exist.
      
      The input is stored in the arrays A, G, and m0, and also in M.
      A[i][j] is the number of monomers of type i in complex j.  G[j] is
        the free energy of complex j IN UNITS OF kT.  m 0[i] is the count
          of monomer species i.  M is the number of solvent particles in the
          box, as calculated from the supplied volume (in liters) and the
            density of water.
            
            The arrays CompIDArray and PermIDArray store the corresponding
            complex IDs and Permutation IDs for the entries loaded in the A
            and G.
            
            This function is used for calculations where only complex
            information is outputed (i.e., permutation information is not
                                     relevant).  ReadInputFilesPerm is used when permutation
                                       information is outputted.
                                       
                                       THE MEMORY FOR THESE ARRAYS IS ALLOCATED IN THIS FUNCTION AND MUST
                                       BE FREED OUTSIDE OF IT.
                                       */
  
  int i,j,k; // Counters
  struct CompStruct *InputStruct; // Struct we store the input in.
  char line[MAXLINE]; // A line from a file
  char *tok; // Token
  char *getResult;
  char tokseps[] = " \t\n"; // Token separators
  int nSS; // Local number of single species
  int cTotal; // Local number of complexes
  int ComplexID; // Complex ID for the line of imput file we are considering
  int *nonzerom0; // Identities of strands that are not zero in ccon
  int *zerom0; // The identities of strands that are set to zero in ccon
  int **newA; // The matrix A for the reformulated problem with zero ccon's taken out
  int *newCompIDArray; // The comp ID's for the reformulated problem
  int *newPermIDArray; // The perm ID's for the reformulated problem
  double *newG; // Free energies for reformulated prob. with zero ccon's taken out
  double *newm0; // Mole fractions of single species with zero ccon's taken out
  long double *Q=0; // Partition functions for complexes
  long double addQ; // A summand in the sum representing Q.
  double Gperm; // free energy of a given permutation
  double m0check; // Used to check input for counts 
  int newnumTotal; // New number of complexes after zero ccon's are removed
  int newnumSS; // New number of single strands after zero ccon's are removed
  int notOK; // Whether or not an entry in A can be kept if there are zero ccon's
  int noPerms; // noPerms = 1 if permutations are not explicitly considered
  int LineOK; // = 1 is the next line in the file is not NULL
  double MolesWaterPerLiter; // Moles of water per liter
  FILE *fp; // Handle for files we open
  FILE *fplog,*fpeq; // Handle for log and eq files
  
  
  // Rename these just so we don't have to use cumbersome pointers
  nSS = *numSS;
  cTotal = *numTotal;
  
  // Record the number of monomer types including those with zero conc.
  *numSS0 = nSS;
  
  // Find out if we need to explicitly consider permutations
  if (sumint(numPermsArray,cTotal) > cTotal) {
    noPerms = 0;
  }
  else {
    noPerms = 1;
  }
  
  // Allocate memory for A, G, and m0
  // THESE ARE NOT FREED UNTIL THE END OF MAIN
  *A = (int **) malloc(nSS * sizeof(int *));
  for (i = 0; i < nSS; i++) {
    (*A)[i] = (int *) malloc(cTotal * sizeof(int));
  }
  *G = (double *) malloc(cTotal * sizeof(double));
  *CompIDArray = (int *) malloc(cTotal * sizeof(int));
  *PermIDArray = (int *) malloc(cTotal * sizeof(int));
  // For this PermIDArray is all -1
  for (j = 0; j < cTotal; j++) {
    (*PermIDArray)[j] = -1;
  }
  *m0 = (int *) malloc(nSS * sizeof(int));
  
  // Allocate memory for the struct
  InputStruct = (struct CompStruct *) malloc(cTotal * sizeof(struct CompStruct));
  for (j = 0; j < cTotal; j++) {
    InputStruct[j].Aj = (int *) malloc (nSS * sizeof(int));
  }
  
  // Allocate memory for the partition functions and initialize
  if (noPerms == 0) {
    Q = (long double *) malloc (cTotal * sizeof(long double));
    for (j = 0; j < cTotal; j++) {
      Q[j] = 0.0;
    }
  }
  
  /* ************ Read in intial concentrations ********************* */
  // Open the count file
  if ((fp = fopen(countFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",countFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_COUNT);
  }
  
  // Blow through comments and blank lines
  while (fgets(line,MAXLINE,fp) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n'));
  
  // Read in the initial concentrations
  if ( (tok = strtok(line,tokseps)) != NULL) {
    // Do a quick check to make sure user isn't using molar concentrations
    m0check = str2double(tok);
    if (m0check < 1) {
      if (quiet == 0) {
        printf("Error in con file: Initial counts must be integers.\n\n");
        printf("Exiting....\n\n");
      }
      exit(ERR_NONINTEGER);
    }
    (*m0)[0] = atoi(tok);
  }
  for (i = 1; i < nSS; i++) {
    while (fgets(line,MAXLINE,fp) != NULL && 
           (line[0] == '%' || line[0] == '\0' || line[0] == '\n'));
    
    if ( (tok = strtok(line,tokseps)) != NULL) {
      // Do a quick check to make sure user isn't using molar concentrations
      m0check = str2double(tok);
      if (m0check < 1) {
        if (quiet == 0) {
          printf("Error in con file: Initial counts must be integers.\n\n");
          printf("Exiting....\n\n");
        }
        exit(ERR_NONINTEGER);
      }
      (*m0)[i] = atoi(tok);
    }
  } 
  // Read in the volume (number of solvent particles)
  fgets(line,MAXLINE,fp);
  if ( (tok = strtok(line,tokseps)) != NULL) {
    (*M) = AVOGADRO * str2double(tok);
  }
  fclose(fp);
  /* *************************************************************** */
  
  /* *********** Write information to log and dist files *********** */
  // First write initial concentrations
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    for (i = 0; i < nSS; i++) {
      fprintf(fplog,"       %d: %d\n",i+1,(*m0)[i]);
    }
    fprintf(fplog,"   --Box volume (liters): ");
    fprintf(fplog,"%8.6e\n",(*M)/AVOGADRO);
    fclose(fplog);
  }
  
  if ((fpeq = fopen(eqFile,"a")) == NULL) {
    if (quiet == 0) {
      printf("Error opening %s.\n\nExiting....\n",eqFile);
    }
    exit(ERR_DIST);
  }
  for (i = 0; i < nSS; i++) {
    fprintf(fpeq,"%%   %d: %d\n",i+1,(*m0)[i]);
  }
  fprintf(fpeq,"%% Box volume (liters): ");
  fprintf(fpeq,"%8.6e\n",(*M)/AVOGADRO);
  if (Toverride == 1) {
    fprintf(fpeq,"%% User supplied temperature of %g\n",(*kT)/kB - ZERO_C_IN_KELVIN);
  }
  fprintf(fpeq,"\n"); // Extra blank line to separate comments from other code
  /* *************************************************************** */
  
  /* ************** Read in A, free energy, and complex IDs ******** */
  // Open the cx file
  if ((fp = fopen(cxFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",cxFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_CX);
  }
  
  // Blow through comments and blank lines and pull out the new kT if necessary
  while (fgets(line,MAXLINE,fp) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n')) {
           
           // Print comment line to output file
           fprintf(fpeq,"%s",line); 
           
           if (Toverride == 0) {
             if (line[0] == '%' && line[1] == ' ' && line[2] == 'T' && line[3] == ' '
                 && line[4] == '=' && line[5] == ' ') { // This is line with temperature data
                   tok = strtok(line,tokseps); // tok = '%'
                   tok = strtok(NULL,tokseps); // tok = 'T'
                   tok = strtok(NULL,tokseps); // tok = '='
                   tok = strtok(NULL,tokseps); // This is the temperature
                   *kT = kB*(str2double(tok) +  ZERO_C_IN_KELVIN);
                 }
           }
         }
  
  // Close eqFile
  fclose(fpeq);
  
  // Build A and Free Energy.
  LineOK = 1;
  while (LineOK == 1 && line[0] != '\0' && line[0] != '\n') {
    // Get the complex ID
    tok = strtok(line,tokseps);  // Complex ID
    ComplexID = atoi(tok) - 1;
    InputStruct[ComplexID].CompID = ComplexID + 1;
    
    // Pull out column of A corresponding to complex (this is done redundantly)
    for (i = 0; i < nSS; i++) {
      if ((tok = strtok(NULL,tokseps)) != NULL) {
        InputStruct[ComplexID].Aj[i] = atoi(tok);
      }
      else {
        if (quiet == 0) {
          printf("Error in input file!\n\nExiting....\n");
        }
        exit(ERR_BADROWINP);
      }
    }
    
    // Enter the free energy
    if ((tok = strtok(NULL,tokseps)) != NULL) {
      Gperm = str2double(tok)/(*kT);
      if (noPerms) {
        InputStruct[ComplexID].FreeEnergy = Gperm;
      }
      else {
        if ((addQ = expl(-((long double) Gperm))) == HUGE_VAL){
          if (quiet == 0) {
            printf("Free energies of structures are too high to allow calculation\n");
            printf("including permuations.  Reformat the problem such that each\n");
            printf("complex has its own free energy (no permuations).\n\n");
            printf("Exiting.....\n");
          }
          exit(ERR_NOPERMS);
        }
        else {
          Q[ComplexID] += addQ;
        }
      }
      
    }
    else {
      if (quiet == 0) {
        printf("Error in input file!\n\nExiting....\n");
      }
      exit(ERR_BADROWINP);
    }
    
    // Put numSS in just because we have to for qsort
    InputStruct[ComplexID].numSS = nSS;
    
    // Read in the next line
    
    do {
      getResult=fgets(line,MAXLINE,fp);
      if (!getResult) {
        LineOK = 0;
        break;
      }
    } while (line[0] == '%' || line[0] == '\0' || line[0] == '\n');
  }
  fclose(fp);
  
  // Compute and enter free energies
  if (noPerms == 0) {
    for (j = 0; j < cTotal; j++) {
      InputStruct[j].FreeEnergy = -logl(Q[j]);
    }
    free(Q);
  }
  /* *************************************************************** */
  
  
  // Sort structure (sorting puts single-stranded entries first)
  qsort(InputStruct,cTotal,sizeof(struct CompStruct),InputCompare);
  
  
  // Make the matrix A and free energy G and the complex ID list
  for (j = 0; j < cTotal; j++) {
    for (i = 0; i < nSS; i++) {
      (*A)[i][j] = InputStruct[j].Aj[i];
    }
    (*G)[j] = InputStruct[j].FreeEnergy;
    (*CompIDArray)[j] = InputStruct[j].CompID;
  }
  
  // Free the struct
  for (j = 0; j < cTotal; j++) {
    free(InputStruct[j].Aj);
  }
  free(InputStruct);
  
  
  /* ************** BEGIN    REFORMATTING PROBLEM *********************** */
  /* 
    This section of the code reformats the problem if there are zero
      entries in the con file.  I.e., if there is a variable whose
      initial concentration is entered as zero, the problem is
      reformulated as if that strand doesn't exist and if a complex
      cannot be formed, the problem is reformulated as if it doesn't
      exist.  Arrays are reallocated to the adjusted size of the
      problem
      */
  
  zerom0 = (int *) malloc(nSS * sizeof(int));
  nonzerom0 = (int *) malloc(nSS * sizeof(int));
  
  // Check to see if any of the concentrations are zero.  If so, 
  // reformulate the problem accordingly.
  j = 0;  // How many entries are zero
  k = 0;  // How many entries are nonzero
  for (i = 0; i < nSS; i++) {
    if ((*m0)[i] == 0) {
      zerom0[j] = i;
      j++;
    }
    else {
      nonzerom0[k] = i;
      k++;
    }
  }
  
  newnumSS = nSS - j;
  
  // First count how many entries we have.  notOK = 1 if complex cannot be made
  // with initial counts
  newnumTotal = 0;
  
  for (j = 0; j < cTotal; j++) {
    notOK = 0;
    for (i = 0; i < nSS; i++) {
      if ((*A)[i][j] > (*m0)[i]) {
        notOK = 1;
      }
    }
    if (notOK == 0) {
      newnumTotal++;
    }
  }
  
  // Allocate memory for new arrays
  newA = (int **) malloc(newnumSS * sizeof(int *));
  for (i = 0; i < newnumSS; i++) {
    newA[i] = (int *) malloc(newnumTotal * sizeof(int));
  }
  newG =  (double *) malloc(newnumTotal * sizeof(double));
  newm0 = (double *) malloc(newnumSS * sizeof(double));
  newCompIDArray = (int *) malloc(newnumTotal * sizeof(int));
  newPermIDArray = (int *) malloc(newnumTotal * sizeof(int));
  
  // Put in the new m0
  for (i = 0; i < newnumSS; i++) {
    newm0[i] = (*m0)[nonzerom0[i]];
  }
  
  // Go through and pick out the entries to keep
  k = 0;
  for (j = 0; j < cTotal; j++) {
    notOK = 0;
    for (i = 0; i < nSS; i++) {
      if ((*A)[i][j] > (*m0)[i]) {
        notOK = 1;
      }
    }
    if (notOK == 0) {
      for (i = 0; i < newnumSS; i++) {
        newA[i][k] = (*A)[nonzerom0[i]][j];
      }
      newG[k] = (*G)[j];
      newCompIDArray[k] = (*CompIDArray)[j];
      newPermIDArray[k] = (*PermIDArray)[j];
      k++;
    }
  }
  
  // Change names of "new" variables
  // Rename newA
  for (i = 0; i < nSS; i++) {
    free((*A)[i]);
  }
  free(*A);    
  *A = (int **) malloc(newnumSS * sizeof(int *));
  for (i = 0; i < newnumSS; i++) {
    (*A)[i] = (int *) malloc(newnumTotal * sizeof(int *));
    for (j = 0; j < newnumTotal; j++) {
      (*A)[i][j] = newA[i][j];
    }
    free(newA[i]);
  }
  free(newA);
  
  // Rename newG
  free(*G);
  *G = (double *) malloc(newnumTotal * sizeof(double));
  for (j = 0; j < newnumTotal; j++) {
    (*G)[j] = newG[j];
  }
  free(newG);
  
  // Rename CompIDArray
  free(*CompIDArray);
  free(*PermIDArray);
  *CompIDArray = (int *) malloc(newnumTotal * sizeof(int));
  *PermIDArray = (int *) malloc(newnumTotal * sizeof(int));
  for (j = 0; j < newnumTotal; j++) {
    (*CompIDArray)[j] = newCompIDArray[j];
    (*PermIDArray)[j] = newPermIDArray[j];
  }
  free(newCompIDArray);
  free(newPermIDArray);
  
  // Rename newm0
  free(*m0);
  (*m0) = (int *) malloc(newnumSS * sizeof(int));
  for (i = 0; i < newnumSS; i++) {
    (*m0)[i] = (int) (newm0[i]+0.1); // Cast it as an int
  }
  free(newm0);
  
  // Rename numTotal and numSS
  *numTotal = newnumTotal;
  *numSS = newnumSS;
  
  free(zerom0);
  free(nonzerom0);
  /* ************** FINISHED REFORMATTING PROBLEM *********************** */
  
  // Calculate molarity of water and convert appropriate quantities to the right units
  MolesWaterPerLiter = WaterDensity((*kT)/kB - ZERO_C_IN_KELVIN);
  (*M) *= MolesWaterPerLiter;
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
void ReadInputFilesPerm(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
                        int **m0, double *M, int *numSS, int *numSS0, 
                        int *newnTotal, int nTotal, char *cxFile, char *countFile, 
                        double *kT, int Toverride, char *logFile, char *eqFile, 
                        int quiet, int WriteLogFile) {
  /*
    This input file reader is the same as above, but used when
      permutation information is to be outputted as well.
      
      THE MEMORY FOR THESE ARRAYS IS ALLOCATED IN THIS FUNCTION AND MUST
      BE FREED OUTSIDE OF IT.
      */
  
  int i,j,k; // Counters
  struct PermStruct *InputStruct; // Struct we store the input in.
  char line[MAXLINE]; // A line from a file
  char *tok; // Token
  char *getResult;
  char tokseps[] = " \t\n"; // Token separators
  int nSS; // Local number of single species
  int *nonzerom0; // Identities of strands that are not zero in ccon
  int *zerom0; // The identities of strands that are set to zero in ccon
  int **newA; // The matrix A for the reformulated problem with zero ccon's taken out
  int *newCompIDArray; // The comp ID's for the reformulated problem
  int *newPermIDArray; // The perm ID's for the reformulated problem
  double *newG; // Free energies for reformulated prob. with zero ccon's taken out
  int *newm0; // Mole fractions of single species with zero ccon's taken out
  double m0check; // Used to check input for counts
  int newnumTotal; // New number of complexes after zero ccon's are removed
  int newnumSS; // New number of single strands after zero ccon's are removed
  int LineOK; // = 1 is the next line in the file is not NULL
  int notOK; // Whether or not an entry in A can be kept if there are zero ccon's
  double MolesWaterPerLiter; // Moles of water per liter
  FILE *fp; // Handle for the cx file
  FILE *fplog,*fpeq; // Handle for log and eq files
  
  // Rename this just so we don't have to use cumbersome pointers
  nSS = *numSS;
  
  // Record the number of monomer types including those with zero conc.
  *numSS0 = nSS;
  
  // Allocate memory for A, G, and m0
  // THESE ARE NOT FREED UNTIL THE END OF MAIN
  *A = (int **) malloc(nSS * sizeof(int *));
  for (i = 0; i < nSS; i++) {
    (*A)[i] = (int *) malloc(nTotal * sizeof(int));
  }
  *G = (double *) malloc(nTotal * sizeof(double));
  *CompIDArray = (int *) malloc(nTotal * sizeof(int));
  *PermIDArray = (int *) malloc(nTotal * sizeof(int));
  *m0 = (int *) malloc(nSS * sizeof(int));
  
  // Allocate memory for the struct
  InputStruct = (struct PermStruct *) malloc(nTotal * sizeof(struct PermStruct));
  for (j = 0; j < nTotal; j++) {
    InputStruct[j].Aj = (int *) malloc (nSS * sizeof(int));
  }
  
  /* ************ Read in intial counts **************************** */
  // Open the count file
  if ((fp = fopen(countFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",countFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_COUNT);
  }
  
  // Blow through comments and blank lines
  while (fgets(line,MAXLINE,fp) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n'));
  
  if ( (tok = strtok(line,tokseps)) != NULL) {
    // Do a quick check to make sure user isn't using molar concentrations
    m0check = str2double(tok);
    if (m0check < 1) {
      if (quiet == 0) {
        printf("Error in con file: Initial counts must be integers.\n\n");
        printf("Exiting....\n\n");
      }
      exit(ERR_NONINTEGER);
    }
    (*m0)[0] = atoi(tok);
  }
  for (i = 1; i < nSS; i++) {
    fgets(line,MAXLINE,fp);
    if ( (tok = strtok(line,tokseps)) != NULL) {
      // Do a quick check to make sure user isn't using molar concentrations
      m0check = str2double(tok);
      if (m0check < 1) {
        if (quiet == 0) {
          printf("Error in con file: Initial counts must be integers.\n\n");
          printf("Exiting....\n\n");
        }
        exit(ERR_NONINTEGER);
      }
      (*m0)[i] = atoi(tok);
    }
  } 
  // Read in the volume (number of solvent particles)
  fgets(line,MAXLINE,fp);
  if ( (tok = strtok(line,tokseps)) != NULL) {
    (*M) = AVOGADRO * str2double(tok);
  }
  fclose(fp);
  /* *************************************************************** */
  
  /* *********** Write information to log and eq files ************* */
  // First write initial concentrations
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    for (i = 0; i < nSS; i++) {
      fprintf(fplog,"       %d: %d\n",i+1,(*m0)[i]);
    }
    fprintf(fplog,"   --Box volume (liters): ");
    fprintf(fplog,"%8.6e\n",(*M)/AVOGADRO);
    fclose(fplog);
  }
  
  if ((fpeq = fopen(eqFile,"a")) == NULL) {
    if (quiet == 0) {
      printf("Error opening %s.\n\nExiting....\n",eqFile);
    }
    exit(ERR_DIST);
  }
  for (i = 0; i < nSS; i++) {
    fprintf(fpeq,"%%   %d: %d\n",i+1,(*m0)[i]);
  }
  fprintf(fpeq,"%% Box volume (liters): ");
  fprintf(fpeq,"%8.6e\n",(*M)/AVOGADRO);
  if (Toverride == 1) {
    fprintf(fpeq,"%% User supplied temperature of %g\n",(*kT)/kB - ZERO_C_IN_KELVIN);
  }
  fprintf(fpeq,"\n"); // Extra blank line to separate comments from other code
  /* *************************************************************** */
  
  /* ************** Read in A, free energy, and complex IDs ******** */
  // Open the cx file
  if ((fp = fopen(cxFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",cxFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_CX);
  }
  
  // Blow through comments and blank lines and pull out the new kT if necessary
  while (fgets(line,MAXLINE,fp) != NULL && 
         (line[0] == '%' || line[0] == '\0' || line[0] == '\n')) {
           if (Toverride == 0) {
             if (line[0] == '%' && line[1] == ' ' && line[2] == 'T' && line[3] == ' '
                 && line[4] == '=' && line[5] == ' ') { // This is line with temperature data
                   tok = strtok(line,tokseps); // tok = '%'
                   tok = strtok(NULL,tokseps); // tok = 'T'
                   tok = strtok(NULL,tokseps); // tok = '='
                   tok = strtok(NULL,tokseps); // This is the temperature
                   *kT = kB*(str2double(tok) +  ZERO_C_IN_KELVIN);
                 }
           }
         }
  
  // Close eqFile
  fclose(fpeq);
  
  // Build A and Free Energy.
  LineOK = 1;
  k = 0; // Counter for entries in input file
  while (LineOK == 1 && line[0] != '\0' && line[0] != '\n') {
    // Get the complex ID
    tok = strtok(line,tokseps);  // Complex ID
    InputStruct[k].CompID = atoi(tok);
    
    // Permutation number
    tok = strtok(NULL,tokseps);
    InputStruct[k].PermID = atoi(tok);
    
    // Pull out column of A corresponding to complex (this is done redundantly)
    for (i = 0; i < nSS; i++) {
      if ((tok = strtok(NULL,tokseps)) != NULL) {
        InputStruct[k].Aj[i] = atoi(tok);
      }
      else {
        if (quiet == 0) {
          printf("Error in input file!\n\nExiting....\n");
        }
        exit(ERR_BADROWINP);
      }
    }
    
    // Enter the free energy
    if ((tok = strtok(NULL,tokseps)) != NULL) {
      InputStruct[k].FreeEnergy = str2double(tok)/(*kT);
    }
    else {
      if (quiet == 0) {
        printf("Error in input file!\n\nExiting....\n");
      }
      exit(ERR_BADROWINP);
    }
    
    // Put numSS in just because we have to for qsort
    InputStruct[k].numSS = nSS;
    
    do {
      getResult=fgets(line,MAXLINE,fp);
      if (getResult && !quiet) {
#ifdef DEBUG
        printf("Got line     : %s",line);
#endif
      }
      if (!getResult) {
        LineOK = 0;
#ifdef DEBUG
          printf("LineOK=0\n");
#endif
        break;
      }
    } while (line[0] == '%' || line[0] == '\0' || line[0] == '\n');
    if (getResult) {
#ifdef DEBUG
      printf("ACCEPTED line: %s",line);
#endif
    }
    k++;
  }
  fclose(fp);
  /* *************************************************************** */
  
  // Sort structure (sorting puts single-stranded entries first)
  qsort(InputStruct,nTotal,sizeof(struct PermStruct),InputComparePerm);
  
  // Make the matrix A and free energy G and the complex ID list
  for (j = 0; j < nTotal; j++) {
    for (i = 0; i < nSS; i++) {
      (*A)[i][j] = InputStruct[j].Aj[i];
    }
    (*G)[j] = InputStruct[j].FreeEnergy;
    (*CompIDArray)[j] = InputStruct[j].CompID;
    (*PermIDArray)[j] = InputStruct[j].PermID;
  }
  
  // Free the struct
  for (j = 0; j < nTotal; j++) {
    free(InputStruct[j].Aj);
  }
  free(InputStruct);
  
  
  /* ************** BEGIN    REFORMATTING PROBLEM *********************** */
  /* 
    This section of the code reformats the problem if there are zero
      entries in the con file.  I.e., if there is a variable whose
      initial concentration is entered as zero, the problem is
      reformulated as if that strand doesn't exist and if a complex
      cannot be formed, the problem is reformulated as if it doesn't
      exist.  Arrays are reallocated to the adjusted size of the
      problem
      */
  
  zerom0 = (int *) malloc(nSS * sizeof(int));
  nonzerom0 = (int *) malloc(nSS * sizeof(int));
  
  // Check to see if any of the concentrations are zero.  If so, 
  // reformulate the problem accordingly.
  j = 0;  // How many entries are zero
  k = 0;  // How many entries are nonzero
  for (i = 0; i < nSS; i++) {
    if ((*m0)[i] == 0) {
      zerom0[j] = i;
      j++;
    }
    else {
      nonzerom0[k] = i;
      k++;
    }
  }
  
  newnumSS = nSS - j;
  // First count how many entries we have.  notOK = 1 if complex cannot be made
  // with initial counts
  newnumTotal = 0;
  
  for (j = 0; j < nTotal; j++) {
    notOK = 0;
    for (i = 0; i < nSS; i++) {
      if ((*A)[i][j] > (*m0)[i]) {
        notOK = 1;
      }
    }
    if (notOK == 0) {
      newnumTotal++;
    }
  }
  
  // Allocate memory for new arrays
  newA = (int **) malloc(newnumSS * sizeof(int *));
  for (i = 0; i < newnumSS; i++) {
    newA[i] = (int *) malloc(newnumTotal * sizeof(int));
  }
  newG =  (double *) malloc(newnumTotal * sizeof(double));
  newm0 = (int *) malloc(newnumSS * sizeof(int));
  newCompIDArray = (int *) malloc(newnumTotal * sizeof(int));
  newPermIDArray = (int *) malloc(newnumTotal * sizeof(int));
  
  // Put in the new m0
  for (i = 0; i < newnumSS; i++) {
    newm0[i] = (*m0)[nonzerom0[i]];
  }
  
  // Go through and pick out the entries to keep
  k = 0;
  for (j = 0; j < nTotal; j++) {
    notOK = 0;
    for (i = 0; i < nSS; i++) {
      if ((*A)[i][j] > (*m0)[i]) {
        notOK = 1;
      }
    }
    if (notOK == 0) {
      for (i = 0; i < newnumSS; i++) {
        newA[i][k] = (*A)[nonzerom0[i]][j];
      }
      newG[k] = (*G)[j];
      newCompIDArray[k] = (*CompIDArray)[j];
      newPermIDArray[k] = (*PermIDArray)[j];
      k++;
    }
  }
  
  // Change names of "new" variables
  // Rename newA
  for (i = 0; i < nSS; i++) {
    free((*A)[i]);
  }
  free(*A);    
  *A = (int **) malloc(newnumSS * sizeof(int *));
  for (i = 0; i < newnumSS; i++) {
    (*A)[i] = (int *) malloc(newnumTotal * sizeof(int *));
    for (j = 0; j < newnumTotal; j++) {
      (*A)[i][j] = newA[i][j];
    }
    free(newA[i]);
  }
  free(newA);
  
  // Rename newG
  free(*G);
  *G = (double *) malloc(newnumTotal * sizeof(double));
  for (j = 0; j < newnumTotal; j++) {
    (*G)[j] = newG[j];
  }
  free(newG);
  
  // Rename CompIDArray and PermIDArray
  free(*CompIDArray);
  free(*PermIDArray);
  *CompIDArray = (int *) malloc(newnumTotal * sizeof(int));
  *PermIDArray = (int *) malloc(newnumTotal * sizeof(int));
  for (j = 0; j < newnumTotal; j++) {
    (*CompIDArray)[j] = newCompIDArray[j];
    (*PermIDArray)[j] = newPermIDArray[j];
  }
  free(newCompIDArray);
  free(newPermIDArray);
  
  // Rename newm0
  free(*m0);
  (*m0) = (int *) malloc(newnumSS * sizeof(int));
  for (i = 0; i < newnumSS; i++) {
    (*m0)[i] = newm0[i];
  }
  free(newm0);
  
  // Rename numTotal and numSS
  *newnTotal = newnumTotal;
  *numSS = newnumSS;
  
  free(zerom0);
  free(nonzerom0);
  
  /* ************** FINISHED REFORMATTING PROBLEM *********************** */
  
  // Calculate molarity of water and convert appropriate quantities to the right units
  MolesWaterPerLiter = WaterDensity((*kT)/kB - ZERO_C_IN_KELVIN);
  (*M) *= MolesWaterPerLiter;
  
}
/* ******************************************************************************** */


/* ******************************************************************************** */
int InputCompare(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.

     Used to sort input first by complex size and then
     "alphabetically" within each complex size.
  */

  int i; // Counter
  const struct CompStruct *ps1 = p1;  // Get the right type of pointer
  const struct CompStruct *ps2 = p2;

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

  // We only get here if there's a duplicate row in A (not mistake in input, but bug)
  // Should never get here.
  printf("Input error in complex ID listing.\n");
  exit(ERR_DUPROWA);

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int InputComparePerm(const void *p1, const void *p2) {
  /* 
     Comparison function (in mandatory form) to send to qsort.
     See Prata, C Primer Plus, 4th Ed. p. 654 for description.
  */

  int i; // Counter
  const struct PermStruct *ps1 = p1;  // Get the right type of pointer
  const struct PermStruct *ps2 = p2;

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

  // Duplicate listing
  return 0;

}
/* ******************************************************************************** */


