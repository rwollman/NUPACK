/*
  InputFileReader.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  For use with Concentrations.c. 

  Reads the input from the .cx and .con files.

  After the comment lines, which must begin with a % character, each
  row of the cxFile contains the complex ID and its free energy in
  kcal/mol.  E.g., if permutation 2 of complex ID number 17 is ABBD
  and has a free energy -18.64 kcal/mol, and there are four monomer
  types, the corresponding row in the cxFile would be:
      17  2  1  2  0  1  -18.62

  If the input NoPermID == 1, i.e., the second column above doesn't
  exist, the entry in the file looks like:
      17  1  2  0  1  -18.62

  The conFile has the initial concentrations of each of the monomer
  types in solution.  These are in units of MOLAR.

  For further info about input and output formatting, see associated
  manual.

  WARNING: This program does little format checking, so if there is an
  error in the input, it is likely to result in some strange error
  and/or seg fault.
*/


#include "InputFileReader.h" // File with important definitions
#include <float.h>
#include "constants.h"


// Structures for storing input and subsequent sorting
struct CompStruct { // Struct for complexes (used for output)
  int *Aj; // Array representing column j of A
  int numSS; // number of entries in Aj
  int CompID;
  double FreeEnergy; // Partition function for species
  char *AuxStr; // String containing auxillary information from input file
};


/* ******************************************************************************** */
void getSize(int *numSS, int *numTotal, int *nTotal, int *LargestCompID,
             int **numPermsArray, char *cxFile, char *conFile, int quiet) {

  /*
    Finds the number of single-strands and the number of complexes
    in the input file using the information in cxFile and conFile.

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
  // Open the con file
  if ((fp = fopen(conFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",conFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_CON);
  }

  // Count the lines in the con file to see how many ss species there are
  *numSS = 0;  // We already read the first line of input
  while (fgets(line,MAXLINE,fp) != NULL) {
    if (line[0] != '\0' && line[0] != '\n' && line[0] != '%') {
      (*numSS)++;
    }
  }
  fclose(fp);
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

  *nTotal = 0; // Initialize count of total number of permutations counter
  *numTotal = 1;  // This is the total number of complexes (initialize to one --
                  // assuming .cx file is nonempty).
  while (fgets(line,MAXLINE,fp) != NULL) {
    // skip commented out lines and blank lines
    if (line[0] != '%' && line[0] != '\n' && line[0] != '\0') {
      CompID = atoi(strtok(line,tokseps));
      if (CompID > *numTotal) {
        *numTotal = CompID;
      }
      (*nTotal)++; // Advance total number of permuations counter
    }
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

  while (fgets(line,MAXLINE,fp) != NULL) {
    if (line[0] != '\0' && line[0] != '\n' && line[0] != '%') {
      tok = strtok(line,tokseps);   // Complex ID
      CompID = atoi(tok) - 1;
      ((*numPermsArray)[CompID])++;
    }
  }
  fclose(fp);

  // Check to make sure the complex ID's are sequential
  for (j = 0; j < *numTotal; j++) {
    if ((*numPermsArray)[j] == 0) {
      if (quiet == 0) {
        printf("Input file must contain all complex ID numbers between\n");
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
double ReadInputFiles(int ***A, double **G, int **CompIDArray, int **PermIDArray, 
                      double **x0, int *numSS, int *numSS0, int *numTotal, 
                      int *numPermsArray, char *cxFile, char *conFile, double *kT, 
                      int Toverride, char  *logFile, char  *eqFile, 
                      char *fpairsFile, int quiet, int WriteLogFile, int DoBPfracs,
                      int NoPermID) {
  /*
    If one of the entries in the con file is zero, the problem is
    reformulated as if that strand does not exist.

    The input is stored in the arrays A, G, and x0.  A[i][j] is the
    number of monomers of type i in complex j.  G[j] is the free
    energy of complex j IN UNITS OF kT.  x0[i] is the initial mole
    fraction of monomer species i.

    The arrays CompIDArray and PermIDArray store the corresponding
    complex IDs and Permutation IDs for the entries loaded in the A
    and G.

    THE MEMORY FOR THESE ARRAYS IS ALLOCATED IN THIS FUNCTION AND MUST
    BE FREED OUTSIDE OF IT.
  */

  int i,j,k; // Counters
  struct CompStruct *InputStruct; // Struct we store the input in.
  char line[MAXLINE]; // A line from a file
  char *tok; // Token
  char tokseps[] = " \t\n"; // Token separators
  int nSS; // Local number of single species
  int cTotal; // Local number of complexes
  int ComplexID; // Complex ID for the line of imput file we are considering
  int *nonzerox0; // Identities of strands that are not zero in ccon
  int *zerox0; // The identities of strands that are set to zero in ccon
  int **newA; // The matrix A for the reformulated problem with zero ccon's taken out
  int *newCompIDArray; // The comp ID's for the reformulated problem
  int *newPermIDArray; // The perm ID's for the reformulated problem
  double *newG; // Free energies for reformulated prob. with zero ccon's taken out
  double *newx0; // Mole fractions of single species with zero ccon's taken out
  long double *Q; // Partition functions for complexes
  long double addQ; // A summand in the sum representing Q.
  double Gperm; // free energy of a given permutation
  int newnumTotal; // New number of complexes after zero ccon's are removed
  int newnumSS; // New number of single strands after zero ccon's are removed
  int notOK; // Whether or not an entry in A can be kept if there are zero ccon's
  int noPerms; // noPerms = 1 if permutations are not explicitly considered
  int LineOK; // = 1 is the next line in the file is not NULL
  double MolesWaterPerLiter; // Moles of water per liter
  FILE *fp; // Handle for files we open
  FILE *fpfpairs=0, *fplog=0, *fpeq=0; // file handles for fpairs, log and eq files

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

  // Allocate memory for A, G, and x0
  // THESE ARE NOT FREED UNTIL THE END OF MAIN
  *A = (int **) malloc(nSS * sizeof(int *));
  for (i = 0; i < nSS; i++) {
    (*A)[i] = (int *) malloc(cTotal * sizeof(int));
  }
  *G = (double *) malloc(cTotal * sizeof(double));
  *CompIDArray = (int *) malloc(cTotal * sizeof(int));
  *PermIDArray = (int *) malloc(cTotal * sizeof(int));
  // For this PermIDArray is all zeros
  for (j = 0; j < cTotal; j++) {
    (*PermIDArray)[j] = 0;
  }
  
  *x0 = (double *) malloc(nSS * sizeof(double));

  // Allocate memory for the struct
  InputStruct = (struct CompStruct *) malloc(cTotal * sizeof(struct CompStruct));
  for (j = 0; j < cTotal; j++) {
    InputStruct[j].Aj = (int *) malloc (nSS * sizeof(int));
  }
  
  // Allocate memory for the partition functions and initialize
  // We do this even if noPerms == 1 so the compiler doesn't give a warning when
  // optimization if turned on.
  Q = (long double *) malloc (cTotal * sizeof(long double));
  for (j = 0; j < cTotal; j++) {
    Q[j] = 0.0;
  }

  /* ************ Read in intial concentrations ********************* */
  // Open the con file
  if ((fp = fopen(conFile,"r")) == NULL) {
    if (quiet == 0) {
      printf("Error in opening file %s!\n",conFile);
      printf("\nExiting....\n\n");
    }
    exit(ERR_CON);
  }

  i = 0;
  while (fgets(line,MAXLINE,fp) != NULL) {
    if (line[0] != '%' && line[0] != '\0' && line[0] != '\n') {
      // Read in the initial concentration
      tok = strtok(line,tokseps);
      (*x0)[i] = str2double(tok);
      i++;
    }
  } 

  fclose(fp);
  /* *************************************************************** */


  /* *********** Write information to log and eq files ************* */
  // Write to log file
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    for (i = 0; i < nSS; i++) {
      fprintf(fplog,"       %d: %8.6e Molar\n",i+1,(*x0)[i]);
    }
    fprintf(fplog,"%%\n");
    fprintf(fplog,"%% Following is the header from the input file (%s):\n",cxFile);
    fprintf(fplog,"%%\n"); // Extra blank comment line to separate comments
  }

  // Write to eq file
  if ((fpeq = fopen(eqFile,"a")) == NULL) {
    if (quiet == 0) {
      printf("Error opening %s.\n\nExiting....\n",eqFile);
    }
    exit(ERR_EQ);
  }
  for (i = 0; i < nSS; i++) {
    fprintf(fpeq,"%%   %d: %8.6e Molar\n",i+1,(*x0)[i]);
  }

  if (Toverride == 1) {
    fprintf(fpeq,"%% User supplied temperature of %g\n",(*kT)/kB - ZERO_C_IN_KELVIN);
  }
  fprintf(fpeq,"%%\n");
  fprintf(fpeq,"%% Following is the header from the input file (%s):\n%%\n",cxFile);
  fprintf(fpeq,"%%\n"); // Extra blank line to separate comments

  // Write to fpairs file
  if (DoBPfracs) {
    if ((fpfpairs = fopen(fpairsFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",fpairsFile);
      }
      exit(ERR_FPAIRS);
    }
    for (i = 0; i < nSS; i++) {
      fprintf(fpfpairs,"%%   %d: %8.6e Molar\n",i+1,(*x0)[i]);
    }
    fprintf(fpfpairs,"%%\n");
    fprintf(fpfpairs,"%% Following is the header from the input file (%s):\n",cxFile);
    fprintf(fpfpairs,"%%\n"); // Extra blank comment line to separate comments
  }
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

    // Print comment line to output files
    if (WriteLogFile) {
      fprintf(fplog,"%s",line); 
    }
    fprintf(fpeq,"%s",line); 
    if (DoBPfracs) {
      fprintf(fpfpairs,"%s",line); 
    }

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

  // Close output files
  if (WriteLogFile) {
    fclose(fplog);
  }
  if (DoBPfracs) {
    fclose(fpfpairs);
  }

  // Build A and Free Energy.
  LineOK = 1;
  while (LineOK == 1) {
    if (line[0] == '%') { // If it's a comment, print it to output file
      fprintf(fpeq,"%s",line); 
    }
    else if (line[0] != '%' && line[0] != '\0' && line[0] != '\n') {
      // Get the complex ID
      tok = strtok(line,tokseps);  // Complex ID
      ComplexID = atoi(tok) - 1;
      InputStruct[ComplexID].CompID = ComplexID + 1;
      
      // Permutation number
      if (NoPermID == 0) {
        tok = strtok(NULL,tokseps);
      }
      
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
          if ((addQ = expl(-((long double) Gperm))) >= HUGE_VAL){
            if (quiet == 0) {
              printf("Free energies of complexes are too high for calculation\n");
              printf("including permutations.  Reformat the problem such that each\n");
              printf("complex has its own free energy (no permutations).\n\n");
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
    }  
    // Read in the next line
    if (fgets(line,MAXLINE,fp) == NULL) {
      LineOK = 0;
    }
    
  }
  fclose(fp);

  // Close eq file
  fclose(fpeq);

  // Compute and enter free energies
  if (noPerms == 0) {
    for (j = 0; j < cTotal; j++) {
      InputStruct[j].FreeEnergy = -(double) logl(Q[j]);
    }
  }
  /* *************************************************************** */


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

  // Free the partition functions
  free(Q);
  
  // Do a quick check of the free energies.  If any are > 0, it's likely there's
  // an input error.  Let the user know if this is the case.
  if (quiet == 0) {
    j = 0;
    while (j < cTotal && (*G)[j] <= 0.0001) {
      j++;
    }
    if (j < cTotal) {
      printf("\n\nWarning: At least one free energy is > 0. %lf\n", (*G)[j]);
      printf("It is likely there is an input error.\n");
      printf("If there is such an error, the the program will still run\n");
      printf("normally and give results, which may be nonsensical.\n");
      printf("If your input file suffix is .cx, .cx-epairs, or .cx-mfe,\n");
      printf("there should be no ordered complex identifier in the input file.\n");
      printf("\n\n");
    }
  }


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
  
  zerox0 = (int *) malloc(nSS * sizeof(int));
  nonzerox0 = (int *) malloc(nSS * sizeof(int));

  // Check to see if any of the concentrations are zero.  If so, 
  // reformulate the problem accordingly.
  j = 0;  // How many entries are zero
  k = 0;  // How many entries are nonzero
  for (i = 0; i < nSS; i++) {
    if ((*x0)[i] <= DBL_MIN) {
      zerox0[j] = i;
      j++;
    }
    else {
      nonzerox0[k] = i;
      k++;
    }
  }

  if (j > 0) { // Have to reformulate
  
    newnumSS = nSS - j;

    // First count how many entries we have.  notOK = 1 if complex contains something
    // with zero initial concentration
    newnumTotal = 0;
    
    for (j = 0; j < cTotal; j++) {
      notOK = 0;
      for (i = 0; i < nSS-newnumSS; i++) {
        if ((*A)[zerox0[i]][j] > 0) {
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
    newx0 = (double *) malloc(newnumSS * sizeof(double));
    newCompIDArray = (int *) malloc(newnumTotal * sizeof(int));
    newPermIDArray = (int *) malloc(newnumTotal * sizeof(int));
    
    // Put in the new x0
    for (i = 0; i < newnumSS; i++) {
      newx0[i] = (*x0)[nonzerox0[i]];
    }
    
    // Go through and pick out the entries to keep
    k = 0;
    for (j = 0; j < cTotal; j++) {
      notOK = 0;
      for (i = 0; i < nSS-newnumSS; i++) {
        if ((*A)[zerox0[i]][j] > 0) {
          notOK = 1;
        }
      }
      if (notOK == 0) {
        for (i = 0; i < newnumSS; i++) {
          newA[i][k] = (*A)[nonzerox0[i]][j];
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
    
    // Rename newx0
    free(*x0);
    (*x0) = (double *) malloc(newnumSS * sizeof(double));
    for (i = 0; i < newnumSS; i++) {
      (*x0)[i] = newx0[i];
    }
    free(newx0);
    
    // Rename numTotal and numSS
    *numTotal = newnumTotal;
    *numSS = newnumSS;
  
  }
  free(zerox0);
  free(nonzerox0);
  /* ************** FINISHED REFORMATTING PROBLEM *********************** */

  // Calculate molarity of water and convert appropriate quantities to the right units
  MolesWaterPerLiter = WaterDensity((*kT)/kB - ZERO_C_IN_KELVIN);
  for (i = 0; i < (*numSS); i++) {
    (*x0)[i] /= MolesWaterPerLiter;
  }

  return MolesWaterPerLiter;

}
/* ******************************************************************************** */

