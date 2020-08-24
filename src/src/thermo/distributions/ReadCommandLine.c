/*
  ReadCommandLine.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  Reads command line input for Distributions.c.  The argument is the
  prefix for the files that contain the input data, e.g., prefix.cx
  and prefix.count.  Uses the package getopt.h to retrieve option
  flags.

  Justin Bois, Caltech, 3 September 2006
  bois@caltech.edu
*/

#include "ReadCommandLine.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <getopt.h> // Takes options from the command line

#include "constants.h"

/* ******************************************************************************** */
void ReadCommandLine(int nargs, char **args, char *cxFile, char *countFile, 
                     char *logFile, char *distFile, char *LambdaFile, 
                     int *SortOutput, int *WriteLambda, double *MaxSizeLambda,
                     double *kT, int *quiet, int *WriteLogFile, int *Toverride,
                     int *NoPermID, int *NUPACK_VALIDATE, int *v3) {

  int options;  // Counters used in getting flags
  int ShowHelp; // ShowHelp = 1 if help option flag is selected
  char prefix[MAXLINE]; // The name of the prefix for the cmat, Fmat, and ccon files
  char InputStr[MAXLINE]; // Dummy string for storing flag options from command line
  FILE *fp; // The cx file, used to check if we can open it.

  if (nargs == 1) {
    printf("For instructions on running this program, run it with the ");
    printf("-help flag.\n\nExiting....\n\n");
    exit(ERR_NOINPUT);
  }

  *NoPermID = 0; // Default is to use .ocx file
  *SortOutput = 1;  // Default is to sort the output by concentration of complexes
  *WriteLogFile = 0; // Default is not to write log file
  *WriteLambda = 0; // Default is not to write out Lambda
  *MaxSizeLambda = 1e6; // Default maximum lambda size
  *kT = kB*(37.0 + ZERO_C_IN_KELVIN); // Default temperature is 37 deg. C
  *quiet = 0; // Default is to show messages on the screen
  *Toverride = 0; // Default is to either use T = 37 or that specified in input file
  *NUPACK_VALIDATE = 0;
  ShowHelp = 0;

  /* version 3 output */
  *v3 = 0;
  int prev_ordered = 1; // used to keep -ordered off if not specified before -v3

  SetExecutionPath(nargs, args);
  
  // Get the option flags
  while (1)
    {
      static struct option long_options [] =
      {
        {"sort",          required_argument,  0, 'a'},
        {"T",             required_argument,  0, 'b'},
        {"maxstates",     required_argument,  0, 'c'},
        {"quiet",         no_argument,        0, 'd'},
        {"help",          no_argument,        0, 'e'},
        {"writestates",   no_argument,        0, 'f'},
        {"writelogfile",  no_argument,        0, 'g'},
        {"ordered",       no_argument,        0, 'h'},
        {"validate",      no_argument,        0, 'i'},
        {"v3.0",          no_argument,        0, 'z'},
        {0, 0, 0, 0}
      };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      options = getopt_long_only (nargs, args,
                "a:b:c:defghi", long_options,&option_index);

      // Detect the end of the options.
      if (options == -1)
        break;

      switch (options)
        {
        case 'a':
	        strcpy(InputStr,optarg);
	        *SortOutput = atoi(InputStr);
          break;

	      case 'b':
	        strcpy(InputStr,optarg);
	        *kT = kB*(str2double(InputStr) + ZERO_C_IN_KELVIN);
	        *Toverride = 1;
	        break;

        case 'c':
	        strcpy(InputStr,optarg);
	        (*MaxSizeLambda) = str2double(InputStr);
          break;

	      case 'd':
	        *quiet = 1;
	        break;

	      case 'e':
	        ShowHelp = 1;
	        break;

	      case 'f':
	        *WriteLambda = 1;
	        break;

	      case 'g':
	        *WriteLogFile = 1;
	        break;

	      case 'h':
	        *NoPermID = 0;
          prev_ordered = 0;
	        break;

        case 'i':
          *NUPACK_VALIDATE = 1;
          *SortOutput = 3;
          break;

        case 'z':
          *v3 = 1;
          *NoPermID = prev_ordered;
          break;

        case '?':
          // getopt_long already printed an error message.
          break;

        default:
          abort ();
        }
    }

  /* version 3 output */
  if (!*quiet) {
    print_deprecation_info(stdout);
  }

  if (ShowHelp) {
    DisplayDistributionsHelp(1);
    exit(ERR_HELP);
  }

  if ((*MaxSizeLambda) > INT_MAX) {
    if (*quiet == 0) {
      printf("The maximum size of Lambda must be less than %d\n",INT_MAX);
      printf("This is the maximal allowed size of an integer.\n");
      printf("Using MaxSizeLambda = %d.\n",INT_MAX-1);
    }
    (*MaxSizeLambda) = INT_MAX-1;
  }

  if (*SortOutput > 4 || *SortOutput < 1) {
    if (*quiet == 0) {
      printf("Sorting option is an integer from 1 to 4.  Output will be sorted by\n");
      printf("average count of complexes.\n\n");
    }
    *SortOutput = 1;
  }

  // Get the the input file
  if (optind == nargs) { // There's no input from the user
    if (*quiet == 0) {
      printf("You must have a prefix or an input file on the command line.\n");
      printf("For instructions on running this program, run it with the ");
      printf("-help flag.\n\nExiting....\n\n");
    }
    exit(ERR_NOINPUT);
  }
  else {
    strcpy(prefix,args[optind]);
  }



  // Adjust sorting option from 2 to 1 if NoPermID == 1, since they're the same
  if (*NoPermID == 1 && *SortOutput == 2) {
    *SortOutput = 1;
  }

  // Name the files
  strcpy(cxFile,prefix);
  strcpy(countFile,prefix);
  strcpy(logFile,prefix);
  strcpy(distFile,prefix);
  strcpy(LambdaFile,prefix);
  if (*NoPermID) {
    strcat(cxFile,".cx");
  }
  else {
    strcat(cxFile,".ocx");
  }
  strcat(countFile,".count");
  strcat(logFile,".log");
  strcat(distFile,".dist");
  strcat(LambdaFile,".states");

  // Do a quick check to make sure the cx file exists before we proceed
  if ((fp = fopen(cxFile,"r")) == NULL) {
    if (*quiet == 0) {
      printf("Error opening %s!\n\nExiting....\n",cxFile);
    }
    exit(ERR_CX);
  }
  fclose(fp);

} 
/* ******************************************************************************** */


/* ******************************************************************************** */
void DisplayDistributionsHelp(int DummyArgument) {
  /*
    Displays the contents of the Distributions help file.
  */

  printf("Please read the NUPACK User Guide for detailed instructions.\n");
  printf("Usage: distributions [OPTIONS] PREFIX\n");
  printf("Calculate concentrations for each complex specified\n");
  printf("Options:\n");
  printf(" -ordered             perform the calculation on ordered complexes\n");
  printf(" -maxstates BIG       maximum number of states to enumerate\n");
  printf(" -writestates         write an output file describing properties for\n");
  printf("                      all population states of the system\n");
  printf(" -sort METHOD         change the sort method for the .eq output file\n");
  printf("                      0: same as input\n");
  printf("                      1: sort by concentration of complex; if -ordered\n");
  printf("                         specified, then sort by concentration of ordered\n");
  printf("                         complex\n");
  printf("                      2: sort by concentration of complex and then, if\n");
  printf("                         if -ordered is selected, by the concentration\n");
  printf("                         each constituent ordered complex\n");
  printf("                      3: output is sorted by complex/ordered complex id\n");
  printf("                      4: output is sorted by strands in the complex,\n");
  printf("                         then by strand quantities, then by ordered\n");
  printf("                         complex identifier\n");
  printf(" -quiet               suppress output to the screen\n");
  printf("\n");


  // Future versions will have detailed header file.
  /*
  char *NUPACKHome; // Home directory of NUPACK, were the help file is
  char HelpFileName[MAXLINE];
  char HelpBuffer[MAXLINE];
  FILE *fp;

  if ((NUPACKHome = getenv("NUPACKHOME")) != NULL) {

    strcpy(HelpFileName,NUPACKHome);
    strcat(HelpFileName,"/");
    strcat(HelpFileName,DISTRIBUTIONS_HELP_FILE);
    
    // Open the cx file
    if ((fp = fopen(HelpFileName,"r")) == NULL) {
      printf("Error in opening %s file!\n",HelpFileName);
      printf("\nExiting....\n\n");
      exit(ERR_HELP_OPEN);
    }
    
    while (fgets(HelpBuffer,MAXLINE,fp) != NULL) {
      printf("%s",HelpBuffer);
    }
  }
  else {
    printf("\nError: Environment variable NUPACKHOME not properly set.\n\n");
    printf("Exiting....\n\n");
    exit(ERR_HELP_OPEN);
  }

  fclose(fp);
  */
}
/* ******************************************************************************** */
void print_deprecation_info(FILE *out) {
  char *dep_mess = 
  "Relative to NUPACK 3.0, the following change was introduced to the\n"
  "distributions executable: the -ordered option is on by default.\n"
  "Use the -v3.0 option to revert to NUPACK 3.0 behavior.\n\n";

  fprintf(out, "%s", dep_mess);
  
  return;
}
/* ******************************************************************************** */
