/*
  ReadCommandLine.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  Reads command line input for Concentrations.c.  The argument is the
  either the prefix for the files that contain the input data, e.g.,
  prefix.cx and prefix.con, or it is the full file name of the input,
  e.g., prefix.ox-epairs, in which case the concentration information
  is if prefix.con.  Uses the package getopt.h to retrieve option
  flags.

  The command line option flags are
  -sort [required argument]
     sort = 0: the output is not sorted, but listed in the
       order in which it was inputted.
     sort = 1 (default for method = 1): the output is sorted first by 
       permutation concentration.
     sort = 2: the output is sorted first by complex ID number and then
       by permutation ID number.
     sort = 3: the output is sorted first by complex concentration.  The 
       permutation concentrations within each complex are then sorted.
     sort = 4: the output is sorted first by size of complexes and then 
       "alphabetically" within each complex size.  I.e., say there are 
       two strands, A, and B.  The ordering for complex size 2 is AA, AB, 
       BB.  The permutations are sorted within each complex but permutation
       ID number.
  -T [required argument] 
     The temperature, in degrees C, at which the calculation will be
     done.  Note that the partition functions used in the calculation
     should be calculated at the same temperature.  The default is T =
     37.
  -quiet [no argument]
     Selecting this flag will supress output to the screen.
  -cutoff [required argument]
     The cutoff value for reporting ensemble pair fractions (if -pairs is selected).
     Default is 0.001.
  -maxiters [required argument]
     The maximum number of iterations allowed in the trust region
     minimization algorithm.  Default is maxiters = 10000.
  -tol [required argument]
     The tolerance for converging to the equilibrium concentrations.
     This is entered as a fraction of the input concentration of the 
     single-stranded species.  I.e., if we have two strands A and B 
     that form complexes, and the input concentration of A is 1e-6 
     molar and B is 2e-6 molar, a value of tol = 0.001 means that 
     the maximal error in the conservation of mass for strand A is 
     0.001*1e-6 molar and 0.001*2e-6 for strand B.  The default value 
     is tol = 0.0000001.
  -maxtrial [required argument]
     The maximum number of initial conditions to try before stopping.
     The default is maxtrial = 100.  Not relevant for method = 2.
  -maxnostep [required argument]
     The maximum number of iterations in the trust region optimization
     that may be tried without a step being taken.  This corresponds
     to a very small trust region due to precision issues.  If
     maxnostep is exceeded, a new initial condition is tried.  The
     default is maxnostep = 50.  Not relevant for method = 2.
  -perturbscale [required argument]
     When new initial conditions are generated, they are perturbed
     from the standard initial condition by a prefactor multiplied
     by a random number from -1 to 1.  This prefactor is perturbscale.
     Its value is adjusted in the program in order to make sure there
     are no overflow errors.  The default is perturbscale = 100.  Not
     relevant for method = 2.
  -writelogfile [no argument]
     When this flag is selected, a log file is written which contains
     information about the trust region steps and convergence.  Default
     is not to write a log file.
  -pairs [no argument]
     Computes pair probabilities as well.
  -seed [required argument]
     Give a seed to the random number generator
  -help [no argument]
     Prints the help file to the screen.
  -validate [no argument]
     Prints all floats to 14 decimal places

  Justin Bois, Caltech, 2 September 2006
*/

#include "ReadCommandLine.h"
#include "constants.h"

/* ******************************************************************************** */
void ReadCommandLine(int nargs, char **args, char *cxFile, char *conFile, 
		     char *logFile, char *eqFile, char *pairsFile, char *fpairsFile,
		     int *SortOutput, int *MaxIters, double *tol, double *kT,
		     int *MaxNoStep, int *MaxTrial, double *PerturbScale, int *quiet,
		     int *WriteLogFile, int *Toverride, int *NoPermID, 
		     int *DoBPfracs, unsigned long *seed, double *cutoff,int * NUPACK_VALIDATE, int *v3) {

  int options;  // Counters used in getting flags
  int ShowHelp; // ShowHelp = 1 if help option flag is selected
  // int NoSortOutputOption; // No sorting option given from command line
  char prefix[MAXLINE]; // The prefix for the input and output files
  char InputStr[MAXLINE]; // Dummy string for storing flag options from command line
  FILE *fp; // The cx file, used to check if we can open it.

  if (nargs == 1) {
    printf("For instructions on running this program, run it with the ");
    printf("-help flag.\n\nExiting....\n\n");
    exit(ERR_NOINPUT);
  }

  // NoSortOutputOption = 1; // No sorting option given from command line
  *SortOutput = 1;  // Default is to sort the output by concentration
  *MaxIters = 10000; // Default maxiters
  *tol = 0.0000001; // Default tolerance is 0.00001% of minimum of the
                    // minimum count among the single-strands
  *kT = kB*(37.0 + ZERO_C_IN_KELVIN); // Default temperature is 37 deg. C
  *quiet = 0; // Default is to show messages on the screen
  *MaxNoStep = 50; // Default is 50 iterations with no step
  *MaxTrial = 100000; // Default is maximum of 100,000 trials
  *PerturbScale = 100; // Default is a scale of 100 on the perturb scale
  *WriteLogFile = 0; // Default is not to write a log file
  *Toverride = 0; // Default is to either use T = 37 or that specified in input file
  *NoPermID = 0; // Default is to use .ocx file => permutation IDs in file
  *DoBPfracs = 0; // Default is not to generate fpairs file.
  *seed = 0; // Default is to seed off the clock.
  *cutoff = 0.001; // Default cutoff
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
	        {"maxiters",      required_argument,  0, 'b'},
	        {"tol",           required_argument,  0, 'c'},
	        {"T",             required_argument,  0, 'd'},
	        {"quiet",         no_argument,        0, 'e'},
	        {"maxtrial",      required_argument,  0, 'f'},
	        {"maxnostep",     required_argument,  0, 'g'},
	        {"help",          no_argument,        0, 'h'},
	        {"perturbscale",  required_argument,  0, 'i'},
	        {"writelogfile",  no_argument,        0, 'j'},
	        {"ordered",       no_argument,        0, 'k'},
	        {"pairs",         no_argument,        0, 'l'},
	        {"seed",          required_argument,  0, 'm'},
	        {"cutoff",        required_argument,  0, 'n'},
          {"validate",      no_argument,        0, 'o'},
          {"v3.0",          no_argument,        0, 'z'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;


      options = getopt_long_only (nargs, args, 
				  "a:b:c:d:ef:g:hi:jklm:n:o", long_options, 
				  &option_index);

      // Detect the end of the options.
      if (options == -1)
        break;

      switch (options)
        {
        case 'a':
	        strcpy(InputStr,optarg);
	        *SortOutput = atoi(InputStr);
	        // NoSortOutputOption = 0; // Record that we've selected a sorting option
          break;

	      case 'b':
	        strcpy(InputStr,optarg);
	        *MaxIters = atoi(InputStr);
	        break;

	      case 'c':
	        strcpy(InputStr,optarg);
	        *tol = str2double(InputStr);
	        break;

	      case 'd':
	        strcpy(InputStr,optarg);
	        *kT = kB*(str2double(InputStr) + ZERO_C_IN_KELVIN);
	        *Toverride = 1;
	        break;

	      case 'e':
	        *quiet = 1;
	        break;

	      case 'f':
	        strcpy(InputStr,optarg);
	        (*MaxTrial) = atoi(InputStr);
	        break;

        case 'g':
	        strcpy(InputStr,optarg);
	        (*MaxNoStep) = atoi(InputStr);
          break;

        case 'h':
	        ShowHelp = 1;
          break;

        case 'i':
	        strcpy(InputStr,optarg);
	        (*PerturbScale) = str2double(InputStr);
          break;

	      case 'j':
	        *WriteLogFile = 1;
	        break;

	      case 'k':
	        *NoPermID = 0;
          prev_ordered = 0;
	        break;

	      case 'l':
	        *DoBPfracs = 1;
	        break;

	      case 'm':
	        strcpy(InputStr,optarg);
	        (*seed) = (unsigned long)atol(InputStr);
	        break;

	      case 'n':
	        strcpy(InputStr,optarg);
	        (*cutoff) = str2double(InputStr);
	        break;

        case 'o':
          *NUPACK_VALIDATE = 1;
          *tol = 0.0000000000001;
	        *SortOutput = 3;
	        // NoSortOutputOption = 0; // Record that we've selected a sorting option
          *cutoff = 0;
          break;

        case 'z':
          *v3 = 1;
          *NoPermID = prev_ordered;
          break;

        case '?':
          // getopt_long already printed an error message.
          break;
        
        default:
          abort();
        }
    }

  /* version 3 output */
  if (!*quiet) {
    print_deprecation_info(stdout);
  }

  if (ShowHelp) {
    DisplayHelpConc();
    exit(ERR_HELP);
  }

  if (*SortOutput > 4) {
    if (*quiet == 0) {
      printf("Sorting option is an integer from 0 to 4.  Output will be sorted by\n");
      printf("concentration.\n\n");
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

  // Name the files
  strcpy(cxFile,prefix);
  strcpy(conFile,prefix);
  strcpy(logFile,prefix);
  strcpy(eqFile,prefix);
  strcpy(pairsFile,prefix);
  strcpy(fpairsFile,prefix);
  if (*NoPermID) {
    strcat(cxFile,".cx");
    strcat(pairsFile,".cx-epairs");
  }
  else {
    strcat(cxFile,".ocx");
    strcat(pairsFile,".ocx-epairs");
  }
  strcat(conFile,".con");
  strcat(logFile,".log");
  strcat(eqFile,".eq");
  strcat(fpairsFile,".fpairs");

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
void DisplayHelpConc() {
  /*
    Displays the contents of the Concentrations help file.
  */

  printf("Please read the NUPACK User Guide for detailed instructions.\n");
  printf("Usage: concentrations [OPTIONS] PREFIX\n");
  printf("Calculate concentrations for each complex specified\n");
  printf("Options:\n");
  printf(" -ordered             perform the calculation on ordered complexes\n");
  printf(" -pairs               compute base-pairing information for the solution\n");
  printf(" -cutoff CUTOFFVALUE  only store ensemble pair fractions above CUTOFFVALUE\n");
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

}
/* ******************************************************************************** */
void print_deprecation_info(FILE *out) {
  char *dep_mess = 
  "Relative to NUPACK 3.0, the following change was introduced to the\n"
  "concentrations executable: the -ordered option is on by default.\n"
  "Use the -v3.0 option to revert to NUPACK 3.0 behavior.\n\n";

  fprintf(out, "%s", dep_mess);
  
  return;
}
/* ******************************************************************************** */
