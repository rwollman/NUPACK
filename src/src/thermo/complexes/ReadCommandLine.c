/*
  ReadCommandLine.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks, 3/2006 and Justin Bois 1/2007

  Reads command line input for multiWrapper.c.
  Uses the package getopt.h to retrieve option flags.
*/

#include "ReadCommandLine.h"
#include <string.h>
#include <stdlib.h>
#include <getopt.h> // Takes options from the command line
#include <time.h>
#include <ctype.h>

#include <shared.h>
#include <thermo/core.h>
#include "complexesStructs.h"

extern globalArgs_t globalArgs;
extern char PARAM_FILE[];


/* ********** SCRIPT FOR DISPLAYING HELP ************************** */
void DisplayHelpComplexes(void);

/* ****************************************************************** */
int ReadCommandLine(int nargs, char **args) {

  // Returns 1 if input file is specified and 0 otherwise.

  int options;  // Counters used in getting flags
  int ShowHelp=0; // ShowHelp = 1 if help option flag is selected
  int option_index;
  char line[MAXLINE];
  char param[MAXLINE];
  float temp;

  /* version 3 output */
  int prev_ordered = 0; // used to keep -ordered off if not specified before -v3

  static struct option long_options [] =
    {
      {"quiet",         no_argument,        NULL, 'a'},
      {"ordered",       no_argument,        NULL, 'b'},
      {"out",           required_argument,  NULL, 'c'},
      {"pairs",         no_argument,        NULL, 'd'},
      {"T",             required_argument,  NULL, 'e'},
      {"dangles",       required_argument,  NULL, 'f'},
      {"material",      required_argument,  NULL, 'g'},
      {"help",          no_argument,        NULL, 'h'},
      {"timeonly",      no_argument,        NULL, 'i'},
      {"listonly",      no_argument,        NULL, 'w'},
      {"debug",         no_argument,        NULL, 'v'},
      {"echo",          no_argument,        NULL, 'j'},
      {"mfe",           no_argument,        NULL, 'k'},
      {"cutoff",        required_argument,  NULL, 'l'},
      {"progress",      no_argument,        NULL, 'm'},
      {"degenerate",    no_argument,        NULL, 'n'},
      {"sodium",        required_argument,  NULL, 'o'},
      {"magnesium",     required_argument,  NULL, 'p'},
      {"longhelixsalt", no_argument,        NULL, 'q'},
      {"validate",      no_argument,        NULL, 'r'},
      {"defect",        no_argument,        NULL, 's'},
      {"v3.0",          no_argument,        NULL, 'z'},
      {0, 0, 0, 0}
    };

  SetExecutionPath(nargs, args);


  NUPACK_VALIDATE=0;
  // Get the option flags
  while (1)
    {

      /* getopt_long stores the option index here. */
      option_index = 0;
      options = getopt_long_only (nargs, args, 
				  "abcde:f:g:h:ijkl:mno:p:qrs", 
                                  long_options, 
				  &option_index);
		
      // Detect the end of the options.
      if (options == -1)
        break;

      switch (options) {
      case 'a':
        globalArgs.quiet = 1;
        break;

      case 'b':
        globalArgs.permsOn = 1;
        prev_ordered = 1;
        break;

      case 'c':
        strcpy( line, optarg);
        if( sscanf(line, "%d", &(globalArgs.out)) != 1) {
          printf("Invalid out value\n");
          exit(1);
        }
        break;

      case 'd':
        globalArgs.dopairs = 1;
        break;

      case 'e':
        strcpy( line, optarg);
        if( sscanf(line, "%f", &(temp)) != 1) {
          printf("Invalid T value\n");
          exit(1);
        }
        globalArgs.T = (long double) temp;
        break;

      case 'f':
        strcpy( line, optarg);
        if( sscanf(line, "%s", param) != 1) {
          printf("Invalid parameters value\n");
          exit(1);
        }
        if (isdigit(param[0])) {
          globalArgs.dangles = atoi(param);
        }
        else {
          if (!strcmp(param,"none")) {
            globalArgs.dangles = 0;
          }
          else if (!strcmp(param,"some")) {
            globalArgs.dangles = 1;
          }
          else if (!strcmp(param,"all")) {
            globalArgs.dangles = 2;
          }
          else {
            printf("Invalid dangles value\n");
            exit(1);
          }
        }
        break;

      case 'g':
        strcpy( line, optarg);
        if( sscanf(line, "%s", param) != 1) {
          printf("Invalid parameters value\n");
          exit(1);
        }
        else if( !strcmp(param, "dna") || !strcmp(param, "dna1998")) {
          globalArgs.parameters = DNA;

        } else if( !strcmp(param, "rna") || !strcmp(param, "rna1995") ) {
          globalArgs.parameters = RNA;
        } else if( !strcmp(param, "rna37") || !strcmp(param, "rna1999") ) {
          if(strcmp(PARAM_FILE,"rna37") == 0) {
            printf("Parameter specification using rna37 has been deprecated. Please use rna1999 instead\n");
          }
          globalArgs.parameters = RNA37;
        } else {
          globalArgs.parameters = USE_SPECIFIED_PARAMETERS_FILE;
          strcpy( PARAM_FILE, param);
        }
        break;

      case 'h':
        ShowHelp = 1;
        break;

      case 'i':
        globalArgs.timeonly = 1;
        break;

      case 'w':
        globalArgs.listonly = 1;
        break;

      case 'v':
        globalArgs.debug = 1;
        break;

      case 'j':
        globalArgs.echo = 1;
        break;

      case 'k':
        globalArgs.mfe = 1;
        break;

      case 'l':
        strcpy( line, optarg);
        if( sscanf(line, "%f", &(temp)) != 1) {
          printf("Invalid cutoff value\n");
          exit(1);
        }
        globalArgs.cutoff = (long double) temp;
        break;

      case 'm':
        globalArgs.progress = 1;
        break;

      case 'n':
        globalArgs.onlyOneMFE = 0;
        break;

      case 'o':
	      strcpy( line, optarg);
	      globalArgs.sodiumconc = str2double(line);
	      break;
	
      case 'p':
	      strcpy( line, optarg);
	      globalArgs.magnesiumconc = str2double(line);
	      break;
	
      case 'q':
	      globalArgs.uselongsalt = 1;
	      break;

      case '?':
        // getopt_long already printed an error message.
        break;

      case 'r':
        NUPACK_VALIDATE=1;
        globalArgs.permsOn = 1;
        globalArgs.cutoff = 0.0;
        
        break;

      case 's':
        globalArgs.dodefect = 1;
        break;
      
      case 'z':
        globalArgs.v3 = 1;
        globalArgs.permsOn = prev_ordered;
        break;
      
      default:
        abort ();
      }
    }
      

  /* version 3 output */
  if (!globalArgs.quiet) {
    print_deprecation_info(stdout);
  }

  if (ShowHelp) {
    DisplayHelpComplexes();
    exit(1);
  }


  // Check salt inputs to make sure we're ok
  if ((globalArgs.sodiumconc != 1.0 || globalArgs.magnesiumconc != 0.0) && globalArgs.parameters != DNA) {
    printf("WARNING: No salt corrections availabe for RNA.  Using 1 M Na and 0 M Mg.\n");
    globalArgs.sodiumconc = 1.0;
    globalArgs.magnesiumconc = 0.0;
  }

  if (globalArgs.sodiumconc  <= 0.0) {
    printf("ERROR: In valid sodium concentration.  Must have [Na+] > 0.\n");
    exit(1);
  }

  if (globalArgs.magnesiumconc  < 0.0) {
    printf("ERROR: In valid magnesium concentration.  Must have [Mg2+] >= 0.\n");
    exit(1);
  }

  if (globalArgs.sodiumconc < 0.05 || globalArgs.sodiumconc > 1.1) {
    printf("WARNING: Salt correction only verified for 0.05 M < [Na+] < 1.1 M.\n");
    printf("         [Na+] = %g M may give erroneous results.\n",(float) globalArgs.sodiumconc);
  }

  if (globalArgs.magnesiumconc > 0.2) {
    printf("WARNING: Salt correction only verified for [Mg2+] <= 0.2 M.\n");
    printf("         [Mg2+] = %g M may give erroneous results.\n",(float) globalArgs.magnesiumconc);
  }
  // The range of validity of the magnesium correction is unknown

  if (globalArgs.uselongsalt && globalArgs.magnesiumconc > 0.0) {
    printf("WARNING: No magnesium correction parameters are available for the long\n");
    printf("         helix mode of salt correction.  Using [Mg2+] = 0.\n");
    globalArgs.magnesiumconc = 0.0;
  }
  

  // Get the the input file
  if (optind == nargs) { // There's no input from the user
    return 0;
  }
  else {
    strcpy(globalArgs.inputFilePrefix,args[optind]);
  }

  return 1;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void DisplayHelpComplexes() {
  printf("Please read the NUPACK User Guide for detailed instructions.\n");
  printf("Usage: complexes [OPTIONS] PREFIX\n");
  printf("Calculate equilibrium properties of all possible unpseudoknotted complexes\n");
  printf("of the input strands up to user-defined size L_max\n");
  PrintNupackThermoHelp();
  printf("Additional options:\n");
  printf(" -ordered         store properties for ordered complexes\n");
  printf(" -pairs           store base-pairing observables\n");
  printf(" -cutoff CUTOFF   set the minimum stored probability/expected value\n");
  printf(" -mfe             compute and store MFEs for all ordered complexes\n");
  printf(" -quiet           suppress output to the screen\n");
  printf("\n");
  // printf("Common options are:\n");
  // printf("\t-material [parameters]\n");
  // printf("\t\tSee the manual and references for parameter sources\n");
  // printf("\t\trna1995\n");
  // printf("\t\tdna1998\n");
  // printf("\t\trna1999\n"); 
  // printf("\t-T [temperature]\n");
  // printf("\t\tTemperature specified in degrees Celsius\n\n");
  // printf("\t-ordered\n");
  // printf("\t\tStore properties for ordered complexes.\n\n");
  // printf("\t-pairs\n");
  // printf("\t\tCalculate base-pairing observables as for the pairs executable.\n\n");
  // printf("\t-mfe\n");
  // printf("\t\tCalculate all minimum free energy structures for each ordered complex.\n");
  // printf("\t\tMust be used in conjunction with the -ordered flag.\n\n");
  // printf("\t-quiet\n");
  // printf("\t\tSuppress output to the screen\n");
}
/* ******************************************************************************** */

int ReadInputFileComplexes( char *filePrefix, int *nStrands,
                            char ***seqs, int **seqlength,
                            int *maxLength, int *maxComplexSize) {

  // Returns 1 is input file is ok, 0 if it has an error, and 2 if it is not found.

  FILE *F_inp;
  char *token;
  char line[MAXLINE];
  char line2[MAXLINE];
  int i;
  char filename[MAXLINE];

  strcpy( filename, globalArgs.inputFilePrefix);
  strcat( filename, ".in");

  strcpy( filePrefix, globalArgs.inputFilePrefix);

  F_inp = fopen( filename, "r");
  if( !F_inp) {
    return 2;
  }

  fgets( line, MAXLINE, F_inp);
  while( line[0] == '%' || line[0] == '>') {
    fgets( line, MAXLINE, F_inp);
  }

  //Read in # of strands
  token = strtok( line, "%>");
  if( sscanf( token, "%d", nStrands ) != 1) {
    printf("Error in %s: %s\n", filename, token);
    fclose( F_inp);
    return 0;
  }

  //allocate function variables
  (*seqs) = (char **) malloc( *nStrands * sizeof( char*));
  (*seqlength) = (int *) malloc( *nStrands * sizeof( int));

  //read in sequences
  *maxLength = 0;
  for( i = 0; i < *nStrands; i++) {
    fgets( line, MAXLINE, F_inp);
    while( line[0] == '%' || line[0] == '>') {
      fgets( line, MAXLINE, F_inp);
    }
    token = strtok( line, "%>");

    if( sscanf( token, "%s", line2 ) != 1) {
      printf("Error in %s: %s\n", filename, token);
      fclose( F_inp);
      return 0;
    }
    (*seqlength)[i] = strlen( line2);
    if( (*seqlength)[i] > *maxLength) *maxLength = (*seqlength)[i];


    (*seqs)[i] = (char*) malloc( ((*seqlength)[i]+1)*sizeof( char));

    strcpy( (*seqs)[i], line2);
  }

  //read in Lmax
  fgets( line, MAXLINE, F_inp);
  while( line[0] == '%' || line[0] == '>') {
    fgets( line, MAXLINE, F_inp);
  }
  token = strtok( line, "%>");

  if( sscanf( token, "%d", maxComplexSize ) != 1) {
    printf("Error in %s: %s\n", filename, token);
    fclose( F_inp);
    return 0;
  }

  //check if complex list size is included
  char * check = fgets(line, MAXLINE, F_inp);
  if(!check) {
    fclose( F_inp);
    return 1;
  }
  while(check && (line[0] == '%' || line[0] == '>') ) {
    check = fgets( line, MAXLINE, F_inp);
  }
  if (!check) {
    fclose( F_inp);
    return 1;
  }


  fclose( F_inp);
  return 1;
}

/* ******** */
void printHeader( int nStrands, char **seqs, int maxComplexSize,
                  int nTotalOrders, int nNewPerms, int nSets, int nNewComplexes,
                  FILE *F_cx, int nargs, char **argv, int isPairs) {

  char comment[2] = "%";
  char timestr[25];
  time_t curtime;
  struct tm *loctime;
  int i;
  int initial_perms = nTotalOrders - nNewPerms;
  
  curtime = time(NULL); //current time
  loctime = localtime( &curtime);


  fprintf(F_cx,"%s NUPACK %s\n",comment, NUPACK_VERSION);
  fprintf(F_cx,"%s Program: complexes\n",comment);

  strncpy(timestr,asctime(loctime),24);
  timestr[24] = '\0';

  fprintf( F_cx, "%s Start time: %s PST\n%s\n", comment,
           timestr, comment);

  fprintf( F_cx, "%s Command: ", comment);
  for( i = 0; i < nargs; i++) {
    fprintf( F_cx, "%s ", argv[ i]);
  }

  fprintf( F_cx, "\n");

  fprintf( F_cx, "%s Maximum complex size to enumerate: %d\n",
           comment, maxComplexSize);

  // Show that we're using a cutoff if this is a pairs file
  if (isPairs && globalArgs.cutoff > 0.0) {
    fprintf(F_cx, "%s Minimum output pair probability: %Lg\n",
            comment, globalArgs.cutoff);
  }
  if (globalArgs.v3) {
    fprintf( F_cx, "%s Number of complexes from enumeration: %d\n",
             comment, nSets);
    fprintf( F_cx, "%s Additional complexes from .list file: %d\n",
             comment, nNewComplexes);
    fprintf( F_cx, "%s Total number of permutations to calculate: %d\n",
             comment, nTotalOrders );
  }
  else {
    fprintf( F_cx, "%s Number of complexes from enumeration: %d\n",
             comment, initial_perms);
    fprintf( F_cx, "%s Additional complexes from .list file: %d\n",
             comment, nNewPerms);
    fprintf( F_cx, "%s Total number of complexes: %d\n",
             comment, nTotalOrders );
  }

  fprintf( F_cx, "%s Parameters: ", comment);

  if( globalArgs.parameters == DNA) {
    fprintf( F_cx, "DNA, 1998\n");
  }
  else if(  globalArgs.parameters == RNA) {
    fprintf( F_cx, "RNA, 1995\n");
  }
  else if( globalArgs.parameters == RNA37) {
    fprintf( F_cx, "RNA, 1999\n");
  }

  fprintf( F_cx, "%s Dangles setting: %d\n", comment,
           globalArgs.dangles);

  fprintf( F_cx, "%s Temperature (C): %.1f\n", comment,
           (float) globalArgs.T);

  fprintf(F_cx,"%s Sodium concentration: %.4f M\n",comment, (float) globalArgs.sodiumconc);
  fprintf(F_cx,"%s Magnesium concentration: %.4f M\n",comment,(float) globalArgs.magnesiumconc);
  fprintf( F_cx, "%s\n", comment);
  fprintf( F_cx, "%s Do not change the comments below this line, as they may be read by other programs!",
           comment);
  fprintf( F_cx, "\n%s\n%s Number of strands: %d\n", comment, comment,
           nStrands);
  fprintf( F_cx, "%s id sequence\n", comment);
  for( i = 0; i < nStrands; i++) {
    fprintf( F_cx, "%s %2d %s\n", comment, i+1, seqs[i]);
  }

  fprintf( F_cx, "%s T = %.1f\n", comment, (float) globalArgs.T);
}

/* ******** */
void print_deprecation_info(FILE *out) {
  char *dep_mess = 
  "Relative to NUPACK 3.0, the following changes were introduced to\n"
  "the complexes executable:\n"
  "  -ordered is on by default\n"
  "  output files .cx and .cx-epairs are not written\n"
  "  the comment lines in .ocx-epairs and .ocx-mfe employ updated terminology\n"
  "Use the -v3.0 option to revert to NUPACK 3.0 behavior.\n\n";

  fprintf(out, "%s", dep_mess);
  
  return;
}
