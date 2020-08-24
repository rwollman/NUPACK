/*
  ReadCommandLineNPK.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Robert Dirks 9/2006, Justin Bois 1/2007

  Uses the package getopt.h to retrieve option flags.

  The command line option flags are

  -T [float] (default = 37 C)
  Sets the temperature to the argument

  -dangles [int]
  Sets the dangle option to the argument.
  0 = no dangles
  1 = normal dangles (mfold) (default)
  2 = vienna D2 dangles (always include dangle energy regardless of nearby structure)

  -params [string]
  look for the energy files with the specified prefix, first looking locally, then in the
  NUPACKHOME directory (if the NUPACKHOME environment variable has been set)

  -mfe [no argument]
  Enable mfe calculations for permutation.
  Warning, there is no guarantee that this can be accomplished in polynomial time.
  Due to the necessity of symmetry corrections.

  -degenerate [no argument]
  Output all degenerate mfe structures

  -samples [int]
  number of samples to take in sample executable

  -seed [int]
  the seed to use for sampling

  -validate
  print everything to 14 decimal places. print all pairs.
*/

#include "ReadCommandLineNPK.h"

double CUTOFF; // = 0.001 by default and can be changed by user
int Multistranded; //  = 1 if this is a multistranded calculation
int seqlengthArray[MAXSTRANDS]; // Length of sequences
int perm[MAXSTRANDS];  // Perm IDs
int nUniqueSequences; // Number of unique sequences entered
#ifdef NUPACK_SAMPLE
int nupack_random_seed; // Seed to use for sampling
#endif // NUPACK_SAMPLE

/* ****************************************************************** */
int ReadCommandLineNPK(int nargs, char **args, char *inputFile) {
  // Returns 1 if input file is specified and 0 otherwise

  int options;  // Counters used in getting flags
  int ShowHelp=0; // ShowHelp = 1 if help option flag is selected
  int batch=0; // Whether to prompt for input if incomplete options supplied
  int option_index;
  char line[MAXLINE];
  char param[MAXLINE];
  float temp;

  static struct option long_options [] =
    {

      {"multi", no_argument,          NULL, 'b'},
      {"pseudo", no_argument,          NULL, 'd'},
      {"degenerate", no_argument,     NULL,'c'},
      {"T", required_argument,          NULL, 'e'},
      {"dangles", required_argument,    NULL, 'f'},
      {"material", required_argument,     NULL, 'g'},
      {"help", no_argument,           NULL, 'h'},
      {"cutoff", required_argument,   NULL,'i'},
      {"batch", no_argument,   NULL,'t'},
      {"sodium", required_argument,   NULL, 'j'},
      {"magnesium", required_argument,   NULL, 'k'},
      {"longhelixsalt", no_argument,    NULL, 'l'},
      {"mfe", no_argument,    NULL, 'm'},
      #ifdef NUPACK_SAMPLE
      {"samples",required_argument, NULL, 'n'},
      {"seed",required_argument,NULL,'o'},
      #endif //NUPACK_SAMPLE
      {"sort",required_argument,NULL,'p'},
      {"validate",no_argument,NULL,'q'},
      {0, 0, 0, 0}
    };

  SetExecutionPath(nargs, args);

  //initialize global parameters
  DNARNACOUNT = RNA;
  DANGLETYPE = 1;
  TEMP_K = 37.0 + ZERO_C_IN_KELVIN;

  DO_PSEUDOKNOTS = 0;
  ONLY_ONE_MFE = 1;
  USE_MFE = 0;
  mfe_sort_method = 0;
  Multistranded = 0;
  CUTOFF = 0.001;
  SODIUM_CONC = 1.0;
  MAGNESIUM_CONC = 0.0;
  USE_LONG_HELIX_FOR_SALT_CORRECTION = 0;
  NUPACK_VALIDATE=0;
  EXTERN_QB = NULL;
  EXTERN_Q = NULL;


  // Get the option flags
  while (1) {
    /* getopt_long stores the option index here. */
    option_index = 0;
    options = getopt_long_only (nargs, args,
                                "bcde:f:g:hti:j:k:lm", long_options,
                                &option_index);

    // Detect the end of the options.
    if (options == -1)
      break;

    switch (options) {

    case 'b':
      Multistranded = 1;
      break;

    case 'c':
      ONLY_ONE_MFE = 0;
      break;

    case 'd':
      DO_PSEUDOKNOTS = 1;
      break;

    case 'e':
      strcpy( line, optarg);
      if( sscanf(line, "%f", &(temp)) != 1) {
        printf("Invalid T value\n");
        exit(1);
      }
      TEMP_K = (DBL_TYPE) (temp + ZERO_C_IN_KELVIN);
      break;

    case 'f':
      strcpy( line, optarg);
      if( sscanf(line, "%s", param) != 1) {
        printf("Invalid parameters value\n");
        exit(1);
      }
      if (isdigit(param[0])) {
        DANGLETYPE = atoi(param);
      }
      else {
        if (!strcmp(param,"none")) {
          DANGLETYPE = 0;
        }
        else if (!strcmp(param,"some")) {
          DANGLETYPE = 1;
        }
        else if (!strcmp(param,"all")) {
          DANGLETYPE = 2;
        }
        else {
          printf("Invalid dangles value\n");
          exit(1);
        }
      }
      break;

    case 'g':
      strcpy( line, optarg);
      if( sscanf(line, "%s", PARAM_FILE) != 1) {
        printf("Invalid parameters value\n");
        exit(1);
      }
      if( strcmp( PARAM_FILE, "rna") == 0 || strcmp( PARAM_FILE, "rna1995") == 0) {
        DNARNACOUNT = RNA;
      }
      else if( strcmp( PARAM_FILE, "dna") == 0 || strcmp( PARAM_FILE, "dna1998") == 0) {
        DNARNACOUNT = DNA;
      }
      else if( strcmp( PARAM_FILE, "rna37") == 0 || strcmp( PARAM_FILE, "rna1999") == 0) {
        if(strcmp(PARAM_FILE,"rna37") == 0) {
          printf("Parameter specification rna37 has been deprecated, please use the name rna1999.\n");
        }
        DNARNACOUNT = RNA37;
      }
      else {
        DNARNACOUNT = USE_SPECIFIED_PARAMETERS_FILE;
      }
      break;

    case 'h':
      ShowHelp = 1;
      NupackShowHelp = ShowHelp;
      return 0;
      break;

    case 'i':
      strcpy( line, optarg);
      if( sscanf(line, "%f", &(temp)) != 1) {
        printf("Invalid cutoff\n");
        exit(1);
      }
      CUTOFF = temp;
      break;

    case 't':
      batch=1;
      break;

    case 'j':
      strcpy( line, optarg);
      SODIUM_CONC = str2double(line);
      break;

    case 'k':
      strcpy( line, optarg);
      MAGNESIUM_CONC = str2double(line);
      break;

    case 'l':
      USE_LONG_HELIX_FOR_SALT_CORRECTION = 1;
      break;

    case 'm':
      USE_MFE = 1;
      break;
#ifdef NUPACK_SAMPLE
    case 'n':
      strcpy(line,optarg);
      if (isdigit(line[0])) {
        nupack_num_samples = atoi(line);
      } else {
        printf("Invalid Number of Samples Specified");
      }
      break;
    case 'o':
      strcpy(line,optarg);
      if (isdigit(line[0])) {
        nupack_random_seed = atoi(line);
      } else {
        printf("Invalid seed provided");
      }
      break;
  #endif //NUPACK_SAMPLE
    case 'p':
      strcpy(line,optarg);
      if (isdigit(line[0])) {
        mfe_sort_method = atoi(line);
      } else {
        printf("Invalid Sort Method\n");
        mfe_sort_method = 0;
      }
      if(mfe_sort_method != 0 && mfe_sort_method != 1) {
        printf("Invalid Sort Method\n");
        mfe_sort_method = 0;
      }
      break;
    case 'q':
      NUPACK_VALIDATE = 1;
      mfe_sort_method = 1;
      CUTOFF=0.0;
      break;
    default:
      abort ();
    }
  }

  // Check salt inputs to make sure we're ok
  if ((SODIUM_CONC != 1.0 || MAGNESIUM_CONC != 0.0) && DNARNACOUNT != DNA) {
    printf("%% ************************************************************************  %%\n");
    printf("%% WARNING: No salt corrections availabe for RNA.  Using 1 M Na and 0 M Mg.  %%\n");
    printf("%% ************************************************************************  %%\n");
    SODIUM_CONC = 1.0;
    MAGNESIUM_CONC = 0.0;
  }

  if (SODIUM_CONC  <= 0.0) {
    printf("ERROR: Invalid sodium concentration.  Must have [Na+] > 0.\n");
    exit(1);
  }

  if (MAGNESIUM_CONC  < 0.0) {
    printf("ERROR: Invalid magnesium concentration.  Must have [Mg2+] >= 0.\n");
    exit(1);
  }

  if (SODIUM_CONC < 0.05 || SODIUM_CONC > 1.1) {
    printf("%% ************************************************************************  %%\n");
    printf("%%    WARNING: Salt correction only verified for 0.05 M < [Na+] < 1.1 M.     %%\n");
    printf("%%                  [Na+] = %.2f M may give erroneous results.               %%\n",(float) SODIUM_CONC);
    printf("%% ************************************************************************  %%\n");
  }

  if (MAGNESIUM_CONC > 0.2) {
    printf("%% ************************************************************************  %%\n");
    printf("%%    WARNING: Salt correction only verified for [Mg2+] < 0.2 M.     %%\n");
    printf("%%                  [Mg2+] = %.2f M may give erroneous results.               %%\n",(float) MAGNESIUM_CONC);
    printf("%% ************************************************************************  %%\n");
  }
  // The range of validity of the magnesium correction is unknown

  if (USE_LONG_HELIX_FOR_SALT_CORRECTION && MAGNESIUM_CONC > 0.0) {
    printf("%% ************************************************************************  %%\n");
    printf("%%  WARNING: No magnesium correction parameters are available for the long   %%\n");
    printf("%%              helix mode of salt correction.  Using [Mg2+] = 0.            %%\n");
    printf("%% ************************************************************************  %%\n");
    MAGNESIUM_CONC = 0.0;
  }

  // Get the the input file
  if (optind == nargs) { // There's no input from the user
    printf("No input file specified.\n");
    if (!batch) {
      printf("Requesting input manually.\n");
      return 0;
    }
    else {
      printf("Exiting in batch mode.\n");
      abort();
    }
  }
  else {
    strcpy(inputFile,args[optind]);
    strcat(inputFile,".in");
  }


  return 1+batch;

}
/* ******************************************************************************** */


void PrintNupackThermoHelp() {
  printf("Thermodynamic parameters:\n");
  printf(" -material MATERIAL           set the material to MATERIAL. Standard\n");
  printf("                              choices are rna1995, dna1998, rna1999.\n");
  printf(" -sodium CONCENTRATION        set the concentration of sodium to\n");
  printf("                              CONCENTRATION\n");
  printf(" -dangles TREATMENT           specify treatment of dangle energies\n");
  printf("                              none, some, or all\n");
  printf(" -T TEMPERATURE               set the temperature to TEMPERATURE\n");
  printf("\n");
}

/* ******************************************************************************** */

void PrintNupackUtilitiesHelp() {
  printf("Calculation specifications:\n");
  printf(" -multi                       specify calculation involving complexes\n");
  printf("                              of multiple strands\n");
  printf(" -pseudo                      include a subclass of pseudoknots\n");
  printf("\n");
}

/* ******************************************************************************** */

int ReadInputFile( char *inputFile, char *theseq, int *v_pi, float *gap, char *structure,
       int *thepairs) {

  FILE *F_inp;
  char *token = NULL;
  char line[MAXLINE];
  char line2[MAXLINE];
  int i, j, permSize = 0;
  char **seqs = NULL;
  int base1 = 0;
  int base2 = 0;
  int nStrands;
  char *lastline = NULL;


  F_inp = fopen( inputFile, "r");
  if( !F_inp) {
    printf("Failed to open input file %s\nRequesting input manually.\n",inputFile);
    return 0;
  }


  lastline = fgets( line, MAXLINE, F_inp);
  while( lastline && (line[0] == '%' || line[0] == '>') ) {
    //ignore comments
    lastline = fgets( line, MAXLINE, F_inp);
  }
  if( !lastline) {
    printf("Error in %s: %s\nRequesting input manually.\n", inputFile, token);
    fclose( F_inp);
    return 0;
  }


  token = strtok( line, "%>");

  if( !Multistranded) {
    if( sscanf( token, "%s", theseq) == 1)  {
      if (isalpha(theseq[0]) == 0) {
        printf("Error in %s: Perhaps you need the -multi flag.\n",inputFile);
        printf("Requesting input manually.\n");
        fclose(F_inp);
        return 0;
      }
      *v_pi = 1;
    }
    else {
      printf("Error in %s: %s\nRequesting input manually.\n", inputFile, token);
      fclose( F_inp);
      return 0;
    }
  }
  else {
    if( sscanf( token, "%d", &nStrands ) != 1) {
      printf("Error in %s: %s\nRequesting input manually.\n", inputFile, token);
      fclose( F_inp);
      return 0;
    }

    //allocate function variables
    seqs = (char **) malloc( nStrands * sizeof( char*));

    //read in sequences
    for( i = 0; i < nStrands; i++) {
      lastline = fgets( line, MAXLINE, F_inp);
      while( lastline && (line[0] == '%' || line[0] == '>') ) {
        lastline = fgets( line, MAXLINE, F_inp);
      }
      token = strtok( line, "%>");

      if( !lastline || sscanf( token, "%s", line2 ) != 1) {
        printf("Error in %s: %s\nRequesting input manually.\n", inputFile, token);
        fclose( F_inp);
        for( j = 0; j < i; j++) free( seqs[j]);
        free( seqs);
        return 0;
      }

      seqs[i] = (char*) malloc( (strlen(line2)+1)*sizeof( char));

      strcpy( seqs[i], line2);

      seqlengthArray[i] = strlen(seqs[i]);

    }

    //read in Permutation
    lastline = fgets( line, MAXLINE, F_inp);
    while( lastline && (line[0] == '%' || line[0] == '>')) {
      lastline = fgets( line, MAXLINE, F_inp);
    }
    token = strtok( line, ",+ ");

    permSize = 0;

    if( !lastline || sscanf( token, "%d", &(perm[ 0])) != 1 || perm[0] > nStrands) {
      printf("Error in %s: %s\nRequesting input manually.\n", inputFile, token);
      fclose( F_inp);
      for( j = 0; j < nStrands; j++) free( seqs[j]);
      free( seqs);
      return 0;
    }

    strcpy( theseq, seqs[ perm[ 0]-1]);
    permSize = 1;
    token = strtok( NULL, ",+ ");
    while( token && sscanf( token, "%d", &(perm[ permSize])) == 1 && perm[ permSize] <= nStrands) {

      strcat( theseq, "+");
      strcat( theseq, seqs[ perm[ permSize]-1] );

      token = strtok( NULL, ",+ ");
      permSize++;
    }

    if( permSize == 0) {
      printf("Error in %s: %s\nRequesting input manually.\n", inputFile, token);
      for( j = 0; j < nStrands; j++) free( seqs[j]);
      free( seqs);
      return 0;
    }
  }

  if( gap != NULL) { //read in energy gap if applicable
    lastline = fgets( line, MAXLINE, F_inp);
    while( lastline && (line[0] == '%' || line[0] == '>') ) {
      lastline = fgets( line, MAXLINE, F_inp);
    }
    token = strtok( line, "%>");

    if( !lastline || sscanf( token, "%f", gap) != 1 || *gap < 0) {
      printf("Error in %s.\nRequesting input manually.\n", inputFile);

      if( Multistranded) {
        for( j = 0; j < nStrands; j++) free( seqs[j]);
        free( seqs);
      }
      return 0;
    }

  }

  if( structure != NULL) {
    lastline = fgets( line, MAXLINE, F_inp);
    while( lastline && (line[0] == '%' || line[0] == '>') ) {
      lastline = fgets( line, MAXLINE, F_inp);
    }
    token = strtok( line, "%");

    if (!lastline) {
      printf("Error in %s.\nRequesting input manually.\n", inputFile);
      printf("Structures must be in dot-paren format for manual input.\n");

      if( Multistranded) {
        for( j = 0; j < nStrands; j++) free( seqs[j]);
        free( seqs);
      }
      return 0;
    }
    else {

      if (isdigit(token[0])) { // Input is pairs
        structure[0] = '\0'; // Flag that the pairs is input

        //  Initialize thepairs
        for (i = 0; i < MAXSEQLENGTH; i++) {
          thepairs[i] = -1;
        }

        token = strtok(token,", \t\n");
        base1 = atoi(token) - 1;
        token = strtok(NULL,", \t\n");
        base2 = atoi(token) - 1;
        thepairs[base1] = base2;
        thepairs[base2] = base1;
        while (fgets(line,MAXLINE,F_inp) != NULL) {
          token = strtok(line,", \t\n");
          base1 = atoi(token) - 1;
          token = strtok(NULL,", \t\n");
          base2 = atoi(token) - 1;
          thepairs[base1] = base2;
          thepairs[base2] = base1;
        }
      }
      else { // input is list of pairs
        thepairs[0] = -1;  // The flag that no pairs information
        if( !sscanf( token, "%s", structure) || strlen( theseq) != strlen( structure) ) {
          printf("Error in %s.\nRequesting input manually.\n", inputFile);

          if( Multistranded) {
            for( j = 0; j < nStrands; j++) free( seqs[j]);
            free( seqs);
          }
          return 0;
        }
      }
    }
  }

  //Calculate v_pi
  if( Multistranded) {
    *v_pi = calculateVPi( perm, permSize);

    for( j = 0; j < nStrands; j++) free( seqs[j]);
    free( seqs);
  }

  fclose( F_inp);

  // Number of unique strands is nStrands
  nUniqueSequences = nStrands;

  return 1;
}

/* ********* */

int calculateVPi( int *perm, int permSize) {
  int k, j;
  int isSymmetric;
  int isSymJ; // is the rotational symmetry permSize/j?
  int v_pi;

  v_pi = 1;
  isSymmetric = 0;
  j = 1;
  while (j <= permSize-1 && !isSymmetric) {
    // j is the distance between symmetric positions
    if (permSize % j == 0) { // only possible if j | permSize
      isSymJ = 1; // Assume symmetric until shown otherwise
      k = 0;
      while (k <= permSize-1 && isSymJ) {
        if (perm[k] != perm[(k+j) % permSize]) {
          isSymJ = 0;
        }
        k++;
      }
      if (isSymJ) {
        v_pi = permSize/j;
        isSymmetric = 1;
      }
      // else check next case
    }
    // else j is not a legal symmetry step
    j++;
  }

  return v_pi;
}

/* ******* */


/* ************ */

void getUserInput( char *theseq, int *v_pi, float *gap, char *structure) {

  char *token;
  char line[ MAXLINE];

  int i, j, permSize = 0;
  char **seqs = NULL;

  int nStrands = 0;

  if( !Multistranded) {
    printf("Enter sequence: ");
    scanf( "%s", theseq);
    *v_pi = 1;
  }
  else {
    printf("Enter number of strands: ");
    scanf( "%d", &nStrands);

    //allocate function variables
    seqs = (char **) malloc( nStrands * sizeof( char*));

    //read in sequences
    for( i = 0; i < nStrands; i++) {
      printf("Enter sequence for strand type %d: \n", i+1);
      scanf( "%s", line);
      seqs[i] = (char*) malloc( (strlen(line)+1)*sizeof( char));
      strcpy( seqs[i], line);

      seqlengthArray[i] = strlen(seqs[i]);

    }


    printf("Enter strand permutation (e.g. 1 2 4 3 2): \n");
    fgets( line, MAXLINE, stdin);
    while( sscanf( line, "%d", &(perm[0])) != 1) {
      fgets( line, MAXLINE, stdin);
    }

    token = strtok( line, ",+ ");

    permSize = 0;
    if( sscanf( token, "%d", &(perm[ 0])) != 1 || perm[0] > nStrands) {
      printf("Illegal permutation index %s\n", token);
      exit(1);
    }

    strcpy( theseq, seqs[ perm[ 0]-1]);
    permSize = 1;
    token = strtok( NULL, ",+ ");
    while( token && sscanf( token, "%d", &(perm[ permSize])) == 1 && perm[ permSize] <= nStrands) {

      strcat( theseq, "+");
      strcat( theseq, seqs[ perm[ permSize]-1] );

      token = strtok( NULL, ",+ ");
      permSize++;
    }
  }

  if( gap != NULL) {
    printf("Enter energy gap (kcal/mole): ");
    if( scanf( "%f", gap) != 1) {
      *gap = -1;
    }
    while( *gap < 0) {
      printf("Reenter energy gap (must be nonnegative): ");
      if( scanf( "%f", gap) != 1) {
        *gap = -1;
      }
    }
  }

  if( structure != NULL) {
    do {
      printf("Enter structure:\n%s\n", theseq);
      scanf("%s", structure);
      printf("%s\n%s\n", theseq, structure);
    } while( strlen( structure) != strlen( theseq) );
  }

  //Calculate v_pi
  if( Multistranded) {
    *v_pi = calculateVPi( perm, permSize);

    for( j = 0; j < nStrands; j++) free( seqs[j]);
    free( seqs);
  }

  // Number of unique strands is nStrands
  nUniqueSequences = nStrands;

}

/* ********* */

void header( int argc, char **argv, char *name, char *outputFile) {

  int i;
  time_t curtime;
  struct tm *loctime;
  FILE *fp;

  curtime = time(NULL); //current time
  loctime = localtime( &curtime);

  if (strcmp(outputFile,"screen") == 0) {
    printf("%s NUPACK %s\n", COMMENT_STRING, NUPACK_VERSION);
    printf("%s Program: %s\n",COMMENT_STRING, name);
    printf("%s Start time: %s", COMMENT_STRING, asctime( loctime));
    printf("%s Command: ", COMMENT_STRING);
    for( i = 0; i < argc; i++) {
      printf( "%s ", argv[ i]);
    }
    printf( "\n");
  }
  else {
    if ((fp = fopen(outputFile,"w")) == 0) {
      printf("Error opening output file %s!\n",outputFile);
      exit(1);
    }
    fprintf(fp,"%s NUPACK %s\n", COMMENT_STRING, NUPACK_VERSION);
    fprintf(fp,"%s Program: %s\n",COMMENT_STRING,name);
    fprintf(fp,"%s Start time: %s", COMMENT_STRING, asctime( loctime));
    fprintf(fp,"%s Command: ", COMMENT_STRING);
    for( i = 0; i < argc; i++) {
      fprintf( fp,"%s ", argv[ i]);
    }
    fprintf( fp, "\n");

    // Include information about the cutoff if this is pairs
    if (CUTOFF > 0.0 && strcmp("pairs",name) == 0) {
      fprintf( fp, "%s Minimum printed pair probability: %g\n",
               COMMENT_STRING,CUTOFF);
    }
    fclose(fp);
  }

}

/* ******* */

void printInputs( int argc, char **argv, const char *seq, int vs,
                  const float * gap, const char *structure, char *outputFile) {

  FILE  *fp = NULL;

  if (strcmp(outputFile,"screen") == 0) {

    if( seq != NULL) {
      printf( "%s Sequence:  %s\n", COMMENT_STRING, seq);
    }

    if( structure != NULL) {
      printf("%s Structure: %s\n", COMMENT_STRING, structure);
    }

    printf( "%s v(pi): %d\n",COMMENT_STRING, vs);

    if( DNARNACOUNT != COUNT) {
      if( DNARNACOUNT == DNA) {
        printf( "%s Parameters: DNA, 1998\n", COMMENT_STRING);
      }
      else if(  DNARNACOUNT == RNA) {
        printf( "%s Parameters: RNA, 1995\n", COMMENT_STRING);
      }
      else if( DNARNACOUNT == RNA37) {
        printf( "%s Parameters: RNA, 1999\n", COMMENT_STRING);
      }
      else if( DNARNACOUNT == USE_SPECIFIED_PARAMETERS_FILE) {
        printf( "%s Parameters: Custom, (%s)\n", COMMENT_STRING, PARAM_FILE);
      }

      printf( "%s Dangles setting: %d\n", COMMENT_STRING, DANGLETYPE);

      printf( "%s Temperature (C): %.1f\n", COMMENT_STRING,
              (float) (TEMP_K - ZERO_C_IN_KELVIN) );

      // Say what salt concentrations were used in the calculation
      printf("%s Sodium concentration: %.4f M\n",COMMENT_STRING,(float) SODIUM_CONC);
      printf("%s Magnesium concentration: %.4f M\n",COMMENT_STRING,(float) MAGNESIUM_CONC);

    }


    if( gap != NULL)
      printf( "%s Energy gap: %.1f\n",COMMENT_STRING, *gap);


    if( DO_PSEUDOKNOTS)
      printf("%s Pseudoknots enabled.\n", COMMENT_STRING);

  }
  else {
    if ((fp = fopen(outputFile,"a")) == 0) {
      printf("Error opening output file %s!\n",outputFile);
      exit(1);
    }

    if(seq != NULL) {
      fprintf(fp, "%s Sequence:  %s\n", COMMENT_STRING, seq);
    }

    if( structure != NULL) {
      fprintf(fp,"%s Structure: %s\n", COMMENT_STRING, structure);
    }

    fprintf(fp, "%s v(pi): %d\n",COMMENT_STRING, vs);

    if( DNARNACOUNT != COUNT) {
      if( DNARNACOUNT == DNA) {
        fprintf(fp, "%s Parameters: DNA, 1998\n", COMMENT_STRING);
      }
      else if(  DNARNACOUNT == RNA) {
        fprintf(fp, "%s Parameters: RNA, 1995\n", COMMENT_STRING);
      }
      else if( DNARNACOUNT == RNA37) {
        fprintf(fp, "%s Parameters: RNA, 1999\n", COMMENT_STRING);
      }
      else if( DNARNACOUNT == USE_SPECIFIED_PARAMETERS_FILE) {
        fprintf(fp, "%s Parameters: Custom, (%s)\n", COMMENT_STRING, PARAM_FILE);
      }

      fprintf(fp, "%s Dangles setting: %d\n", COMMENT_STRING, DANGLETYPE);

      fprintf(fp, "%s Temperature (C): %.1f\n", COMMENT_STRING,
              (float) (TEMP_K - ZERO_C_IN_KELVIN) );

      // Say what salt concentrations were used in the calculation
      fprintf(fp,"%s Sodium concentration: %.4f M\n",COMMENT_STRING, (float) SODIUM_CONC);
      fprintf(fp,"%s Magnesium concentration: %.4f M\n",COMMENT_STRING,(float) MAGNESIUM_CONC);

    }

    if( gap != NULL)
      fprintf(fp, "%s Energy gap: %.1f\n",COMMENT_STRING, *gap);


    if( DO_PSEUDOKNOTS)
      fprintf(fp,"%s Pseudoknots: enabled.\n", COMMENT_STRING);

    fclose(fp);
  }

}
