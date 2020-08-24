
#include "read_command_line.h"

void displayHelp(void) {
  #ifdef DEBUG
    printf("Please read the provided nupack manual instructions.\n");
  #endif

}

void makeUCase(char *line) {
  int i = 0;
  while (line[i]) {
    line[i] = (char)toupper(line[i]);
    i++;
  }
  
}

void print_design_help() {
  printf("Usage: design [OPTIONS] PREFIX\n");
  printf("Design the sequence of one or more interacting strands to adopt a\n");
  printf("target secondary structure at equilibrium\n");
  printf("Example: design -T 25 -material dna1998 -sodium 0.9 -magnesium 0.1 example\n");
  printf("\n");
  PrintNupackThermoHelp();
  printf("\n");
  printf("Initialization:\n");
  printf(" -init initmode   select the initialization mode from AU, CG, RND, SSM\n");
  printf(" -loadinit        load the initial sequence from PREFIX.init\n");
  printf(" -outputinit      output the initial sequence to PREFIX.init\n");
  printf(" -loadseed        load the random seed from PREFIX.seed\n");
  printf(" -outputseed      output the random number generator seed to PREFIX.seed\n");
  printf("\n");
  printf("Objective:\n");
  printf(" -fstop FSTOPVALUE        set the stop condition to FSTOPVALUE\n");
  printf("                          (default 0.01)\n");
  printf(" -prevent PFILE           set the prevented patterns using PFILE\n");
  printf("\n");
  printf("Design process:\n");
  printf(" -mreopt VALUE    set the maximum # of parent node reoptimizations\n");
  printf("                  to VALUE\n");
  printf(" -mleafopt VALUE  set the maximum # of leaf reoptimizations no VALUE\n");
  printf("\n");
  printf("Output control:\n");
  printf(" -pairs           save the pair probabilities to PREFIX.ppairs\n");
  printf(" -cutoff CUTOFF   only save pair probabilities above CUTOFF in\n");
  printf("                  the pair probabilities file (default 0.001)\n");
  printf(" -json            generate a JSON formatted summary of the design\n");
  exit(0);
}


/* ****************************************************************** */
void readCommandLine(int nargs, char **args, char* inputPrefix, char *psFile, int *initMode, 
      int *bypassDesign, int *bypassHierarchy, int *bypassGuidance, int *designMode, 
      int *loadSeeds, int *mReopt,int *mLeafopt, DBL_TYPE *nRatio,int *loadInit, int *hMin, DBL_TYPE *pcut, int *quick,int * output_init,int * output_seed, int * output_ppairs, int * output_json) {//, char *inputFile) {
  
  int options;  // Counters used in getting flags
  int showHelp=0; // ShowHelp = 1 if help option flag is selected
  int option_index;
  char line[MAXLINE];
  char param[MAXLINE];
  char prevent[MAXLINE];
  char initMethod[MAXLINE];
  char programMode[MAXLINE];
  float temp;

  #ifdef DEBUG
  #endif
  
  static struct option long_options [] = 
  {
    
    {"T", required_argument,          NULL, 'a'},
    {"dangles", required_argument,    NULL, 'b'},
    {"material", required_argument,     NULL, 'c'},
    {"prevent", required_argument,    NULL, 'd'},
    {"mleafopt",required_argument, NULL,'e'},
    {"init", required_argument, NULL, 'f'},
    {"nodesign", no_argument, NULL, 'g'},
    {"help",no_argument,NULL,'h'},
    {"noguidance", no_argument, NULL, 'i'},
    {"designmode", required_argument, NULL, 'j'},
    {"loadseed", no_argument, NULL, 'k'},
    {"mreopt", required_argument, NULL, 'l'},
    {"fstop", required_argument, NULL, 'm'}, 
    {"loadinit", no_argument, NULL, 'n'},
    {"cutoff", required_argument, NULL, 'o'},
    {"quick", no_argument, NULL,'p'},
    {"sodium",required_argument,NULL,'q'},
    {"magnesium",required_argument,NULL,'r'},
    // output options
    {"outputinit",no_argument,NULL,'s'},
    {"outputseed",no_argument,NULL,'t'},
    {"nohierarchy", no_argument, NULL, 'u'},
    {"pairs",no_argument,NULL,'v'},
    {"json",no_argument,NULL,'w'},
    {0, 0, 0, 0}
  };
  
  //initialize global parameters
  SODIUM_CONC = 1.0;
  MAGNESIUM_CONC = 0.0;
  //USE_LONG_SALT = 0;
  DNARNACOUNT = RNA;
  DANGLETYPE = 1;
  TEMP_K = 37.0 + ZERO_C_IN_KELVIN;
  *initMode = RANDOM_INIT;
  *bypassDesign = 0;
  *bypassHierarchy = 0;
  *bypassGuidance = 0;
  *designMode = N_OPTIMIZATION;
  psFile[0] = 0;
  *loadSeeds = 0;
  *loadInit = 0;
  *mLeafopt = M_LEAFOPT_DEFAULT;
  *mReopt = M_REOPT_DEFAULT;
  *nRatio = N_RATIO_DEFAULT;
  *hMin = H_MIN_RNA;
  *pcut = 0.001;
  *quick = 0;
  *output_init = 0;
  *output_seed = 0;
  *output_ppairs = 0;
  *output_json = 0;
  
  // Get the option flags
  while (1) {
    /* getopt_long stores the option index here. */
    option_index = 0;
    options = getopt_long_only (nargs, args, 
      "a:b:c:d:e:f:ghij:kl:m:no:pq:r:stuvw", long_options, 
      &option_index);
    
    // Detect the end of the options.
    if (options == -1)
      break;
    
      switch (options) {
        
        
      // T
      case 'a':
        strcpy( line, optarg);
        if( sscanf(line, "%f", &(temp)) != 1) {
          printf("Invalid T value\n");
          exit(1);
        }
        TEMP_K = (DBL_TYPE) (temp + ZERO_C_IN_KELVIN);
        break;
        
      // dangles
      case 'b':
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
        
      // material
      case 'c':
        strcpy( line, optarg);
        if( sscanf(line, "%s", PARAM_FILE) != 1) {
          printf("Invalid parameters value\n");
          exit(1);
        }
        if( strcmp( PARAM_FILE, "rna") == 0 || strcmp( PARAM_FILE, "rna1995") == 0) {
          DNARNACOUNT = RNA;
          *hMin = H_MIN_RNA;
        }
        else if( strcmp( PARAM_FILE, "dna") == 0 || strcmp( PARAM_FILE, "dna1998") == 0) {
          DNARNACOUNT = DNA;
          *hMin = H_MIN_DNA;
        }
        else if( strcmp( PARAM_FILE, "rna37") == 0 || strcmp( PARAM_FILE, "rna1999") == 0) {
          if(strcmp(PARAM_FILE,"rna37") == 0) {
            printf("Parameter specification using rna37 has been deprecated. Please use rna1999 instead.\n");
          }
          DNARNACOUNT = RNA37;
          *hMin = H_MIN_RNA;
        }
        else {
          // default assumes RNA
          DNARNACOUNT = USE_SPECIFIED_PARAMETERS_FILE;
          *hMin = H_MIN_RNA;
        }
        break;
      // prevent
      case 'd':
        strcpy( line, optarg);
        if( sscanf(line, "%s", prevent) != 1) {
          printf("Invalid prevented strings value\n");
          exit(1);
        }
        strcpy(psFile,prevent);
        break;
      // mleafopt
      case 'e':
        strcpy( line, optarg);
        if( sscanf(line, "%d", mLeafopt) != 1) {
          printf("Invalid mreopt value\n");
          exit(1);
        }
        break;
      // init
      case 'f':
        strcpy( line, optarg);
        if( sscanf(line, "%s", initMethod) != 1) {
          printf("Invalid initialization method\n");
          exit(1);
        }
        makeUCase(initMethod);
        
  
        if (!strcmp(initMethod, "CG") 
          || !strcmp(initMethod, "GC")) {
          *initMode = CG_INIT;
        }
        else if (!strcmp(initMethod, "AU") 
          || !strcmp(initMethod, "UA") 
          || !strcmp(initMethod, "TA") 
          || !strcmp(initMethod, "AT")) {
          *initMode = AU_INIT;
        
        }
        else if (!strcmp(initMethod, "RND") 
          || !strcmp(initMethod, "RANDOM")) {
          *initMode = RANDOM_INIT;        
        }
        else if (!strcmp(initMethod, "SSM")) {
          *initMode = SSM_INIT;
        }
        else {
          printf("Invalid initialization method: %s\n", initMethod);
          exit(1);
        }
        
        break;
      // nodesign
      case 'g':
        *bypassDesign = 1;
        break;
      case 'h':
        print_design_help();
        break;
      // noguidance
      case 'i':
        *bypassGuidance = 1;
        break;
      // designmode
      case 'j':
        strcpy( line, optarg);
        if( sscanf(line, "%s", programMode) != 1) {
          printf("Invalid design mode\n");
          exit(1);
        }
        makeUCase(programMode);
        
  
        if (!strcmp(programMode, "MFE")) {
          // printf("MFE optimization selected\n");
          *designMode = MFE_OPTIMIZATION;
        }
        else if (!strcmp(programMode, "P")) {
          // printf("P optimization selected\n");
          *designMode = P_OPTIMIZATION;
        } else {
          // printf("N optimization selected\n");
        }
        break;
      // loadseed
      case 'k':
        *loadSeeds = 1;
        break;
      // mreopt
      case 'l':
        strcpy( line, optarg);
        if( sscanf(line, "%d", mReopt) != 1) {
          printf("Invalid mreopt value\n");
          exit(1);
        }
        break;
      // fstop
      case 'm':
        strcpy( line, optarg);
        if( sscanf(line, "%Lf", nRatio) != 1) {
          printf("Invalid fstop value\n");
          exit(1);
        }
        break;    
      // loadinit
      case 'n':
        *loadInit = 1;
        break;    
      // cutoff
      case 'o':
        strcpy( line, optarg);
        if( sscanf(line, "%Lf", pcut) != 1) {
          printf("Invalid probability cutoff value\n");
          exit(1);
        }
        // *output_ppairs = 1;
        break;
      // quick
      case 'p':
        *quick = 1;
        break;
      // sodium
      case 'q':
        strcpy(line,optarg);
        if( sscanf(line, "%Lf", &SODIUM_CONC) != 1) {
          printf("Invalid sodium concentration\n");
          exit(1);
        }
        break;
      // magnesium
      case 'r':
        strcpy(line,optarg);
        if( sscanf(line, "%Lf", &MAGNESIUM_CONC) != 1) {
          printf("Invalid magnesium concentration\n");
          exit(1);
        }
        break;
      // outputinit
      case 's':
        *output_init = 1;
        break;
      // outputseed
      case 't':
        *output_seed = 1;
        break;
      // nohierarchy
      case 'u':
        *bypassHierarchy = 1;
        break;
      // ppairs
      case 'v':
        *output_ppairs = 1;
        break;
      // json
      case 'w':
        *output_json = 1;
        break;

      default:
        abort ();
      }
  }
  
  if (*bypassDesign) {
    *mReopt = -1;
  }
  
  if (optind == nargs) { // There's no input from the user
    printf("You must have a prefix or an input file on the command line.\n");
    printf("For instructions on running this program, run it with the ");
    printf("-help flag.\n\nExiting....\n\n");
    // exit(ERR_NOINPUT);
    exit(1);
  }
  else {
    strcpy(inputPrefix,args[optind]);
  }
  
  /*
  if( nargs >= 2) {
    strcpy( inputFile, args[ nargs-1]);
    strcat(inputFile,".in");
  }
  */
  
  
  if (showHelp) {
    displayHelp();
    exit(1);
  }
} 

    
    
    
