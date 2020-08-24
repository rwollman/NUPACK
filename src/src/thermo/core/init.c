#include "init.h"

/* ************************************************* */
void ReadSequence( int *seqlength, char **seq, char filename[ MAXLINE] ) {
  FILE *fp;
  char line[MAXLINE];
  int i; //position in line
  int linenumber; //line number
  char tempSeq[ MAXSEQLENGTH];

  *seqlength = 0;

#ifndef PRINTRESULTSONLY  
  printf("Reading Input File...\n");
#endif
  fp = fopen(filename, "r");

  if( fp == NULL) {  // Make sure input file exits 
    fprintf(stderr, "Error opening file %s\n", filename);
    exit(1);  
  }

  linenumber = 0;
  while( fgets(line, MAXLINE, fp)!= NULL ) {  // Read lines
    linenumber++;
    i = 0;

    while( line[0] != '>' && line[i] != '\n') {
      line[i] = (char)toupper( line[i] );
      
      if( line[i] != 'A' && line[i] != 'T' && line[i] != 'C'
         && line[i] != 'G' && line[i] != 'U') {
           if( line[i] != ' ' && line[i] != '\t') {
             fprintf(stderr, "Invalid base at line %d, position %d\n", linenumber, i+1);
             fclose( fp);
             exit(1);
           }
         }
      else {
        if( line[i] != ' ' && line[i] != '\t' && line[i] != '\n') {
          tempSeq[ *seqlength] = line[i];
          // seq[ *seqlength] = line[i];
          (*seqlength)++;
        }
      }
      i++;
    }
  }

  *seq = (char *) calloc( *seqlength + 1, sizeof( char) );
  if( *seq == NULL) {
    fprintf(stderr, "ReadSequence: unable to allocate %lu bytes for '*seq'\n",
        (unsigned long)*seqlength+1 * sizeof( char));
    exit(1);
  }

  for( i = 0; i < *seqlength; i++) {
    (*seq)[ i] = tempSeq[ i];
  }
  (*seq)[*seqlength] = '\0';

  fclose( fp);
}

/* **************************** */
int getSequenceLength( char *seq, int *nStrands /*, seq2, nicks*/) {
  
  int i; //position in sequence
  int done = FALSE;
  int seqlength = 0;
  char tmpC;

  i = 0;
  *nStrands = 1;
  while( done == FALSE) {
    tmpC = (char)toupper( seq[i]);
    if( tmpC == '+') {
      (*nStrands)++;
      //nicks[ seqlength] = 1;
    }
    else if( tmpC != 'A' && tmpC != 'T' && tmpC != 'C'
            && tmpC != 'G' && tmpC != 'U') {
              done = TRUE; 
            }
    else {
      //seq2[ seqlength] = tmpC;
      seqlength++;
    }
    i++;
  }

  if( seqlength > MAXSEQLENGTH) {
    fprintf(stderr, "Sequences longer than maximum of %d\n", MAXSEQLENGTH);
    exit(1);
  }

  return seqlength;
}

/* **************************** */
int getSequenceLengthInt( int seq[], int *nStrands /*, seq2, nicks*/) {

  int i; //position in sequence
  int done = FALSE;
  int seqlength = 0;
  int tmpC;

  i = 0;
  *nStrands = 1;
  while( done == FALSE) {
    tmpC = seq[i];
    if( tmpC == STRAND_PLUS) {
      (*nStrands)++;
      //nicks[ seqlength] = 1;
    }
    else if( tmpC != BASE_A && tmpC != BASE_T && tmpC != BASE_C
            && tmpC != BASE_G && tmpC != BASE_U) {
              done = TRUE; 
            }
    else {
      //seq2[ seqlength] = tmpC;
      seqlength++;
    }
    i++;
  }

  if( seqlength > MAXSEQLENGTH) {
    fprintf(stderr, "Sequences longer than maximum of %d\n", MAXSEQLENGTH);
    exit(1);
  }

  return seqlength;
}

/* *********** */
void processMultiSequence( int inputSeq[], int seqlength, int nStrands,
                           int seq[], int nicks[]) {
  
  int i, j;
  int nNick = 0;
  int done;
  int *tmpSeq;


  tmpSeq = (int*) malloc( (seqlength+nStrands)*sizeof( int) );
  memcpy( tmpSeq, inputSeq, (seqlength+nStrands)*sizeof( int));

  j = 0;
  for( i = 0; i < seqlength; i++) {
    done = FALSE;
    while( !done) {
      if( tmpSeq[j] == BASE_A || tmpSeq[j] == BASE_C || tmpSeq[j] == BASE_G ||
        tmpSeq[j] == BASE_T || tmpSeq[j] == BASE_U) {
        done = TRUE;
        seq[i] = tmpSeq[j];
      }
      else if( tmpSeq[j] == STRAND_PLUS) {
        nicks[nNick++] = i-1;
      }
      j++;

      if( j >= seqlength + nStrands && !done) {
        fprintf(stderr, "\nError in processing sequence:\n%d\n", inputSeq[0]);
        fprintf(stderr, "seqlength = %d, nStrands = %d\n", seqlength, nStrands);
        exit(1);
      }
    }
  }
  seq[ seqlength] = -1;
  free(tmpSeq);
}


/* ********************************************* */
void InitLDoublesMatrix( DBL_TYPE **Q, int size, char name[]) {
  // Allocate cleared memory for a DBL_TYPEs matrix.
  *Q =  (DBL_TYPE *) calloc( size, sizeof( DBL_TYPE));
  if( *Q == NULL) {
    fprintf(stderr, "InitLDoublesMatrix: unable to allocate %lu bytes for"
       " %s!\n", (unsigned long) size * sizeof( DBL_TYPE),  name);
    exit(1);
  }
}

void ClearLDoublesMatrix(DBL_TYPE **Q, int size, char name[]) {
  (void)name;
  memset(*Q, 0,size * sizeof(DBL_TYPE));
}

/* ******************************************** */
void nonZeroInit( DBL_TYPE Q[], int seq[], int seqlength) {
  // Set Q[i, i-1] = 1.
  int i;

  for( i = 0; i <= seqlength; i++) {
    Q[ pf_index(i, i-1, seqlength)] = ExplDangle(i,i-1,seq,seqlength);
  }
}

/* ***************************** */
void manageQgIx( DBL_TYPE **QgIx, DBL_TYPE **QgIx_1, DBL_TYPE **QgIx_2, int d, int seqlength) {
  // Allocate and deallocate QgIx matrices

  long int maxStorage;
  long int dTest;
  DBL_TYPE *temp;

  if( d == 16) { //First time
    maxStorage = 0;
    for( dTest = d; dTest < d+3; dTest++) {
      maxStorage = MAX( maxStorage, (seqlength-dTest)*(dTest-5)*
                       ((dTest-1)*(dTest-2)/2) );
    }

    *QgIx = (DBL_TYPE *) calloc( maxStorage, 
                                sizeof( DBL_TYPE));
    *QgIx_1 = (DBL_TYPE *) calloc( maxStorage, 
                                  sizeof( DBL_TYPE));
    *QgIx_2 = (DBL_TYPE *) calloc( maxStorage, 
                                  sizeof( DBL_TYPE));

    if( *QgIx == NULL || *QgIx_1 == NULL || *QgIx_2 == NULL) {
      fprintf(stderr, "manageQgIx: unable to allocate memory for QgIx, QgIx_1 or QgIx_2\n");
      exit(1);
    }
  }
  else if( d > 16) { // every case beyond the first
    maxStorage = 0;
    for( dTest = d; dTest < d+3; dTest++) {
      maxStorage = MAX( maxStorage, (seqlength-dTest)*(dTest-5)*
                       ((dTest-1)*(dTest-2)/2) );
    }

    temp = *QgIx;
    *QgIx = *QgIx_1;
    *QgIx_1 = *QgIx_2;
    free(temp);
    temp = NULL;

    *QgIx_2 = (DBL_TYPE *) calloc( maxStorage, sizeof(DBL_TYPE) );
    if( *QgIx_2 == NULL) {
      fprintf(stderr, "Error in QgIx_2 allocation! %ld %d %d\n", 
             maxStorage, d, seqlength);
      exit(1);
    }
  }
}

/* **************************************** */
void manageFgIx( DBL_TYPE **FgIx, DBL_TYPE **FgIx_1, DBL_TYPE **FgIx_2, int d, int seqlength) {
  // Allocate and deallocate FgIx matrices

  int i;
  long int maxStorage;
  long int dTest;
  DBL_TYPE *temp;

  if( d > 16) { // every case beyond the first
    maxStorage = 0;
    for( dTest = d; dTest < d+3; dTest++) {
      maxStorage = MAX( maxStorage, (seqlength-dTest)*(dTest-5)*
                       ((dTest-1)*(dTest-2)/2) );
    }

    temp = *FgIx;
    *FgIx = *FgIx_1;
    *FgIx_1 = *FgIx_2;
    free(temp);
    temp = NULL;

    *FgIx_2 = (DBL_TYPE *) calloc( maxStorage, sizeof(DBL_TYPE) );
    if( *FgIx_2 == NULL) {
      fprintf(stderr, "Error in QgIx_2 allocation! %ld %d %d\n", 
             maxStorage, d, seqlength);
      exit(1);
    }
    for( i = 0; i < maxStorage; i++) {
      (*FgIx_2)[i] = NAD_INFINITY;
    }
  }
  else if( d == 16) { //First time

    maxStorage = 0;
    for( dTest = d; dTest < d+3; dTest++) {
      maxStorage = MAX( maxStorage, (seqlength-dTest)*(dTest-5)*
                       ((dTest-1)*(dTest-2)/2) );
    }

    *FgIx = (DBL_TYPE *) calloc( maxStorage, 
                                sizeof( DBL_TYPE));
    *FgIx_1 = (DBL_TYPE *) calloc( maxStorage, 
                                  sizeof( DBL_TYPE));
    *FgIx_2 = (DBL_TYPE *) calloc( maxStorage, 
                                  sizeof( DBL_TYPE));

    if( *FgIx == NULL || *FgIx_1 == NULL || *FgIx_2 == NULL) {
      fprintf(stderr, "Error in FgIx, FgIx_1, FgIx_2 allocation\n");
      exit(1);
    }

    for( i = 0; i < maxStorage; i++) {
      (*FgIx_2)[i] = (*FgIx_1)[i] = (*FgIx)[i] = NAD_INFINITY;
    }
  }
}

/* ******************** */
void PrecomputeValuesN5( int seqlength) {
  int i;

  PrecomputeValuesN5f( seqlength);
  for( i = 8; i < seqlength; i++) {
    sizeTerm[i] = EXP_FUNC( -sizeTerm[i]/(kB*TEMP_K) );    
  }
}

/* *********************** */
void PrecomputeValuesN5f( int seqlength) {
  int i;

  sizeTerm = (DBL_TYPE *) calloc( seqlength, sizeof( DBL_TYPE) );
  if( sizeTerm == NULL) {
    fprintf(stderr, "UnaPrecomputeValuesN5f: ble to allocate memory for sizeTerm!\n");
    exit(1);
  }

  for( i = 8; i < seqlength; i++) {
    if( i <= 30) {
      sizeTerm[i] = -loop37[ i - 1];
    }
    else {
      sizeTerm[i] = -loop37[ 30 - 1];
      sizeTerm[i] -= sizeLog (i); //1.75*kB*TEMP_K*logl( i/30.0); 
    }

    if( i + 2 <= 30) {
      sizeTerm[i] += loop37[ i+2 - 1];
    }
    else {
      sizeTerm[i] += loop37[ 30 - 1];
      sizeTerm[i] += sizeLog (i+2); //1.75*kB*TEMP_K*LOG_FUNC( (i+2)/30.0); 
    } 
  }
}

/* *************************************************** */
void manageQx( DBL_TYPE **Qx, DBL_TYPE **Qx_1, DBL_TYPE **Qx_2, int len, int seqlength) {
  // Allocate and deallocate QbIx matrices

  int i;
  int maxStorage;
  DBL_TYPE *temp;

  maxStorage = (seqlength - len)*(len - 1);
  for( i = len+1; i <= len + 2; i++) {
    maxStorage = MAX( maxStorage, ( seqlength - i)*(i - 1) );
  }

  if( len == 11) { //first use of these matrices
    *Qx = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Qx_1 = 
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Qx_2 = 
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    
    if( *Qx == NULL || *Qx_1 == NULL || *Qx_2 == NULL) {
      fprintf(stderr, "Error in Qx, Qx_1, Qx_2 allocation\n");
    }
    
  }

  //else if( len == seqlength - 1) {
  //new
  //free( *Qx); *Qx = NULL;
  //free( *Qx_1); *Qx_1 = NULL;
  //free( *Qx_2); *Qx_2 = NULL;
  
  /* old
  temp = *Qx;
  *Qx = *Qx_1;
  *Qx_1 = temp; //this is needed for pairs calculations
  */
  
  /*
  *Qx_2 = //this array is worthless
  (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
  */
  //}

  else if( len > 11) {

    temp = *Qx;
    *Qx = *Qx_1;
    *Qx_1 = *Qx_2;
    free( temp);
    *Qx_2 = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));

    if( *Qx_2 == NULL) {
      fprintf(stderr, "Error in Qx_2 allocation\n");
    }
  }
}

/* *************************************************** */
void manageFx( DBL_TYPE **Fx, DBL_TYPE **Fx_1, DBL_TYPE **Fx_2, int len, int seqlength) {
  // Allocate and deallocate Fx matrices
  
  int i;
  int maxStorage;
  DBL_TYPE *temp;
  
  maxStorage = (seqlength - len)*(len - 1);
  for( i = len+1; i <= len + 2; i++) {
    maxStorage = MAX( maxStorage, ( seqlength - i)*(i - 1) );
  }
  
  
  if( len == 11) { //first use of these matrices
    *Fx = (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Fx_1 = 
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    *Fx_2 = 
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    
    if( *Fx == NULL || *Fx_1 == NULL || *Fx_2 == NULL) {
      fprintf(stderr, "Error in Fx, Fx_1, Fx_2 allocation\n");
    }
    
    for( i = 0; i < maxStorage; i++) 
      (*Fx)[i] = (*Fx_1)[i] = (*Fx_2)[i] = NAD_INFINITY;
    
  }
  
  else if( len > 11) {
    
    temp = *Fx;
    *Fx = *Fx_1;
    *Fx_1 = *Fx_2;
    free( temp);
    *Fx_2 =
      (DBL_TYPE *) calloc( maxStorage, sizeof( DBL_TYPE));
    
    if( *Fx_2 == NULL) {
      fprintf(stderr, "Error in Fx_2 allocation\n");
    }
    for( i = 0; i < maxStorage; i++) 
      (*Fx_2)[i] = NAD_INFINITY;
    
  }
}


/* ************************************** */
DBL_TYPE computeSaltCorrection(DBL_TYPE sodiumConc, DBL_TYPE magnesiumConc,
			       int useLongHelix) {

  // No correction for RNA since we don't have parameters
  if (DNARNACOUNT != DNA || (sodiumConc == 1.0 && magnesiumConc == 0.0)) { 
    return 0.0;
  }

  // Ignore magnesium for long helix mode (not cited why, for consistency with Mfold)
  if (useLongHelix) { 
    return -(0.2 + 0.175*log(sodiumConc)) * TEMP_K / 310.15;
  }

  return -0.114*log(sodiumConc + 3.3*sqrt(magnesiumConc)) * TEMP_K / 310.15;
}



/* ************************************** */
void LoadEnergies( void) {
  
  const char *default_param_files[] = { "dna1998", "rna1995", "rna1999"};
  
  DBL_TYPE H_loop37[90];  
  DBL_TYPE H_tloop_energy[4096];
  DBL_TYPE H_triloop_energy[2048];
  DBL_TYPE H_MMEnergiesHP[6*16];
  DBL_TYPE H_MMEnergiesIL[256];
  DBL_TYPE H_IL_SInt2[16*36];
  DBL_TYPE H_IL_SInt4[256*36];
  DBL_TYPE H_IL_AsInt1x2[64*36];
  DBL_TYPE H_dangle_energy[48];
  DBL_TYPE H_asymmetry_penalty[4];
  DBL_TYPE H_max_asymmetry;
  DBL_TYPE H_Stack[36];
  DBL_TYPE H_ALPHA_1, H_ALPHA_2, H_ALPHA_3, H_BETA_1, H_BETA_2, 
  H_BETA_3, H_BETA_1M, H_BETA_1P;
  DBL_TYPE H_POLYC3, H_POLYCINT, H_POLYCSLOPE;
  DBL_TYPE H_AT_PENALTY;
  DBL_TYPE H_BIMOLECULAR;
  
  DBL_TYPE G_loop37[90];  
  DBL_TYPE G_tloop_energy[4096];
  DBL_TYPE G_triloop_energy[2048];
  DBL_TYPE G_MMEnergiesHP[6*16];
  DBL_TYPE G_MMEnergiesIL[256];
  DBL_TYPE G_IL_SInt2[16*36];
  DBL_TYPE G_IL_SInt4[256*36];
  DBL_TYPE G_IL_AsInt1x2[64*36];
  DBL_TYPE G_dangle_energy[48];
  DBL_TYPE G_asymmetry_penalty[4];
  DBL_TYPE G_max_asymmetry;
  DBL_TYPE G_Stack[36];
  DBL_TYPE G_ALPHA_1, G_ALPHA_2, G_ALPHA_3, G_BETA_1, G_BETA_2,
  G_BETA_3, G_BETA_1M, G_BETA_1P;
  DBL_TYPE G_POLYC3, G_POLYCINT, G_POLYCSLOPE;
  DBL_TYPE G_AT_PENALTY;
  DBL_TYPE G_BIMOLECULAR;
  DBL_TYPE water_conc;

  //temporary storage of data
  int nRead;
  int array[MAXLINE];
  char *token;
  char tetraloop[6];
  char triloop[5];
  int indexL, tmpIndex, index4;
  
  FILE *fp;	
  char line[MAXLINE];
  int i, j, k;
  char fileG[MAX_FILENAME_LEN];
  char fileH[MAX_FILENAME_LEN];
  const char *cur_loc = NULL;
  char fileNameRoot[MAX_FILENAME_LEN];
  
  static DBL_TYPE temp = 0;
  static int energySet = FALSE;
  static int params = -1;
  static int dtype = -1;
  static char parameterFileName[MAX_FILENAME_LEN] = "";
  int n_param_locations = 2;
  int i_param_location;
  const char *default_param_locations[] = {
    "/usr/local/share/nupack",
    "/usr/share/nupack"
  };
  int n_env_prefixes = 2;
  int n_env_suffixes = 2;
  const char *default_env_var_prefixes[] = {
    "NUPACKINSTALL",
    "NUPACKHOME"
  };
  const char *default_env_var_suffixes[] = {
    "share/nupack",
    "parameters",
  };
  int i_env_prefix;
  int i_env_suff;

  /*
  //check if invalid temperature and parameters are used
  if( DNARNACOUNT == RNA37 && TEMP_K != 310.15) {
  printf("The rna1999 RNA parameters can only be used at 37 C.  Please try again!");
  exit(1);
  }
  */
  
  //check if parameters have been loaded (or have changed)
  if( temp != TEMP_K || temp == 0 || params != DNARNACOUNT || 
     (params == USE_SPECIFIED_PARAMETERS_FILE && strcmp( parameterFileName, PARAM_FILE) != 0) 
     || dtype != DANGLETYPE ) {
    energySet = FALSE;
    temp = TEMP_K;
    params = DNARNACOUNT;
    dtype = DANGLETYPE;
  }
  
  if( energySet == TRUE) return; //only load energy once
  energySet = TRUE;
  
  if( DNARNACOUNT == COUNT) {
    setParametersToZero();
    return;
  }
  
  // Get density of water
  water_conc = (DBL_TYPE) WaterDensity(TEMP_K - ZERO_C_IN_KELVIN);

  // Compute the salt correction
  SALT_CORRECTION = computeSaltCorrection(SODIUM_CONC,MAGNESIUM_CONC,
					  USE_LONG_HELIX_FOR_SALT_CORRECTION);
  
  
  // Parameter file input.  If a path was given as a command line parameter,
  // follow this path and return an error if the parameter files cannot be
  // found.
  //
  // Otherwise, check the following locations for the parameter files in the
  // following:
  //  1) Try the current directory
  //  2) Get the directory of the program name from argv[0]
  //      - try BINDIR/../share/nupack/
  //      - try BINDIR/../parameters/
  //  3) Try all BINDIR directories in the PATH
  //      - try BINDIR/../share/nupack/
  //      - try BINDIR/../parameters/
  
  if( DNARNACOUNT < USE_SPECIFIED_PARAMETERS_FILE) {
    strcpy( fileNameRoot, default_param_files[ DNARNACOUNT]);
    strcpy( parameterFileName, ""); //Used to check if parameter files need to be reloaded.
  }
  else if( DNARNACOUNT == USE_SPECIFIED_PARAMETERS_FILE) {
    strcpy( fileNameRoot, PARAM_FILE);
    strcpy( parameterFileName, PARAM_FILE); //store this to check if parameter reload is needed.
  }
  
  //check first for .dG parameter file using current directory as home
  strcpy( fileG, fileNameRoot);
  strcat( fileG, ".dG");

  strcpy( fileH, fileNameRoot);
  strcat( fileH, ".dH");
  
  if( !fileExists( fileG) ) {
    //if files not found, use environment variable NUPACKINSTALL as root
    
    i_env_prefix = 0;
    i_env_suff = 0;
    while (!fileExists(fileG) && i_env_prefix < n_env_prefixes) {
      cur_loc = getenv(default_env_var_prefixes[i_env_prefix]);
      if (cur_loc != NULL) {
        for (i_env_suff = 0; i_env_suff < n_env_suffixes; i_env_suff++) {
          snprintf(fileG, MAX_FILENAME_LEN, "%s/%s/%s.dG", cur_loc, 
              default_env_var_suffixes[i_env_suff], fileNameRoot);
          snprintf(fileH, MAX_FILENAME_LEN, "%s/%s/%s.dH", cur_loc,
              default_env_var_suffixes[i_env_suff], fileNameRoot);
          if (fileExists(fileG)) {
            break;
          }
        }
      }
      i_env_prefix ++;
    }

    if (!fileExists(fileG)) {
      for (i_env_suff = 0; i_env_suff < n_env_suffixes; i_env_suff++) {
        snprintf(fileG, MAX_FILENAME_LEN, "%s/../%s/%s.dG", COMMAND_PATH,
            default_env_var_suffixes[i_env_suff], fileNameRoot);
        snprintf(fileH, MAX_FILENAME_LEN, "%s/../%s/%s.dH", COMMAND_PATH,
            default_env_var_suffixes[i_env_suff], fileNameRoot);
        if (fileExists(fileG)) {
          break;
        }
      }
    }

    i_param_location = 0;
    while (!fileExists(fileG) && i_param_location < n_param_locations) {
      cur_loc = default_param_locations[i_param_location];
      snprintf(fileG, MAX_FILENAME_LEN, "%s/%s.dG", cur_loc,
          fileNameRoot);
      snprintf(fileH, MAX_FILENAME_LEN, "%s/%s.dH", cur_loc,
          fileNameRoot);
      i_param_location ++;
    }

    if (!fileExists(fileG)) {
      fprintf(stderr, "Unable to find %s.dG in any lookup location.\n", fileNameRoot);
      fprintf(stderr, "Ensure that parameter files are in the current directory,\n");
      fprintf(stderr, "$NUPACKINSTALL/parameters, or another default lookup directory\n");
      exit(1);
    } 
  }

  fp = fopen( fileG, "r");
  
  if( fp == NULL) {  // Make sure input file exits 
    fprintf(stderr, "Error opening loop data file: %s\n", fileG);
    exit(1);  
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Stacking data
  for( i = 0; i < 6; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    int tmp_array;
    while( token != NULL) {
      
      //brutal hack to fix mac platform upgrade errors
      strerror(0);
      strerror(1);

      tmp_array=0;
      if( sscanf( token, "%d", &tmp_array ) == 1) {
        array[nRead]=tmp_array;
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in stacking data format\n");
      exit(1);
    }
    for( j = 0; j < 6; j++) {
      Stack[i*6+j] = G_Stack[i*6+j] = (DBL_TYPE) array[j]/100.0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  for( i = 0; i < 3; i++) {
    while( line[0] == '>') {
      fgets( line, MAXLINE, fp);
    }
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead > 30) {
      fprintf(stderr, "Error in Loop energies data\n");
      exit(1);
    }
    for( j = 0; j < 30; j++) {
      if( nRead-1 >= j) {
        loop37[30*(2-i)+j] = G_loop37[30*(2-i)+j] = (DBL_TYPE) array[j]/100.0;
      }
      else {
        loop37[30*(2-i)+j] = G_loop37[30*(2-i)+j] = G_loop37[30*(2-i)+nRead-1]+
          1.75*kB*TEMP_K*LOG_FUNC( (j+1)/(1.0*nRead));
      }
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "Error in asymmetry terms!\n");
    exit(1);
  }
  
  for( j = 0; j < 4; j++) {
    asymmetry_penalty[j] = G_asymmetry_penalty[j] = (DBL_TYPE) array[j]/100.0;
  }
  max_asymmetry = G_max_asymmetry = (DBL_TYPE) array[4]/100.0;
  
  //Triloops
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 2048; i++) {
    triloop_energy[i] = G_triloop_energy[i] = 0;
  }
  while( line[0] != '>') {
    
    if( sscanf( line, "%s %d", triloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 5; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( triloop[4-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( triloop[4-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( triloop[4-i] == 'U' || triloop[4-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( triloop[4-i] != 'A') {
          fprintf(stderr, "Error in triloop indexing %s\n", triloop);
        }
      }
      triloop_energy[ indexL] = G_triloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
    }
    else {
      fprintf(stderr, "Error in triloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  //Tetraloops
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 4096; i++) {
    tloop_energy[i] = G_tloop_energy[i] = 0;
  }
  while( line[0] != '>') {
    if( sscanf( line, "%s %d", tetraloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 6; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( tetraloop[5-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( tetraloop[5-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( tetraloop[5-i] == 'U' || tetraloop[5-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( tetraloop[5-i] != 'A') {
          fprintf(stderr, "Error in tetraloop indexing %s\n", tetraloop);
        }
      }
      tloop_energy[ indexL] = G_tloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
    }
    else {
      fprintf(stderr, "Error in tetraloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Hairpin: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch hairpin format! %d\n", nRead);
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      MMEnergiesHP[ 6*i + j] = G_MMEnergiesHP[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Interior: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch Interior format!\n");
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      MMEnergiesIL[ 6*i + j] = G_MMEnergiesIL[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "1. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      dangle_energy[i*4+j]  = G_dangle_energy[i*4+j] = (DBL_TYPE) array[j]/100.0;
      
      if( DANGLETYPE == 0)  dangle_energy[i*4+j] = G_dangle_energy[i*4+j] = 0.0;

    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "2. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      dangle_energy[24+ i*4+j] = G_dangle_energy[24+ i*4+j] = (DBL_TYPE) array[j]/100.0;
      if( DANGLETYPE == 0) dangle_energy[24+ i*4+j] = G_dangle_energy[24+ i*4+j] = 0;
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }


  //Multiloop parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) ) ==1) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  
  if( nRead != 3) {
    fprintf(stderr, "3. Error in dangle data format!\n");
    exit(1);
  }
  
  ALPHA_1 = G_ALPHA_1 = (DBL_TYPE) array[0]/100.0;
  ALPHA_2 = G_ALPHA_2 = (DBL_TYPE) array[1]/100.0;
  ALPHA_3 = G_ALPHA_3 = (DBL_TYPE) array[2]/100.0;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //AT PENALTY
  if( sscanf(line, "%d", &(array[0]) ) == 1) {
    AT_PENALTY = G_AT_PENALTY = (DBL_TYPE) array[0]/100.0;
  }
  else {
    fprintf(stderr, "Error in AT PENALTY data\n");
    exit(1);
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x1 interior loop
  for( i = 0; i < 36; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        IL_SInt2[ i*16 + j*4 + k] = G_IL_SInt2[ i*16 + j*4 + k] = (DBL_TYPE) array[ k]/100.0;
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //2x2 interior loop
  for( i = 0; i < 36*16; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 2x2 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        IL_SInt4[ 1536*(i/96)+256*((i%96)/16) +
                 64*((i%16)/4) + 4*(i%4) + k*16 + j] = 
          G_IL_SInt4[ 1536*(i/96)+256*((i%96)/16) +
                     64*((i%16)/4) + 4*(i%4) + k*16 + j] = 
          (DBL_TYPE) array[ k]/100.0;
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x2 interior loop
  for( i = 0; i < 144; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x2 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        IL_AsInt1x2[ 384*(i/24) + 4*((i%24)/4) + 24*(i%4) +
                    96*j + k] = G_IL_AsInt1x2[ 384*(i/24) + 4*((i%24)/4) + 24*(i%4) +
                                              96*j + k] = (DBL_TYPE) array[ k]/100.0;
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //polyC hairpin parameters
  nRead = 0;

  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 3) {
    fprintf(stderr, "4. Error in polyC hairpin parameters!\n");
    exit(1);
  }
  
  POLYC3 = G_POLYC3 = (DBL_TYPE) array[0]/100.0;
  POLYCSLOPE = G_POLYCSLOPE = (DBL_TYPE) array[1]/100.0;
  POLYCINT = G_POLYCINT = (DBL_TYPE) array[2]/100.0;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  
  //Pseudoknot parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "5. Error in dangle data format!\n");
    exit(1);
  }
  
  BETA_1 = G_BETA_1 = (DBL_TYPE) array[0]/100.0;
  BETA_2 = G_BETA_2 = (DBL_TYPE) array[1]/100.0;
  BETA_3 = G_BETA_3 = (DBL_TYPE) array[2]/100.0;
  BETA_1M = G_BETA_1M = (DBL_TYPE) array[3]/100.0;
  BETA_1P = G_BETA_1P = (DBL_TYPE) array[4]/100.0;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //BIMOLECULAR TERM
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 1) {
    fprintf(stderr, "Error in BIMOLECULAR format!\n");
    exit(1);
  }
  
  G_BIMOLECULAR = (DBL_TYPE) array[0]/100.0;
  BIMOLECULAR = G_BIMOLECULAR - kB*TEMP_K*LOG_FUNC( water_conc);
  
  fclose( fp);
  
  /* ****************************** */
  //Load Enthalpies and calculate modified G
  
  //If Temperature == 37 C, add the salt correction andskip he dH file.
  //  Else, make sure it is present
  if( TEMP_K > 37.0 + ZERO_C_IN_KELVIN - 0.001 && TEMP_K < 37.0 + ZERO_C_IN_KELVIN + 0.001) {

    // Make the salt corrections
    // Stacked bases
    for (i = 0; i < 36; i++) {
      Stack[i] += SALT_CORRECTION;
    }
    
    // Loop correction.  Covers all hairpins, and bulges, but overcounts for
    // 1-base bulges.  This is corrected in the function InteriorEnergyFull.
    // This also covers large (non-tabulated) interior loops.
    for (i = 0; i < 90; i++) {
      loop37[i] += SALT_CORRECTION;
    }
    
    // Corrections for tabulated interior loops
    for (i = 0; i < 16*36; i++) {
      IL_SInt2[i] += SALT_CORRECTION;
    }
    for (i = 0; i < 256*36; i++) {
      IL_SInt4[i] += SALT_CORRECTION;
    }
    for (i = 0; i < 64*36; i++) {
      IL_AsInt1x2[i] += SALT_CORRECTION;
    }
    
    // Multiloop
    ALPHA_1 += SALT_CORRECTION;
    
    return;
  }
  
  // Set fileH to be adjacent to fileG up above
  if( ! fileExists( fileH)) {
    fprintf(stderr, "Unable to find file %s.dH\n", fileNameRoot);
    fprintf(stderr, "Consequently, your temperature must be set to 37.0 C, not %.1f.  Job Aborted.\n",
           (float) (TEMP_K - ZERO_C_IN_KELVIN) );
    exit(1);
  }
  
  
  fp = fopen( fileH, "r");
  if( fp == NULL) {  // Make sure input file exits 
    fprintf(stderr, "Error opening loop data file: %s\n", fileH);
    exit(1);  
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Stacking data
  for( i = 0; i < 6; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) ) == 1) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in stacking data format\n");
      exit(1);
    }
    for( j = 0; j < 6; j++) {
      H_Stack[i*6+j] = (DBL_TYPE) array[j]/100.0;
      Stack[i*6+j] = (G_Stack[i*6+j] - H_Stack[i*6+j])*TEMP_K/310.15
        + H_Stack[i*6+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  for( i = 0; i < 3; i++) {
    while( line[0] == '>') {
      fgets( line, MAXLINE, fp);
    }
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead > 30) {
      fprintf(stderr, "Error in Loop energies data\n");
      exit(1);
    }
    for( j = 0; j < 30; j++) {
      if( nRead-1 >= j) {
        H_loop37[30*(2-i)+j] = (DBL_TYPE) array[j]/100.0;
      }
      else {
        H_loop37[30*(2-i)+j] = H_loop37[30*(2-i)+nRead-1]+
          1.75*kB*TEMP_K*LOG_FUNC( (j+1)/(1.0*nRead));
      }
      
      loop37[30*(2-i)+j] = (G_loop37[30*(2-i)+j] - H_loop37[30*(2-i)+j])*
        TEMP_K/310.15 + H_loop37[30*(2-i)+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "Error in asymmetry terms!\n");
    exit(1);
  }
  
  for( j = 0; j < 4; j++) {
    H_asymmetry_penalty[j] = (DBL_TYPE) array[j]/100.0;
    asymmetry_penalty[j] = (G_asymmetry_penalty[j] - 
                            H_asymmetry_penalty[j])*
      TEMP_K/310.15 + H_asymmetry_penalty[j];
  }
  H_max_asymmetry = (DBL_TYPE) array[4]/100.0;
  max_asymmetry = (G_max_asymmetry - 
                   H_max_asymmetry)* TEMP_K/310.15 + H_max_asymmetry;
  
  
  //Triloops
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 2048; i++) {
    H_triloop_energy[i] = 0;
    triloop_energy[i] = (G_triloop_energy[i] - 
                         H_triloop_energy[i])* TEMP_K/310.15 + 
      H_triloop_energy[i];
  }
  while( line[0] != '>') {
    
    if( sscanf( line, "%s %d", triloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 5; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( triloop[4-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( triloop[4-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( triloop[4-i] == 'U' || triloop[4-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( triloop[4-i] != 'A') {
          fprintf(stderr, "Error in triloop indexing %s\n", triloop);
        }
      }
      H_triloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
      triloop_energy[indexL] = (G_triloop_energy[indexL] - 
                                H_triloop_energy[indexL])* TEMP_K/310.15 + 
        H_triloop_energy[indexL];
    }
    else {
      fprintf(stderr, "Error in triloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  //Tetraloops
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  for( i = 0; i < 4096; i++) {
    H_tloop_energy[i] = 0;
    tloop_energy[i] = (G_tloop_energy[i] - 
                       H_tloop_energy[i])* TEMP_K/310.15 + 
      H_tloop_energy[i];
    //printf("%f ", TEMP_K);
  }
  while( line[0] != '>') {
    
    if( sscanf( line, "%s %d", tetraloop, &(array[0]) ) == 2) {
      indexL = 0;
      for( i = 0 ; i < 6; i++) {
        tmpIndex = 1;
        for( j = 0; j < i; j++) {
          tmpIndex *= 4;
        }
        if( tetraloop[5-i] == 'C') {
          indexL += tmpIndex;
        }
        else if( tetraloop[5-i] == 'G') {
          indexL += tmpIndex*2;
        }
        else if( tetraloop[5-i] == 'U' || tetraloop[5-i] == 'T') {
          indexL += tmpIndex*3;
        }
        else if( tetraloop[5-i] != 'A') {
          fprintf(stderr, "Error in tetraloop indexing %s\n", tetraloop);
        }
      }
      
      H_tloop_energy[ indexL] = (DBL_TYPE) array[0]/100.0;
      tloop_energy[indexL] = (G_tloop_energy[indexL] - 
                              H_tloop_energy[indexL])* TEMP_K/310.15 + 
        H_tloop_energy[indexL];
    }
    else {
      fprintf(stderr, "Error in tetraloop data\n%s\n",line);
    }
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Hairpin: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch hairpin format! %d\n", nRead);
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      H_MMEnergiesHP[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
      MMEnergiesHP[6*i+j] = (G_MMEnergiesHP[6*i+j] - 
                             H_MMEnergiesHP[6*i+j])* TEMP_K/310.15 + 
        H_MMEnergiesHP[6*i+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //fprintf(fp, "Mismatch Interior: \n");
  for( i = 0; i < 16; i++) {
    
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 6) {
      fprintf(stderr, "Error in mismatch Interior format!\n");
      exit(1);
    }
    
    for( j = 0; j < 6; j++) {
      H_MMEnergiesIL[ 6*i + j] = (DBL_TYPE) array[j]/100.0;
      MMEnergiesIL[6*i+j] = (G_MMEnergiesIL[6*i+j] - 
                             H_MMEnergiesIL[6*i+j])* TEMP_K/310.15 + 
        H_MMEnergiesIL[6*i+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "6. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      H_dangle_energy[i*4+j] = (DBL_TYPE) array[j]/100.0;
      if( DANGLETYPE == 0) H_dangle_energy[i*4+j] = 0;
      
      dangle_energy[i*4+j] = (G_dangle_energy[i*4+j] - 
                              H_dangle_energy[i*4+j])* TEMP_K/310.15 + 
        H_dangle_energy[i*4+j];
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Read in Dangles
  for( i = 0; i < 6; i++) {
    nRead = 0;
    token = strtok( line, " ");
    while( token != NULL) {
      if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
        nRead++;
      }
      token = strtok( NULL, " ");
    }
    
    if( nRead != 4) {
      fprintf(stderr, "7. Error in dangle data format!\n");
      exit(1);
    }
    for( j = 0; j < 4; j++) {
      H_dangle_energy[24+ i*4+j] = (DBL_TYPE) array[j]/100.0;
      if( DANGLETYPE == 0) 
        H_dangle_energy[24+ i*4+j] = 0;
      
      dangle_energy[24+i*4+j] = (G_dangle_energy[24+i*4+j] - 
                                 H_dangle_energy[24+i*4+j])* TEMP_K/310.15 + 
        H_dangle_energy[24+i*4+j];
      
    }
    
    fgets( line, MAXLINE, fp);
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  //Multiloop parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) ) ==1) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  
  if( nRead != 3) {
    fprintf(stderr, "8. Error in dangle data format!\n");
    exit(1);
  }
  
  H_ALPHA_1 = (DBL_TYPE) array[0]/100.0;
  H_ALPHA_2 = (DBL_TYPE) array[1]/100.0;
  H_ALPHA_3 = (DBL_TYPE) array[2]/100.0;
  
  ALPHA_1 = (G_ALPHA_1 - H_ALPHA_1)* TEMP_K/310.15 + 
    H_ALPHA_1;
  ALPHA_2 = (G_ALPHA_2 - H_ALPHA_2)* TEMP_K/310.15 + 
    H_ALPHA_2;
  ALPHA_3 = (G_ALPHA_3 - H_ALPHA_3)* TEMP_K/310.15 + 
    H_ALPHA_3;
  
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //AT PENALTY
  if( sscanf(line, "%d", &(array[0]) ) == 1) {
    H_AT_PENALTY = (DBL_TYPE) array[0]/100.0;
    AT_PENALTY = (G_AT_PENALTY - H_AT_PENALTY)* TEMP_K/310.15 + 
      H_AT_PENALTY;
    
  }
  else {
    fprintf(stderr, "Error in AT PENALTY data\n");
    exit(1);
  }
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x1 interior loop
  for( i = 0; i < 36; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x1 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        H_IL_SInt2[ i*16 + j*4 + k] = (DBL_TYPE) array[ k]/100.0;
        IL_SInt2[i*16+j*4+k] = (G_IL_SInt2[i*16+j*4+k] - 
                                H_IL_SInt2[i*16+j*4+k])* TEMP_K/310.15 + 
          H_IL_SInt2[i*16+j*4+k];
        
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //2x2 interior loop
  for( i = 0; i < 36*16; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 2x2 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        index4 = 1536*(i/96)+256*((i%96)/16) +
          64*((i%16)/4) + 4*(i%4) + k*16 + j;
        H_IL_SInt4[index4] = 
          (DBL_TYPE) array[ k]/100.0;
        IL_SInt4[index4] = 
          (G_IL_SInt4[index4] - 
           H_IL_SInt4[index4])* TEMP_K/310.15 + 
          H_IL_SInt4[index4];
        
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //1x2 interior loop
  for( i = 0; i < 144; i++) {
    fgets( line, MAXLINE, fp); //read in label
    for( j = 0; j < 4; j++) {
      
      nRead = 0;
      token = strtok( line, " ");
      while( token != NULL) {
        if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
          
          nRead++;
        }
        token = strtok( NULL, " ");
      }
      
      if( nRead != 4) {
        fprintf(stderr, "Error in 1x2 format!\n");
        exit(1);
      }
      
      for( k = 0; k < 4; k++) {
        index4 = 384*(i/24) + 4*((i%24)/4) + 24*(i%4) +
          96*j + k;
        H_IL_AsInt1x2[index4] = (DBL_TYPE) array[ k]/100.0;
        IL_AsInt1x2[index4] = 
          (G_IL_AsInt1x2[index4] - 
           H_IL_AsInt1x2[index4])* TEMP_K/310.15 + 
          H_IL_AsInt1x2[index4];
      }
      
      fgets( line, MAXLINE, fp);
    }
  }
  
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //polyC hairpin parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 3) {
    fprintf(stderr, "9. Error in dangle data format!\n");
    exit(1);
  }
  
  H_POLYC3 = (DBL_TYPE) array[0]/100.0;
  H_POLYCSLOPE = (DBL_TYPE) array[1]/100.0;
  H_POLYCINT = (DBL_TYPE) array[2]/100.0;
  
  POLYC3 = (G_POLYC3 - H_POLYC3)* TEMP_K/310.15 + 
    H_POLYC3;
  POLYCSLOPE = (G_POLYCSLOPE - H_POLYCSLOPE)* TEMP_K/310.15 + 
    H_POLYCSLOPE;
  POLYCINT = (G_POLYCINT - H_POLYCINT)* TEMP_K/310.15 + 
    H_POLYCINT;
  
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //Pseudoknot parameters
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 5) {
    fprintf(stderr, "Error in pseudoknot data format!\n");
    exit(1);
  }
  
  H_BETA_1 = (DBL_TYPE) array[0]/100.0;
  H_BETA_2 = (DBL_TYPE) array[1]/100.0;
  H_BETA_3 = (DBL_TYPE) array[2]/100.0;
  H_BETA_1M = (DBL_TYPE) array[3]/100.0;
  H_BETA_1P = (DBL_TYPE) array[4]/100.0;
  
  BETA_1 = (G_BETA_1 - H_BETA_1)* TEMP_K/310.15 + 
    H_BETA_1;
  BETA_2 = (G_BETA_2 - H_BETA_2)* TEMP_K/310.15 + 
    H_BETA_2;
  BETA_3 = (G_BETA_3 - H_BETA_3)* TEMP_K/310.15 + 
    H_BETA_3;
  BETA_1M = (G_BETA_1M - H_BETA_1M)* TEMP_K/310.15 + 
    H_BETA_1M;
  BETA_1P = (G_BETA_1P - H_BETA_1P)* TEMP_K/310.15 + 
    H_BETA_1P;
  
  fgets( line, MAXLINE, fp);
  while( line[0] == '>') {
    fgets( line, MAXLINE, fp);
  }
  
  //BIMOLECULAR TERM
  nRead = 0;
  token = strtok( line, " ");
  while( token != NULL) {
    if( sscanf( token, "%d", &(array[ nRead]) )==1 ) {
      
      nRead++;
    }
    token = strtok( NULL, " ");
  }
  
  if( nRead != 1) {
    fprintf(stderr, "Error in bimolecular data format!\n");
    exit(1);
  }
  
  H_BIMOLECULAR = (DBL_TYPE) array[0]/100.0;
  
  BIMOLECULAR = (G_BIMOLECULAR - H_BIMOLECULAR)* TEMP_K/310.15 + 
    H_BIMOLECULAR - kB*TEMP_K*LOG_FUNC( water_conc);
  
  fclose( fp);


  // Make the salt corrections
  // Stacked bases
  for (i = 0; i < 36; i++) {
    Stack[i] += SALT_CORRECTION;
  }

  // Loop correction.  Covers all hairpins, and bulges, but overcounts for
  // 1-base bulges.  This is corrected in the function InteriorEnergyFull.
  // This also covers large (non-tabulated) interior loops.
  for (i = 0; i < 90; i++) {
    loop37[i] += SALT_CORRECTION;
  }

  // Corrections for tabulated interior loops
  for (i = 0; i < 16*36; i++) {
    IL_SInt2[i] += SALT_CORRECTION;
  }
  for (i = 0; i < 256*36; i++) {
    IL_SInt4[i] += SALT_CORRECTION;
  }
  for (i = 0; i < 64*36; i++) {
    IL_AsInt1x2[i] += SALT_CORRECTION;
  }

  // Multiloop
  ALPHA_1 += SALT_CORRECTION;

}

/* ************** */

void setParametersToZero( void) {
  
  int i;
	
  for( i = 0; i < 90; i++) {
    loop37[i] = 0;
  }
  for ( i = 0; i < 4096; i++) {
    tloop_energy[i] = 0;
  }
  for( i = 0; i < 2048; i++) {
    triloop_energy[i] = 0;
  }
  for( i = 0; i < 16*6; i++) {
    MMEnergiesHP[i] = 0;
  }
  for( i = 0; i < 256; i++) {
    MMEnergiesIL[i] = 0;
  }
  for( i = 0; i < 16*36; i++) {
    IL_SInt2[i] = 0;
  }
  for( i = 0; i < 256*36; i++) {
    IL_SInt4[ i] = 0;
  }
  for( i = 0; i < 64*36; i++) {
    IL_AsInt1x2[i] = 0;
  }
  for( i = 0; i < 48; i++) {
    dangle_energy[i] = 0;
  }
  for( i = 0; i < 4; i++) {
    asymmetry_penalty[i] = 0;
  }
  max_asymmetry = 0;
	
  for( i = 0; i < 36; i++) {
    Stack[i] = 0;
  }

  ALPHA_1 = ALPHA_2 = ALPHA_3 = BETA_1 = BETA_2 =
    BETA_3 = BETA_1M = BETA_1P = POLYC3 = POLYCINT = POLYCSLOPE =
    AT_PENALTY = BIMOLECULAR = 0;
}

/* *************** */
void InitEtaN( int **etaN, const int *nicks, int seqlength) {
  
  int i,j,k, nick;
  int indexE;
  
  for( i = 0; i <= seqlength-1; i++) {
    for( j = i-1; j <= seqlength-1; j++) {
      indexE = pf_index( i, j, seqlength);
      etaN[ indexE] = (int *) malloc( 2*sizeof( int));
      etaN[ indexE][0] = 0;
      etaN[ indexE][1] = -1;
      
    }
  }
  
  k = 0;
  nick = nicks[k];
  while( nick != -1) {
    for( i = 0; i <= nick; i++) {
      for( j = nick; j <= seqlength-1; j++) { 
        indexE =  pf_index(i,j,seqlength);
        etaN[ indexE][0]++;
        if( etaN[ indexE][1] == -1) { 
          //assume nicks are assigned in increasing order 
          etaN[indexE][1] = k;
        }
      }
    }
    nick = nicks[++k];
  }
}

/* *************************** */
int EtaNIndex_old( float i, float j, int seqlength) { 
  return pf_index( (int) i, (int) j, seqlength);
}

/* ************************** */
void initPF( int seqlength) {
  
  //N^5 or bigger
  static int oldSeqlength = -1;
  extern long int maxGapIndex; 
  
  if( oldSeqlength != seqlength) {
    maxGapIndex = seqlength*(seqlength-1)*(seqlength-2)*(seqlength-3)/24;
    PrecomputeValuesN5( seqlength);
  }
  oldSeqlength = seqlength;
}

/* ************************** */
void initMfe( int seqlength) {
  
  extern long int maxGapIndex; 
  static int oldSeqlength = -1;
  
  if( oldSeqlength != seqlength) {
    maxGapIndex = seqlength*(seqlength-1)*(seqlength-2)*(seqlength-3)/24;
    PrecomputeValuesN5f( seqlength);
  }
  oldSeqlength = seqlength;
  
}

/* ******************* */
