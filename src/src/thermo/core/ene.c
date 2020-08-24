#include "ene.h"

unsigned int seqHash;
int use_cache;
DBL_TYPE ExplDangleRaw( int i, int j, int seq[], int seqlength);

/* ************************************** */

DBL_TYPE HelixEnergy( int i, int j, int h, int m) {
  // Calculate the energy of the helical region closed by pair
  // i-j and h-m.  Data from Zuker's mfold file stack.dgd

  int shift_ij; // Type of base pair
  int shift_hm; // Type of base pair

  extern DBL_TYPE Stack[];

  shift_ij = GetMismatchShift( i, j);
  shift_hm = GetMismatchShift( h, m);

  if( shift_ij < 4 && shift_hm < 4) {
    return Stack[ ( i - 1)*6 + (h - 1) ];
  }

  if( shift_ij < 4 && shift_hm >= 4) {
    return Stack[ (i - 1)*6 + (h + 1) ];
  }

  if( shift_ij >= 4 && shift_hm < 4) {
    return Stack[ (i + 1)*6 + (h - 1) ];
  }

  if( shift_ij >= 4 && shift_hm >= 4) {
    return Stack[ (i + 1)*6 + (h + 1) ];
  }
  else {
    fprintf(stderr, "Error in HelixEnergy!");
    exit(1);
    return NAD_INFINITY; // This never is returned
  }
}

// *******************************************************************
DBL_TYPE InteriorMM( char a, char b, char x, char y) {
/*
  Interior Mismatch calculation

  This calculates the Mismatch interaction energies between positions
  1 -> 5' a x 3'
  2 -> 3' b y 5'
  Interactions energies taken from file tstacki2.dgd.
*/

  extern DBL_TYPE MMEnergiesIL[];
  int cp_shift;
  DBL_TYPE energy = 0.0;

  cp_shift = GetMismatchShift( a, b );
  energy = MMEnergiesIL[ (4*(( x) - 1) + (( y) - 1) )*6 + cp_shift];

  return energy;
}

/* ********************************************** */


DBL_TYPE HairpinEnergy( int i, int j, int seq[] ) {

  // This gives the energy of the hairpion closed by bases i and j
  DBL_TYPE energy;  //energy of hairpin

  int triloopnumber; //Index for specific triloop
  int tloopnumber; //index for tloops

  int size; //Doesn't include closing pair i-j

  int cp_shift; //Classification of base-pair for energy mismatch

  int polyC = TRUE;  //Is the hairpin a poly-C?
  int k;
  for( k = i+1; k < j; k++) {
    if( seq[k] != BASE_C) {
      polyC = FALSE;
      break;
    }
  }

  size = j - i - 1;

  if( size < 3) {
    return NAD_INFINITY;
  }

  if( CanPair( seq[i], seq[j]) == FALSE ) {
    return NAD_INFINITY;
  }

  if( size <= 30) {
    energy = loop37[ 60 + size - 1];
  }
  else {
    energy = loop37[ 60 + 30 - 1];
    energy += sizeLog (size); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);

    if( DNARNACOUNT == COUNT) {
      energy = 0;
    }

  }

  if( size == 3) {
    //Get Triloop energy

    if( seq[i] != BASE_C && seq[j] != BASE_C) {
      energy += AT_PENALTY;
    }

    triloopnumber = 256*(( seq[i]) - 1) +
      64*(( seq[i + 1]) - 1) +
      16*(( seq[i + 2]) - 1) +
      4*( ( seq[j - 1]) - 1) +
      1*( ( seq[j]) - 1);

    // 0 mismatch energy for triloops
    energy += triloop_energy[ triloopnumber];

    //Poly-C loop
    if( polyC == TRUE) {
      energy += POLYC3;
    }
    
    return energy;
  }
  else if (size == 4) {

    tloopnumber = 1024*(( seq[i]) - 1) +
      256*(( seq[i + 1]) - 1) +
      64*( ( seq[i + 2]) - 1) +
      16*( ( seq[j - 2]) - 1) +
      4*(  ( seq[j - 1]) - 1) +
      1*(  ( seq[j])- 1);
    energy +=  tloop_energy[ tloopnumber];
  }
  
  cp_shift = GetMismatchShift( seq[i], seq[j]);

  energy += MMEnergiesHP[(4*(( seq[i + 1]) - 1) +
                          (( seq[j - 1]) - 1) )*6
                         + cp_shift];

  //Poly-C loop
  if( polyC == TRUE) {
    energy += POLYCSLOPE*size + POLYCINT;
  }
  
  return energy;
}



/* ****************************************** */
DBL_TYPE InteriorEnergy(  int i, int j, int h, int m, int seq[]) {
  return InteriorEnergyFull( i, j, h, m, seq, TRUE);
}

DBL_TYPE InteriorEnergyFull( int i, int j, int h, int m, int seq[],
                             int calcIJ) {

  DBL_TYPE energy = 0.0;
  int L1, L2; //lengths of the 2 single stranded regions
  int size;
  int asymmetry;
  int cp_shift, ip_shift;  // For classifying basepairs

  if( DNARNACOUNT == COUNT) return 0;
#ifdef DEBUG
  if( i >= h || h >= m || m >= j) {
    fprintf(stderr, "Invalid boundary to interior loop! %d %d %d %d\n", i, h, m, j);
    exit(1);
  }
#endif

  L1 = h - i - 1;
  L2 = j - m - 1;
  size = L1 + L2;

  if( size == 0) { //Helical region
    energy = HelixEnergy( seq[i], seq[j], seq[h], seq[m] );
  }

  else if ( L1*L2 == 0) { //Bulge
    if( size <= 30) {
      energy = loop37[ 30 + size - 1];
    }
    else {
      energy = loop37[ 30 + 30 - 1];
      energy += sizeLog (size); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);
    }

    if( L1 + L2 == 1 ) { //single bulge...treat as a stacked region
      energy += HelixEnergy( seq[i], seq[j], seq[h], seq[m] );
      energy -= SALT_CORRECTION;  // Correct for the extra salt correction
                                 // added from the HelixEnergy
    }
    else {
      // Next do AT_Penalty for no GC termination, assuming size >= 2
      if( seq[i] != BASE_C && seq[j] != BASE_C) {
        energy += AT_PENALTY;
      }
      if( seq[h] != BASE_C && seq[m] != BASE_C) {
        energy += AT_PENALTY;
      }
    }
  }
  else if ( L1 > 0 && L2 > 0) {
    asymmetry = abs( L1 - L2);
    if( asymmetry > 1 || size > 4) { //data not tabulated

      energy = asymmetryEfn( L1, L2, size);

      //Stacking Energy
      if( L1 > 1 && L2 > 1) { //Non-GAIL Version
        energy += InteriorMM( seq[m], seq[h], seq[m+1], seq[h-1]);

        if( calcIJ == TRUE)
          energy += InteriorMM( seq[i], seq[j], seq[i+1], seq[j-1]);
      }
      else if( L1 == 1 || L2 == 1) {// GAIL =>assume AA terminal mismatch
#ifndef NO_GAIL
        energy +=
          InteriorMM( seq[m], seq[h], BASE_A, BASE_A);
        if( calcIJ == TRUE)
          energy += InteriorMM( seq[i], seq[j], BASE_A, BASE_A);
#else
        energy += InteriorMM( seq[m], seq[h], seq[m+1], seq[h-1]);
        if( calcIJ == TRUE)
          energy += InteriorMM( seq[i], seq[j], seq[i+1], seq[j-1])
#endif
      }
      else {
        fprintf(stderr, "Error: Unclassified interior loop!\n");
        exit(1);
      }
    }
    else { //get tabulated data
      if( asymmetry == 0 && size == 2) {
        cp_shift = GetMismatchShift( seq[i], seq[j]);
        ip_shift = GetMismatchShift( seq[h], seq[m]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        energy += IL_SInt2[ 96*cp_shift + 16*ip_shift +
                           4*(( seq[i+1]) - 1) +
                           (( seq[ j -1]) - 1) ];
      }
      else if( asymmetry == 0 && size == 4) {
        cp_shift = GetMismatchShift( seq[i], seq[j]);
        ip_shift = GetMismatchShift( seq[h], seq[m]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        energy += IL_SInt4[ cp_shift*256*6 +  ip_shift*256 +
                           (4*(( seq[ i+1])  - 1) +
                            ( seq[ j - 1])   - 1)*16 +
                           (4*( ( seq[ i+2]) - 1) +
                            ( seq[ j - 2])   - 1) ];
      }
      else if( asymmetry == 1 && L1 == 1) {
        cp_shift = GetMismatchShift( seq[i], seq[j]);
        ip_shift = GetMismatchShift( seq[h], seq[m]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        energy += IL_AsInt1x2[ cp_shift*4*24*4 +
                              (( seq[ j - 2]) - 1)*24*4 +
                              (( seq[ i + 1]) - 1)*24 +
                              4*ip_shift +
                              ((( seq[ j - 1]) - 1) % 4) ];
      }
      else if( asymmetry == 1 && L1 == 2) {
        cp_shift = GetMismatchShift( seq[j], seq[i]);
        ip_shift = GetMismatchShift( seq[m], seq[h]);
        if (cp_shift==-1 || ip_shift==-1) return 0.0; //Wrongly called
        //note reversed order of inputs above.
        //This is to comply with the format of asint1x2

        energy += IL_AsInt1x2[ ip_shift*4*24*4 +
                              (( seq[i + 1]) - 1)*24*4 +
                              (( seq[j - 1]) - 1)*24 +
                              4*cp_shift +
                              ((( seq[i + 2]) - 1) % 4) ];
      }
      else {
        fprintf(stderr, "Error in tabulated Interior Loop!\n");
        exit(1);
      }
    }
  }
  else {
    fprintf(stderr, "Improperly classified Interior Loop!\n");
    exit(1);
  }

  return energy;
}



/* ********************************************** */
DBL_TYPE DangleEnergyWithPairs( int i, int j, fold *thefold) {

  DBL_TYPE dangle5 = NAD_INFINITY;
  DBL_TYPE dangle3 = NAD_INFINITY;
  int dangle_shift;
  extern DBL_TYPE dangle_energy[];
  int *pairs = thefold->pairs;
  int *seq = thefold->seq;
  int seqlength = thefold->seqlength;

  int pairi1, pairj1;
  int nick3 = 0;
  int nick5 = 0;

  if( i == 0) {
    pairi1 = -1;
  }
  else {
    pairi1 = pairs[i-1];
  }

  if( j == seqlength - 1) {
    pairj1 = -1;
  }
  else {
    pairj1 = pairs[j+1];
  }

  if( j == seqlength - 1) {
    dangle3 = 0;
  }
  else {
    dangle_shift = GetMismatchShift( seq[ pairj1], seq[ j+1]);
#ifdef MATCH_PF
    if( dangle_shift >= 4) {
#ifdef STRUCTURE_WARNINGS
      fprintf(stderr, "Error! This struture not in PF because of wobble %d %d\n", i, j);
#endif
        return NAD_INFINITY;
      exit(1);
    }
#endif
    if( j != -1)
      dangle3 = dangle_energy[ 24 + dangle_shift*4 + ( seq[ j]) - 1];
  }

  if( i == 0) {
    dangle5 = 0;
  }
  else {
    dangle_shift = GetMismatchShift( seq[ i-1], seq[ pairi1]);

#ifdef MATCH_PF
    if( dangle_shift >= 4) {
#ifdef STRUCTURE_WARNINGS
      fprintf(stderr, "Error! This struture not in PF because of wobble- %d %d\n",i,j);
#endif
        return NAD_INFINITY;
    }
#endif
    if( i != seqlength)
      dangle5 = dangle_energy[ dangle_shift*4 + ( seq[ i]) - 1];
  }

  if( DANGLETYPE != 2 && j == i - 1) {
    return 0;
  }
  else if( DANGLETYPE == 2 && j == i - 1 && (i == 0 || j == seqlength - 1) ) {
    return 0;
  }

  if( i != 0 && thefold->isNicked[i-1]) {
    nick5 = 1;
  }
  if(  j != seqlength - 1 && thefold->isNicked[j]) {
    nick3 = 1;
  }

  if( nick5 && nick3) return 0;
  if( nick5 && !nick3) return dangle3;
  if( !nick5 && nick3) return dangle5;

  if( DANGLETYPE == 1 && i == j && i != 0 && j != seqlength - 1) {
    return MIN(dangle3, dangle5 );
  }

  return dangle3 + dangle5;
}

/* ******************************** */
#ifdef COAXIAL
DBL_TYPE CoaxDangle( int whichDangle, int i, int j, int pairs[], int seq[], int seqlength) {

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;
  extern DBL_TYPE dangle_energy[];

  int pairi1 = pairs[i-1];
  int pairj1 = pairs[j+1];

#ifdef MATCH_PF
  fprintf(stderr, "Coaxially Stacking needs to be off to match PF calculations!\n");
  exit(1);
#endif

#ifndef VIENNA_D2
  if( j == i - 1) {
    return 0;
  }

#else
  if( j == i - 1 && (i == 0 || j == seqlength - 1) ) {
    return 0;
  }
#endif

  if( j == seqlength - 1) {
    dangle3 = 0;
  }
  else {
    dangle_shift = GetMismatchShift( seq[ pairj1], seq[ j+1]);
#ifdef MATCH_PF
    if( dangle_shift >= 4) {
#ifdef STRUCTURE_WARNINGS
      fprintf(stderr, "Error! This struture not in PF because of wobble %d %d\n", i, j);
#endif
        return NAD_INFINITY;
      //exit(1);
    }
#endif
    dangle3 = dangle_energy[ 24 + dangle_shift*4 + ( seq[ j]) - 1];
  }

  if( i == 0) {
    dangle5 = 0;
  }
  else {
    dangle_shift = GetMismatchShift( seq[ i-1], seq[ pairi1]);
#ifdef MATCH_PF
    if( dangle_shift >= 4) {

#ifdef STRUCTURE_WARNINGS
      fprintf(stderr, "Error! This struture not in PF because of wobble- %d %d\n",i,j);
#endif
        return NAD_INFINITY;
    }
#endif

    dangle5 = dangle_energy[ dangle_shift*4 +
                            ( seq[ i]) - 1];
  }

  if( whichDangle == 3) {
    return dangle3;
  }

  if( whichDangle == 5) {
    return dangle5;
  }

  if( whichDangle == 53) {

#ifndef VIENNA_D2
    if( i == j && i != 0 && j != seqlength - 1) {
      return MIN(dangle3, dangle5 );
    }
    else {
      return dangle3 + dangle5;
    }
#else
    return dangle3 + dangle5;
#endif
  }
  else {
    fprintf(stderr, "Invalid whichDangle Value of %d in CoaxDangle\n", whichDangle);
    exit(1);
  }
}
#endif

/* ******************************** */
DBL_TYPE DangleEnergy( int i, int j, int seq[], int seqlength) {
  //0 energy except for dangles

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;

  if( DANGLETYPE != 2) {
    if( j == i - 1) {
      return 0;
    }
  }
  else if( j == i - 1 && (i == 0 || j == seqlength - 1) ) {
    return 0;
  }

  if( j == seqlength - 1) {
    dangle3 = 0;
  }
  else {
    int pt=GetPairType( seq[ j + 1]);
    if (pt==-1) {
      printf("i=%d j=%d seq[%d]=%d\n",i,j,j-1,seq[j-1]);
      exit(-1);
    }
    dangle_shift = 3 - pt;
    dangle3 = dangle_energy[ 24 + dangle_shift*4 +
                            ( seq[ j]) - 1];
  }

  if( i == 0) {
    dangle5 = 0;
  }
  else {
    int pt=GetPairType( seq[ i - 1]);
    if (pt==-1) {
      printf("i=%d j=%d seq[%d]=%d\n",i,j,i-1,seq[i-1]);
      exit(-1);
    }

    dangle_shift = pt;
    dangle5 = dangle_energy[ dangle_shift*4 +
                            ( seq[ i]) - 1];
  }

  if( DANGLETYPE != 2 && i == j && i != 0 && j != seqlength - 1) {
    return MIN(dangle3, dangle5 );
  }

  return dangle3 + dangle5;
}

/* ******************************** */
DBL_TYPE ExplDangleRaw( int i, int j, int seq[], int seqlength) {
  //0 energy except for dangles

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;

  if( (j == i - 1) || (j==-1 && i>0)) {
    return 1.0;
  }
  if( (j==-1 && i>0) || (j == i - 1 && (i == 0 || j == seqlength - 1)) ) {
    return 1.0;
  }

  if( j == seqlength - 1) {
    dangle3 = 0;
  }
  else {
    dangle_shift = 3 - GetPairType( seq[ j + 1]);
    dangle3 = dangle_energy[ 24 + dangle_shift*4 +
                            seq[ j] - 1];
  }

  if( i == 0) {
    dangle5 = 0;
  }
  else {
    dangle_shift = GetPairType( seq[i-1]);
    dangle5 = dangle_energy[ dangle_shift*4 +
                            seq[ i] - 1];
  }

  if(i == j && i != 0 && j != seqlength - 1) {
    if (DANGLETYPE == 2) return EXP_FUNC(-(dangle5 + dangle3)/(TEMP_K*kB));
    return EXP_FUNC( -MIN(dangle3, dangle5)/(TEMP_K*kB) );
  }

  return EXP_FUNC( -(dangle3 + dangle5)/(TEMP_K*kB) );
}

DBL_TYPE ExplDangle( int i, int j, int seq[], int seqlength) {
  static DBL_TYPE *EDCache=NULL;
  static int CacheInd=-1, nCaches=0, dangleTypeCache=2;
  static DBL_TYPE TCache;
  static unsigned int SCache=1;

  if (!use_cache) return ExplDangleRaw(i,j,seq,seqlength);
  if (CacheInd==-1 || SCache!=seqHash || TCache!=TEMP_K
      || dangleTypeCache!=DANGLETYPE){ // We got a new sequence or temp
    int d,e;
    if (CacheInd!=-1 || EDCache)
      free(EDCache);
    EDCache=(DBL_TYPE *)calloc((seqlength+1)*(seqlength+1),sizeof(DBL_TYPE));
    if (!EDCache){
      use_cache=0;
      fprintf(stderr, "ExplDangle: unable to allocate %lu bytes, disabling cache\n",(unsigned long)(seqlength+1)*(seqlength+1)*sizeof(DBL_TYPE));
      return ExplDangleRaw(i,j,seq,seqlength);
      //exit(0);
    }
    CacheInd=nCaches++;

    seqHash=vechash((char *)seq, (unsigned int) seqlength*sizeof(int));
    SCache=seqHash;

    TCache=TEMP_K;
    dangleTypeCache=DANGLETYPE;

    for (d=0;d<seqlength+1;d++){
      for (e=-1;e<seqlength;e++){
        EDCache[d*(seqlength+1)+e+1]=ExplDangleRaw(d,e,seq,seqlength);
      }
    }
  }


  return EDCache[(i*(seqlength+1)+(j+1))];
}


/* *********** */
DBL_TYPE NickDangle(int i, int j, const int *nicks, int **etaN, int hairpin,
                    int seq[], int seqlength) {

  DBL_TYPE dangle5 = 0;
  DBL_TYPE dangle3 = 0;
  int dangle_shift;
  extern DBL_TYPE dangle_energy[];
  int nick;
  int nIndex;

  nick = -5;

  if( i != 0) { //if j == seqlength -1, this is still OK
    nIndex = EtaNIndex( i-0.5, j+0.5, seqlength);
  }
  else {
    nIndex = EtaNIndex( i+0.5, j+0.5, seqlength);
  }

  if( etaN[ nIndex][0] >= 2 ||
     ( etaN[ nIndex][0] == 1 && (i == 0 || j == seqlength -1)) ) {

       return NAD_INFINITY;
     }
  else if( etaN[ nIndex][0] >= 1) {
    nick = nicks[ etaN[ nIndex][1]];
  }

  if( DNARNACOUNT == COUNT)
    return 0;

  if( j == i - 1) {
    return 0;
  }
  if( j == i - 1 && (i == 0 || j == seqlength - 1) ) {
    return 0;
  }

  if( j == seqlength - 1 || j == nick) {
    dangle3 = 0;
  }
  else {
    if( hairpin == FALSE) {
      dangle_shift = 3 - GetPairType( seq[ j + 1]);
    }
    else {
      dangle_shift = GetMismatchShift( seq[i-1], seq[j+1]);
    }

    dangle3 = dangle_energy[ 24 + dangle_shift*4 +
                            ( seq[ j]) - 1];
  }

  if( i == 0 || i-1 == nick) {
    dangle5 = 0;
  }
  else {
    if( hairpin == FALSE) {
      dangle_shift = GetPairType( seq[i-1]);
    }
    else {
      dangle_shift = GetMismatchShift( seq[i-1], seq[j+1]);
    }

    dangle5 = dangle_energy[ dangle_shift*4 +
                            ( seq[ i]) - 1];
  }

  if( nick >= i-1 && nick <= j) {
    return dangle3 + dangle5;
  }
  else {
    if( j > i || j == seqlength - 1 || i == 0) {
      return dangle3 + dangle5;
    }
    if(j == i && i != 0 && j != seqlength - 1) {
      if (DANGLETYPE == 2) return dangle3 + dangle5;
      return MIN(dangle3, dangle5);
    }
    //j == i-1 already handled above
  }
  fprintf(stderr, "Error with function: NickDangle\n");
  exit(-1);
  return -1; //Error!  This should never happen
}

/* ************** */
DBL_TYPE NickedEmptyQ( int i, int j, const int nicks[], int seq[],
                      int seqlength, int **etaN) {

  if( j <= i || etaN[ EtaNIndex( i+0.5, j-0.5, seqlength)][0] == 0) {
    return EXP_FUNC( -1*NickDangle(i, j, nicks, etaN,
                              FALSE, seq, seqlength)
                /(kB*TEMP_K));
  }
  else { //disconnected
    return 0;
  }

}

/* ******** */
DBL_TYPE NickedEmptyF( int i, int j, const int nicks[], int seq[],
                       int seqlength, int **etaN) {
  DBL_TYPE result = NAD_INFINITY;

  if( j <= i || etaN[ EtaNIndex( i+0.5, j-0.5, seqlength)][0] == 0) {
    result = NickDangle(i, j, nicks, etaN, FALSE, seq, seqlength);
  }

  return result;
}

/* ********* */
DBL_TYPE ExplInternal( int i, int j, int h, int m, int seq[]) {
  // Calculates E^(-energy/RT) of interior loop closed by i-j and h-m

  DBL_TYPE energy = InteriorEnergy( i, j, h, m, seq);
  if( energy == NAD_INFINITY) {
    return 0.0;
  }
  return EXP_FUNC( - energy/( kB*TEMP_K));
}

DBL_TYPE sizeLog(int size){
  static DBL_TYPE *slCache[MAXSTRANDS], *edc, tc;
  static int CacheInd=-1, nCaches=0;
  static DBL_TYPE TCache[MAXSTRANDS];

  if (CacheInd==-1 || tc!=TEMP_K){ // We got a new sequence or temp
    static unsigned int keySize=sizeof(int)+sizeof(DBL_TYPE);
    char key[sizeof(int*)+sizeof(DBL_TYPE)];
    static int IndCache[MAXSTRANDS];
    static void *indP=NULL;
    static hash *expHash;
    int d;

    if (CacheInd==-1) { // We need to create a new hash
      expHash=hash_new((unsigned int) MAXSTRANDS);
    }

    // Calculate key
    memcpy((void*)key, (void *)&size,sizeof(int));
    memcpy((void*)(key+sizeof(int)), (void *)&TEMP_K ,sizeof(DBL_TYPE));

    // Search for key
    indP=hash_get (expHash, key, keySize);
    if (!indP){
      CacheInd=nCaches++;
      IndCache[CacheInd]=CacheInd;
      hash_add(expHash, key, keySize, (void *)&IndCache[CacheInd]);

      slCache[CacheInd]=(DBL_TYPE *)malloc(MAXSEQLENGTH*sizeof(DBL_TYPE));
      TCache[CacheInd]=TEMP_K;

      slCache[CacheInd][0]=0.0;
      for (d=1;d<MAXSEQLENGTH;d++){
          slCache[CacheInd][d]=1.75*kB*TEMP_K*LOG_FUNC( d/30.0);
      }
      indP=hash_get (expHash, key, keySize);
    }
    CacheInd=*(int*)indP;
    edc=slCache[CacheInd];
    tc=TCache[CacheInd];
  }

  return edc[size];

}

DBL_TYPE sizeEnergyLog(int size){
  static DBL_TYPE *slCache[MAXSTRANDS], *edc, tc;
  static int CacheInd=-1, nCaches=0;
  static DBL_TYPE TCache[MAXSTRANDS];

  if (CacheInd==-1 || tc!=TEMP_K){ // We got a new sequence or temp
    static unsigned int keySize=sizeof(int)+sizeof(DBL_TYPE);
    char key[sizeof(int*)+sizeof(DBL_TYPE)];
    static int IndCache[MAXSTRANDS];
    static void *indP=NULL;
    static hash *expHash;
    int d;

    if (CacheInd==-1) { // We need to create a new hash
      expHash=hash_new((unsigned int) MAXSTRANDS);
    }

    // Calculate key
    memcpy((void*)key, (void *)&size,sizeof(int));
    memcpy((void*)(key+sizeof(int)), (void *)&TEMP_K ,sizeof(DBL_TYPE));

    // Search for key
    indP=hash_get (expHash, key, keySize);
    if (!indP){
      CacheInd=nCaches++;
      IndCache[CacheInd]=CacheInd;
      hash_add(expHash, key, keySize, (void *)&IndCache[CacheInd]);

      slCache[CacheInd]=(DBL_TYPE *)malloc(MAXSEQLENGTH*sizeof(DBL_TYPE));
      TCache[CacheInd]=TEMP_K;
      DBL_TYPE oldSizeEnergy, newSizeEnergy;
      for (d=10;d<MAXSEQLENGTH;d++){
        if( d <= 30) {
          oldSizeEnergy = loop37[ d - 1];
        }
        else {
          oldSizeEnergy = loop37[ 30 - 1];
          oldSizeEnergy += sizeLog (d); //1.75*kB*TEMP_K*LOG_FUNC( size/30.0);
        }

        if( d - 2 <= 30) {
          newSizeEnergy = loop37[ d-2 - 1];
        }
        else {
          newSizeEnergy = loop37[ 30 - 1];
          newSizeEnergy += sizeLog (d-2); //1.75*kB*TEMP_K*LOG_FUNC( (size-2)/30.0);
        }

        slCache[CacheInd][d]=EXP_FUNC( -(newSizeEnergy - oldSizeEnergy)/(kB*TEMP_K));
      }
      indP=hash_get (expHash, key, keySize);
    }
    CacheInd=*(int*)indP;
    edc=slCache[CacheInd];
    tc=TCache[CacheInd];
  }

  return edc[size];

}


/* ******* */
DBL_TYPE asymmetryEfn( int L1, int L2, int size) {

  int asymmetry_index;
  DBL_TYPE energy;
  int asymmetry = abs( L1 - L2);

  //Loop Size Energy
  if( size <= 30) {
    energy = loop37[ size - 1];
  }
  else {
    energy = loop37[ 30 - 1];
    energy += sizeLog(size);
  }

  //Asymmetry rountine copied from efn.f in Zuker's mfold package.
  asymmetry_index = 4;

  if( L1 < asymmetry_index) {
    asymmetry_index = L1;
  }

  if( L2 < asymmetry_index) {
    asymmetry_index = L2;
  }

  if( asymmetry*asymmetry_penalty[ asymmetry_index - 1] < max_asymmetry ) {
    energy += asymmetry*asymmetry_penalty[ asymmetry_index - 1];
  }
  else {
    energy += max_asymmetry; // MAX asymmetry penalty
  }
  return energy;
}



/* ********************** */
