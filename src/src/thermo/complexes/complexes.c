/*
   complexes.c is part of the NUPACK software suite
   Copyright (c) 2007 Caltech. All rights reserved.
   Coded by: Robert Dirks, 3/2006 and Justin Bois 1/2007

   Multistranded partition function
   Calculates all distinct complexes and
   circular permutations for N strands.

   Computes and outputs all of the partition functions

   See NUPACK manual for usage notes
*/


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

#include <thermo/core.h>

#include "complexesStructs.h"
#include "complexesUtils.h"
#include "permBG.h"
#include "ReadCommandLine.h"

//Global variables
extern int nStrands;

globalArgs_t globalArgs;

int main( int argc, char **argv) {
  //generate subsets for given sequences

  int i, j, k, q, r, l;  // Counters

  int nPercent = 40; // Number of percent signs in a line with all comments

  char **seqs; //list of all seqs
  int *seqlength; //list of all seqlengths

  multiset *allSets;

  char firstChar[2];

  int nSets;
  int totalSets; // including sets in the input list
  int setStart=0;  // index of set to start with, for use with -listonly

  int maxLength;

  int maxComplexSize=0;
  int maxListComplexSize=0;
  int complexSize;
  int nMonomer;

  int nNewComplexes;
  int nNewPerms;

  char *pfSeq;
  int tmpLength; // Stores current sequence length
  int seqNum[ MAXSEQLENGTH+1]; // Stores current sequence as ints
  long double pf;
  int *nicks;

  //int done;

  double totalOrders;
  int nTotalOrders;

  int totalOrders2 = 0;
  int status;

  int seqCode;
  FILE *F_cx = NULL, *F_ocx = NULL, // *.cx, only output with -v3
  *F_list = NULL, *F_perm = NULL,
  *F_permPr=NULL,
  //*F_permAvgP=NULL,
  *F_cxAvg=NULL, // *.cx-epairs, only output with -v3
  //*F_cxAvgP=NULL,
  *F_permMfe = NULL,
  //*F_cxMfe = NULL,
  *F_permAvg = NULL,
  *F_prog = NULL,
  *F_defect = NULL;
  char filePrefix[100], cxName[110], ocxName[110],
  listName[110], permName[110],
  permPrName[110],
  //permAvgPName[110],
  cxAvgName[110],
  //cxAvgPName[110],
  permMfeName[110],
  //cxMfeName[110],
  permAvgName[110],
  progName[110],
  defectName[110];

  char *token;

  double prog = 0.0, progDenom = 0.0; // Used in calculating progress

  long double TEMP_K;

  int strandId, strandId2;
  int strandPos, strandPos2;

  permutation *currentPerm=NULL;
  permutation *tmpPerm;

  //time_t start, end;
  char line[ MAXLINE];
  char line2[MAXLINE];
  int fileRead;
  int inputFileSpecified;

  int permId;
  long double **permPr = NULL;
  int indQ;
  int pf_qr;
  int **indexCx = NULL;
  int totalStrandsLength=0; //length of all the strand types, concatenated
  int row, column;
  long double expectedUnpaired;

  //permAvgPairs and cxAvgPairs store the expected value of each
  //class of base pairs, grouped by permutation or complex, respectively
  long double *permAvgPairs = NULL;
  long double *cxAvgPairs=NULL;
  long double ppf;

  double estimatedTime = 0;
  double N3C = 1.0e-6; //estimatedTime = N3C * seqlength^3

  time_t curtime; //for printing date and time
  struct tm *loctime;

  dnaStructures mfes = {NULL, 0, 0, 0, 10000};
  long double mfeValue;

  int lastCxId = 1; //used to index complex id#,
  //in case some are not used due to no possible secondary structures

  // Set defaults of global args
  globalArgs.quiet = 0;
  globalArgs.permsOn = 1;
  globalArgs.dopairs = 0;
  globalArgs.T = 37.0;
  globalArgs.dangles = 1;
  globalArgs.parameters = RNA;
  globalArgs.out = 1; //.cx file
  globalArgs.timeonly = 0;
  globalArgs.listonly = 0;
  globalArgs.debug = 0;
  globalArgs.echo = 0;
  globalArgs.mfe = 0;
  globalArgs.cutoff = 0.001; // Cutoff bp probability to report
  globalArgs.progress = 0;
  globalArgs.onlyOneMFE = 1;
  globalArgs.sodiumconc = 1.0;
  globalArgs.magnesiumconc = 0.0;
  globalArgs.uselongsalt = 0;
  globalArgs.dodefect = 0;
  strcpy( globalArgs.inputFilePrefix, "NoInputFile");

  /* version 3 output */
  globalArgs.v3 = 0;

  inputFileSpecified = ReadCommandLine( argc, argv);
  

  if (globalArgs.dodefect == 1) {
    globalArgs.mfe = 1;
    globalArgs.permsOn = 1;
    globalArgs.dopairs = 1;
  }
  if (globalArgs.mfe == 1 && globalArgs.permsOn == 0) {
    if( !globalArgs.quiet) {
      printf("Warning, including -mfe requires -ordered. ");
      printf("As a result, -ordered has been enabled.\n");
    }
    globalArgs.permsOn = 1;
  }

  if (!inputFileSpecified) {
    printf("No input file specified.\n");
    fileRead = 0;
  }
  else {
    fileRead = ReadInputFileComplexes( filePrefix, &nStrands,
                                      &seqs, &seqlength, &maxLength,
                                      &maxComplexSize);
    if (fileRead == 2) {
      printf("Input file %s.in not found.\n",filePrefix);
      fileRead = 0;
    }
  }

  if( !fileRead ) {
    printf("Requesting input manually.\n");
    printf("Enter file prefix: ");
    scanf("%s", filePrefix);
  }

  TEMP_K = globalArgs.T + ZERO_C_IN_KELVIN;
  
  /* version 3 output */
  if(globalArgs.v3) {
    sprintf( cxName, "%s.cx", filePrefix);
    sprintf( cxAvgName, "%s.cx-epairs", filePrefix);
  }

  sprintf( ocxName, "%s.ocx", filePrefix);
  sprintf( listName, "%s.list", filePrefix);
  sprintf( permName, "%s.ocx-key", filePrefix);
  sprintf( permPrName, "%s.ocx-ppairs", filePrefix);
  sprintf( permMfeName, "%s.ocx-mfe", filePrefix);
  sprintf( permAvgName, "%s.ocx-epairs", filePrefix);
  sprintf( progName, "%s.prog", filePrefix);
  sprintf( defectName, "%s.ocx-defect", filePrefix);
  
  /* version 3 output */
  char *composition = NULL;
  char *ordering = NULL;
  if (globalArgs.v3) {
    composition = "complex";
    ordering = "order";
  }
  else {
    composition = "composition";
    ordering = "ordering";
  }

  if( globalArgs.out == 1) {
    /* version 3 output */
    if(globalArgs.v3) {
      F_cx = fopen( cxName, "w");
      if( !F_cx) printf("Error: Unable to create %s\n", cxName);
    }

    if( globalArgs.mfe) {
      F_permMfe = fopen( permMfeName, "w");
      if( !F_permMfe )
        printf("Error: Unable to create %s\n", permMfeName);
    }

    if( globalArgs.dopairs) {
      /* version 3 output */
      if(globalArgs.v3) {
        F_cxAvg = fopen( cxAvgName, "w");
        if( !F_cxAvg) printf("Error: Unable to create %s\n", cxAvgName);
      }

      if( globalArgs.permsOn) {
        F_permPr = fopen( permPrName, "w");

        F_permAvg = fopen( permAvgName, "w");
        if( !F_permPr || !F_permAvg)
          printf("Error: Unable to create %s or %s\n",
                 permPrName,
                 permAvgName);
      }
    }
    if( globalArgs.dodefect) {
      F_defect = fopen(defectName,"w");
      if(!F_defect) {
        printf("Error: unable to create %s\n",
                defectName);
      }
    }


    if ( globalArgs.progress) {
      F_prog = fopen(progName,"w");
      if (!F_prog) {
        printf("Error: Unable to create %s.\n",progName);
      }
      else {
        prog = 0.0;
        if(!NUPACK_VALIDATE) {
          fprintf(F_prog,"%.4f\n\n",prog); // Second newline is necessary for webserver
        } else {
          fprintf(F_prog,"%.14f\n\n",prog); // Second newline is necessary for webserver
        }
        fclose(F_prog);
      }
    }
  }

  F_list = fopen( listName, "r");

  if( F_list == NULL && !(globalArgs.quiet)) {  // check if file exists
    printf("There is no input list %s.\n", listName);
  }

  if( globalArgs.permsOn && globalArgs.out) {
    F_perm = fopen( permName, "w");
    if( !F_perm) printf("Error: Unable to create %s\n", permName);

    F_ocx = fopen( ocxName, "w");
    if( !F_ocx) printf("Error: Unable to create %s\n", ocxName);
  }


  if( !fileRead) { //if error in reading input file, get manual input
    printf("Enter Number of Different Sequences: ");
    scanf("%d", &nStrands);

    //allocate function variables
    seqs = (char **) malloc( nStrands*sizeof( char*));
    seqlength = (int *) malloc( nStrands*sizeof( int));

    maxLength = 0;
    for( i = 0; i<=nStrands-1; i++) {
      printf("Enter Sequence %d:\n", i+1);
      scanf("%s", line);
      seqlength[i] = strlen( line);
      if( seqlength[i] > maxLength) maxLength = seqlength[i];

      seqs[i] = (char*) malloc( (seqlength[i]+1)*sizeof( char));
      strcpy( seqs[i], line);
    }

    printf("Enter max complex size to completely enumerate: ");
    scanf("%d", &maxComplexSize);
  }

  if( globalArgs.dopairs) {
    indexCx = (int **) malloc( nStrands*sizeof( int*));

    totalStrandsLength = 0;
    for( i = 0; i < nStrands; i++) {
      if( globalArgs.dopairs) {
        indexCx[i] = (int*) malloc( seqlength[i]*sizeof( int) );
        for( j = 0; j < seqlength[i]; j++) {
          indexCx[i][j] = totalStrandsLength++;
        }
      }
    }

    //these will store the average number of pairs for a particular base type
    cxAvgPairs = (long double *) malloc (totalStrandsLength*
                                         totalStrandsLength*
                                         sizeof( long double));
    permAvgPairs = (long double *) malloc (totalStrandsLength*
                                           totalStrandsLength*
                                           sizeof( long double));

    if( !cxAvgPairs || !permAvgPairs) {
      printf("Unable to allocate cxAvgPairs or permAvgPairs!\n");
      exit(1);
    }
  }

  //time( &start);


  // Read information from .list file
  maxListComplexSize = maxComplexSize;
  //printf("maxComplexSize=%d\n",maxComplexSize);
  nNewComplexes = nNewPerms = 0;
  int CLastLine = 0;
  int * lastCLine = malloc(nStrands * sizeof(int));
  if( F_list != NULL) {
    while( fgets( line, MAXLINE, F_list)) {
      sscanf(line, "%1s", firstChar);

      if( firstChar[0] == '%') {
        continue;
      }
      if( firstChar[0] == 'C') {
        if(CLastLine) {
          nNewPerms += makeFCPermutations(NULL,lastCLine,complexSize,nStrands);
        }
        nNewComplexes++;
        complexSize = 0;
        token = strtok( line, " ");
        CLastLine = 1;
        for( i = 0; i < nStrands; i++) {
          token = strtok( NULL, " ");
          if( !token ||  sscanf( token, "%d", &nMonomer) != 1) {
            printf("Error in list file, line: %s\n", line);
            exit(1);
          }
          lastCLine[i] = nMonomer;

          complexSize += nMonomer;
        }

        // Check to make sure input complex size is greater than max complex size
        if (complexSize <= maxComplexSize) {
          printf("Error in list file.  All complexes in list file must be\n");
          printf("greater than maximum complex size specified in input file.\n");
          exit(1);
        }
        if( complexSize > maxListComplexSize)
          maxListComplexSize = complexSize;
      }
      else if( isdigit( firstChar[0]) ) {
        CLastLine = 0;
        nNewPerms++;
      }
    }
    if(CLastLine && complexSize > maxComplexSize) {
      nNewPerms += makeFCPermutations(NULL,lastCLine,complexSize,nStrands);
    }
  }

  //determine total # of distinct strand orders (lovasz, 3.23b)
  totalOrders = nNewPerms;
  for( i = 1; i <= maxComplexSize; i++) {
    for( j = 1; j <= i; j++) {

      totalOrders += pow( nStrands, gcd( j, i))/i;
    }
  }
  nTotalOrders = totalOrders + 0.1;

  //Next, generate all multisets
  nSets = binomial_coefficient(maxComplexSize + nStrands,maxComplexSize) - 1;
  if(nSets < 1) {
    fprintf(stderr,"Integer overflow occurred while counting permutations!\n");
    exit(1);
  }

  totalSets = nSets + nNewComplexes;

  if( globalArgs.out == 1) {
    
    /* version 3 output */
    if(globalArgs.v3) {
      printHeader( nStrands, seqs, maxComplexSize, nTotalOrders,
                  nNewPerms, nSets, nNewComplexes, F_cx, argc, argv, 0);
    }
    if( globalArgs.permsOn) {
      printHeader( nStrands, seqs, maxComplexSize, nTotalOrders,
                  nNewPerms, nSets, nNewComplexes, F_ocx, argc, argv, 0);
      printHeader( nStrands, seqs, maxComplexSize, nTotalOrders,
                  nNewPerms, nSets,
                  nNewComplexes, F_perm, argc, argv, 0);
    }

    if( globalArgs.mfe) {
      printHeader( nStrands, seqs, maxComplexSize,
                  nTotalOrders, nNewPerms, nSets,
                  nNewComplexes, F_permMfe, argc, argv, 0);
    }

    if( globalArgs.dopairs) {
      /* version 3 output */
      if(globalArgs.v3) {
        printHeader( nStrands, seqs, maxComplexSize, nTotalOrders,
                    nNewPerms, nSets,
                    nNewComplexes, F_cxAvg, argc, argv, 1);
      }
      if( globalArgs.permsOn) {
        printHeader( nStrands, seqs, maxComplexSize, nTotalOrders,
                    nNewPerms, nSets,
                    nNewComplexes, F_permPr, argc, argv, 1);
        printHeader( nStrands, seqs, maxComplexSize, nTotalOrders,
                    nNewPerms, nSets,
                    nNewComplexes, F_permAvg, argc, argv, 1);

      }
    }
    if(globalArgs.dodefect) {
      printHeader(nStrands, seqs, maxComplexSize, nTotalOrders,
                  nNewPerms, nSets, nNewComplexes, F_defect, argc, argv, 1);
      fprintf(F_defect,"%%\n");
      fprintf(F_defect,"%% Ensemble defects n(s,phi) and normalized ensemble defects n(s,phi)/N\n");
      fprintf(F_defect,"%% Comp\tPerm\t");
      int i_strand;
      for(i_strand = 0 ; i_strand < nStrands; i_strand++) {
        fprintf(F_defect, "S%i\t", i_strand+1);
      }
      fprintf(F_defect,"n(s,phi)\tn(s,phi)/N\tStruc Prob\tPart Func\n");
    }
  }

  // Generate all necklaces for each length with order nStrands

  int cursize = 1;
  permutation * allPermutations = (permutation*)
                    malloc(nTotalOrders * sizeof(permutation));
  int added = 0;
  int offset = 0;
  for(cursize = 1; cursize <= maxComplexSize ; cursize++) {
    added = makePermutations(allPermutations + offset,cursize,nStrands);
    offset += added;
  }
  totalOrders2 = offset;

  if( F_list != NULL) {
    rewind( F_list); //go to beginning of prefix.list

    //Get extra sets and perms from prefix.list
    offset = totalOrders2;
    int curStrand;
    int curStrandIndex;
    int curNumStrands;
    int list_file_line = 0;
    int warning_printed = 0;
    CLastLine = 0;

    while(NULL != fgets(line,MAXLINE,F_list)) {
      list_file_line ++;
      strncpy(line2,line,MAXLINE);
      token = strtok(line," ,\t\n");
      if(token == NULL) {
        continue;
      }
      if(token[0] == '%') {
        continue;
      }
      if(token[0] == 'C') {
        if(!warning_printed) {
          printf("Warning: 'C' lines in the list file are deprecated please specify ordered complexes directly\n");
          warning_printed = 1;
        }

        if(CLastLine) {
          curNumStrands = 0;
          for(curStrand = 0 ; curStrand < nStrands ; curStrand++) {
            curNumStrands += lastCLine[curStrand];
          }
          offset += makeFCPermutations(allPermutations + offset, lastCLine, curNumStrands, nStrands);
        }

        curStrand = 0;
        token = strtok(NULL," ,\t\n");
        while(NULL != token && curStrand < nStrands) {
          if(! sscanf(token, "%d", &curNumStrands)) {
            fprintf(stderr,"Error on line %d: encountered %s\n",list_file_line,token);
            exit(1);
          }
          lastCLine[curStrand] = curNumStrands;
          token = strtok(NULL," ,\t\n");
          curStrand++;
        }

        CLastLine = 1;
        continue;
      }
      
      if(0 == sscanf(token, "%d", &curStrand)) {
        fprintf(stderr,"Error on line %d: encountered %s\n",list_file_line,token);
        exit(1);
      }
      
      curNumStrands = 0;
      while(NULL != token) {
        if(0 == sscanf(token, "%d", &curStrand)) {
          fprintf(stderr,"Error on line %d: encountered %s\n",list_file_line,token);
          exit(1);
        }
        curNumStrands++;
        token = strtok(NULL, " ,\t\n");
      }
      
      allPermutations[offset].nSeqs = curNumStrands;
      allPermutations[offset].code = (int *) malloc(curNumStrands * sizeof(int));
      allPermutations[offset].strand_sums = (int *) malloc(nStrands * sizeof(int));
      allPermutations[offset].symmetryFactor = 1;
      for(curStrand = 0; curStrand < nStrands ; curStrand++) {
        allPermutations[offset].strand_sums[curStrand] = 0;
      }
      token = strtok(line2," ,\t\n");

      curStrandIndex = 0;
      CLastLine = 0;
      while(NULL != token) {
        sscanf(token, "%d", &curStrand);
        if(curStrand > nStrands) {
          fprintf(stderr,"Error on line %d: %i > number of strands\n",list_file_line,curStrand);
          exit(1);
        }
        allPermutations[offset].code[curStrandIndex] = curStrand;
        allPermutations[offset].strand_sums[curStrand - 1] ++;
        curStrandIndex ++;
        token = strtok(NULL, " ,\t\n");
      }
      offset++;
    }
    if(CLastLine ) {
      curNumStrands = 0;
      for(curStrand = 0; curStrand < nStrands ; curStrand++) {
        curNumStrands += lastCLine[curStrand];
      }
      offset += makeFCPermutations(allPermutations+offset,lastCLine,curNumStrands,nStrands);
    }
    free(lastCLine); lastCLine=NULL;
  }

  totalOrders2 = offset;
  qsort( allPermutations, totalOrders2, sizeof(permutation), &comparePermutations);

  totalSets = CountSets(allPermutations, totalOrders2, nStrands);

  allSets = (multiset*) malloc( (totalSets)*sizeof( multiset));

  int maxSeqLength = FillSets(allSets, allPermutations,
                          totalSets, totalOrders2,
                          nStrands, seqlength) ;

  maxListComplexSize = GetMaxComplexSize(allSets,totalSets);
  nicks = (int*) malloc( maxListComplexSize*sizeof(int) );

  if( !(globalArgs.quiet))
    printf("Permutation generation complete.\n");
  if (globalArgs.listonly) setStart=nSets; // Don't do all sets

  if( nTotalOrders != totalOrders2) {
    printf("Internal error! Total number of permutations is incorrect! %d != %d\n",
           nTotalOrders, totalOrders2);
    exit(1);
  }

  if( F_list != NULL) {
    fclose( F_list);
  }

  if( !(globalArgs.quiet) && !(globalArgs.timeonly))
    printf("Starting partition function calculations.\n");


  // 2006/03/08, Add in a sorting subroutine, to order sets by total # of strands
  // Don't sort with -listonly so we can skip all i < nSets
  if (!globalArgs.listonly) qsort( allSets, nSets, sizeof( multiset), &compareMultisets);

  //allocate memory for pfSeq;
  pfSeq = (char*) malloc( (maxSeqLength + 1) *sizeof(char) );

  if( globalArgs.permsOn) {
    //for( i = 0; i <= totalSets-1; i++)
    //  printPerms( F_perm, i+1, nStrands, &(allSets[i]));


    permPr = (long double **)
      malloc( nStrands*sizeof( long double*));
    for( j = 0; j < nStrands; j++) { //calloc initialize to zero
      permPr[j] = (long double*)
        calloc( seqlength[j], sizeof( long double));
    }
  }

  //estimate total time needed and denominator for progress
  for( i = setStart; i <= totalSets-1; i++) {
    estimatedTime += allSets[i].nPerms * N3C * pow(allSets[i].totalLength, 3);
    progDenom += pow(allSets[i].totalLength,3);
  }
  if( !(globalArgs.quiet) || globalArgs.timeonly) {
    printf("Rough time estimate for calculation: %.2f seconds\n", estimatedTime);
  }

  for( i = setStart; i <= totalSets-1; i++) {
    //initialize complex mfe variables
    allSets[i].mfe = 10000;
    allSets[i].nMfePerms = 0;
    allSets[i].mfePerms = (int *) calloc( 10, sizeof(int));
  }


  status = setStart;
  if (!globalArgs.timeonly) {
    for( i = setStart; i <= totalSets-1; i++) {

      //time( &end);
      if( !(globalArgs.quiet))
        printf("Status: Set %d / %d: nPerms (%d)  %d / %d\n", i+1,
               totalSets,
               allSets[i].nPerms, status+1, totalOrders2);
      status += allSets[i].nPerms;

      allSets[i].pf = 0; //initialize pf

      if( globalArgs.dopairs) {
        //initialize avgBp
        allSets[i].avgBp = (long double **)
          malloc( nStrands*sizeof( long double*));
        for( j = 0; j < nStrands; j++) { //calloc initialize to zero
          allSets[i].avgBp[j] = (long double*)
            calloc( seqlength[j], sizeof( long double));
        }
        //initialize expectation of pairs for complex
        for( j = 0; j < totalStrandsLength*totalStrandsLength; j++)
          cxAvgPairs[ j] = 0;

      }


      currentPerm = allSets[i].perms;
      permId = 1;
      while( currentPerm != NULL) {
        resetNicks( maxListComplexSize, nicks);

        seqCode = (currentPerm->code)[0] - 1;
        strcpy(pfSeq, seqs[ seqCode]);

        //set sequences and nicks
        if( allSets[i].nSeqs >= 2) nicks[0] = seqlength[ seqCode] - 1;

        for( k = 0; k <= allSets[i].nSeqs-2; k++) {

          seqCode = (currentPerm->code)[k+1] - 1;

          strcat( pfSeq, "+");
          strcat( pfSeq, seqs[ seqCode]);

          if( k != allSets[i].nSeqs-2)
            nicks[k+1] = nicks[k]+
            seqlength[ seqCode];

        }

        strncpy( currentPerm->seq, pfSeq,
                allSets[i].totalLength + allSets[i].nSeqs);


        if( globalArgs.dopairs) {
          pairPr = (DBL_TYPE *) malloc( (allSets[i].totalLength+1)*
                                       (allSets[i].totalLength+1)*
                                       sizeof( DBL_TYPE) );
          if( globalArgs.permsOn) {
            for( j = 0; j < nStrands; j++)  //initialize to zero
              for( k = 0; k < seqlength[j]; k++)
              permPr[j][k] = 0;

            for( j = 0; j < totalStrandsLength*totalStrandsLength; j++)
              permAvgPairs[j] = 0; //initialize to zero

          }
        }

        //call library function to compute pseudoknot-free partition function
        tmpLength = strlen( pfSeq);
        convertSeq(pfSeq, seqNum, tmpLength);
        pf =  pfuncFullWithSym( seqNum, 3, globalArgs.parameters,
                                globalArgs.dangles,
                                globalArgs.T,
                                globalArgs.dopairs,
                                currentPerm->symmetryFactor,
                                globalArgs.sodiumconc,
                                globalArgs.magnesiumconc,
                                globalArgs.uselongsalt);

        if( globalArgs.mfe) {
          //if desired, calculate mfe.  This could potentially be exponentially slow
          //if there is a high degree of symmetry.
          mfeValue = mfeFullWithSym( seqNum, tmpLength, &mfes, 3, globalArgs.parameters,
                                     globalArgs.dangles, globalArgs.T,
                                     currentPerm->symmetryFactor,globalArgs.onlyOneMFE,
                                     globalArgs.sodiumconc, globalArgs.magnesiumconc,
                                     globalArgs.uselongsalt);


          if( mfes.nStructs >= 1) {
            mfeValue = mfes.validStructs[0].correctedEnergy;
          }

          if( mfeValue <= allSets[i].mfe) {
            if( mfeValue < allSets[i].mfe) {
              allSets[i].mfe = mfeValue;
              allSets[i].nMfePerms = 1;
              allSets[i].mfePerms[0] = permId;
            }
            else {
              allSets[i].mfePerms[ (allSets[i].nMfePerms)++] = permId;
              if( allSets[i].nMfePerms % 10 == 0) {
                allSets[i].mfePerms =
                  (int *) realloc( allSets[i].mfePerms,
                                  (allSets[i].nMfePerms + 10)*sizeof( int));
              }
            }
          }
        }

        //print permutation info
        if( globalArgs.permsOn && globalArgs.out == 1) {
          if( pf <= 0.0) fprintf(F_ocx, "%% ");
          fprintf(F_ocx, "%d\t%d\t", lastCxId, permId);

          for( j = 0; j <= nStrands - 1; j++) {
            fprintf( F_ocx, "%d\t", allSets[i].code[j]); //strand composition
          }

          if( pf > 0.0) {
            if(!NUPACK_VALIDATE) {
              fprintf( F_ocx, "%.8Le\n",-1*(kB*TEMP_K)*LOG_FUNC( pf) );
            } else {
              fprintf( F_ocx, "%.14Le\n",-1*(kB*TEMP_K)*LOG_FUNC( pf) );
            }
          } else {
            fprintf( F_ocx, "No legal secondary structures!\n");
          }

          // Print results to ocx.mfe file
          if (globalArgs.mfe && pf > 0.0) {
            // Print a comment line for separation
            fprintf(F_permMfe,"\n%% ");
            for (l = 0; l < nPercent; l++) {
              fprintf(F_permMfe,"%%");
            }
            fprintf(F_permMfe," %%\n");

            // Print the record number /* version 3 output */
            fprintf(F_permMfe,"%% %s%d-%s%d\n",composition,lastCxId,ordering,permId);

            // Print the total strand length
            fprintf(F_permMfe,"%d\n",allSets[i].totalLength);
          }


          // Print results to the ocx-ppairs and ocx-epairs files
          if( globalArgs.dopairs && pf > 0.0) {

            // Print a comment line for separation
            fprintf(F_permPr,"\n%% ");
            fprintf(F_permAvg,"\n%% ");
            for (l = 0; l < nPercent; l++) {
              fprintf(F_permPr,"%%");
              fprintf(F_permAvg,"%%");
            }
            fprintf(F_permPr," %%\n");
            fprintf(F_permAvg," %%\n");

            // Print the record number /* version 3 output */
            fprintf(F_permPr,"%% %s%d-%s%d\n",composition,lastCxId,ordering,permId);
            fprintf(F_permAvg,"%% %s%d-%s%d\n",composition,lastCxId,ordering,permId);

            // Print the total strand length
            fprintf(F_permPr,"%d\n",allSets[i].totalLength);
            fprintf(F_permAvg, "%d\n",totalStrandsLength);

          }

          if( globalArgs.dodefect ) {
            if( pf <= 0.0) {
              fprintf(F_defect, "%% ");
            }
            fprintf(F_defect, "%d\t%d\t", lastCxId, permId);

            for( j = 0; j <= nStrands - 1; j++) {
              fprintf( F_defect, "%d\t", allSets[i].code[j]); //strand composition
            }
          }

          permId++;
        }

        if(globalArgs.dodefect ) {
          if(pf > 0) {
            int * cur_struct = mfes.validStructs[0].theStruct;
            DBL_TYPE cur_mfe = mfes.validStructs[0].correctedEnergy;
            DBL_TYPE cur_defect = 0;
            DBL_TYPE norm_defect = 0 ;
            for(q = 0 ; q < mfes.seqlength ; q++) {
              r = cur_struct[q];
              if(r < 0) {
                r = mfes.seqlength;
              }
              pf_qr = (q * (allSets[i].totalLength + 1)) + r;
              cur_defect += 1 - (pairPr[pf_qr]);
              norm_defect = cur_defect / mfes.seqlength;
            }

            if(NUPACK_VALIDATE) {
              fprintf(F_defect, "%.14Le\t",(long double) cur_defect);
              fprintf(F_defect, "%.14Le\t",(long double) norm_defect);
              fprintf(F_defect, "%.14Le\t",(long double) EXP_FUNC(-cur_mfe/(kB*TEMP_K))/pf);
              fprintf(F_defect, "%.14Le\n",(long double) pf);
            } else {
              fprintf(F_defect, "%.6Le\t",(long double) cur_defect);
              fprintf(F_defect, "%.6Le\t",(long double) norm_defect);
              fprintf(F_defect, "%.6Le\t",(long double) EXP_FUNC(-cur_mfe/(kB*TEMP_K))/pf);
              fprintf(F_defect, "%.6Le\n",(long double) pf);
            }
          } else {
            fprintf(F_defect, "No legal secondary structures!\n");
          }
        }

        if( globalArgs.mfe && globalArgs.permsOn && pf > 0.0) {
          printMfesToFile( &mfes, F_permMfe, nicks);
          // Print a comment line for separation
          fprintf(F_permMfe,"%% ");
          for (l = 0; l < nPercent; l++) {
            fprintf(F_permMfe,"%%");
          }
          fprintf(F_permMfe," %%\n");
        }

        if( globalArgs.mfe) {
          clearDnaStructures( &mfes);
        }

        //avgBp[i][j] initialized to zero; average bp for this permutation

        if( globalArgs.dopairs && pf > 0.0) {
          for( q = 0; q <= allSets[i].totalLength-1; q++) {
            //compute the strand type and relative position for each base in the
            //concatenation of sequences
            strandId = (currentPerm->baseCode)[2*q];
            strandPos = (currentPerm->baseCode)[2*q+1];
            indQ = (1+allSets[i].totalLength)*q;
            for( r = 0; r <= allSets[i].totalLength-1; r++) {
              strandId2 = (currentPerm->baseCode)[2*r];
              strandPos2 = (currentPerm->baseCode)[2*r+1];

              pf_qr = indQ + r;

              ppf = pairPr[ pf_qr]*pf;

              allSets[i].avgBp[ strandId][ strandPos] += ppf;

              cxAvgPairs[ indexCx[ strandId][ strandPos]*totalStrandsLength +
                          indexCx[ strandId2][ strandPos2] ] += ppf;

              if( globalArgs.permsOn)  {
                permPr[ strandId][strandPos] += pairPr[ pf_qr];
                permAvgPairs[ indexCx[ strandId][ strandPos]*totalStrandsLength +
                             indexCx[ strandId2][ strandPos2] ] +=
                  pairPr[ pf_qr];

                if( globalArgs.out == 1 && r > q &&
                   pairPr[ pf_qr] >= globalArgs.cutoff) {
                    if(!NUPACK_VALIDATE) {
                     fprintf( F_permPr, "%d\t%d\t%.3Le\n", q+1, r+1, (long double)pairPr[ pf_qr]);
                    } else {
                     fprintf( F_permPr, "%d\t%d\t%.14Le\n", q+1, r+1, (long double)pairPr[ pf_qr]);
                    }
                   }
              }

            }
          }
          if( globalArgs.out == 1 && globalArgs.permsOn && pf > 0.0) {

            //add in unpaired probabilities to F_permPr (.ocx-ppairs)
            for( q = 0; q <= allSets[i].totalLength-1; q++) {
              r = allSets[i].totalLength;
              pf_qr = (1+allSets[i].totalLength)*q+r;
              if(!NUPACK_VALIDATE) {
                fprintf( F_permPr, "%d\t%d\t%.3Le\n", q+1, r+1, (long double)pairPr[ pf_qr]);
              } else {
                fprintf( F_permPr, "%d\t%d\t%.14Le\n", q+1, r+1, (long double)pairPr[ pf_qr]);
              }
            }
            // Print a comment line for separation
            fprintf(F_permPr,"%% ");
            for (l = 0; l < nPercent; l++) {
              fprintf(F_permPr,"%%");
            }
            fprintf(F_permPr," %%\n");

            //store average number of each type of pair (.ocx-epairs)
            for( q = 0; q < totalStrandsLength*totalStrandsLength; q++) {
              row = (q/totalStrandsLength) + 1;
              column = (q%totalStrandsLength) + 1;

              if( permAvgPairs[q] >= globalArgs.cutoff &&
                 row < column ) {
                  if(!NUPACK_VALIDATE) {
                   fprintf( F_permAvg, "%d\t%d\t%.3Le\n", row, column, permAvgPairs[q]);
                  } else {
                   fprintf( F_permAvg, "%d\t%d\t%.14Le\n", row, column, permAvgPairs[q]);
                  }
                 }
            }

          }

          free( pairPr); pairPr = NULL;

          //Also store average unpaired for each base (.ocx-epairs)
          if( globalArgs.permsOn && globalArgs.out == 1 && pf > 0.0) {
            for( q = 0; q < nStrands; q++) {
              for( r = 0; r < seqlength[q]; r++) {
                row = indexCx[q][r] + 1;
                column = totalStrandsLength + 1;
                expectedUnpaired = allSets[i].code[q] - permPr[q][r];
                if (expectedUnpaired >= globalArgs.cutoff) {
                  if(!NUPACK_VALIDATE) {
                  fprintf( F_permAvg, "%d\t%d\t%.3Le\n", row, column, expectedUnpaired);
                } else {
                  fprintf( F_permAvg, "%d\t%d\t%.14Le\n", row, column, expectedUnpaired);
                }
                }
              }
            }
            // Print a comment line for separation
            fprintf(F_permAvg,"%% ");
            for (l = 0; l < nPercent; l++) {
              fprintf(F_permAvg,"%%");
            }
            fprintf(F_permAvg," %%\n");
          }
        }
        currentPerm->pf = pf;
        allSets[i].pf += pf;

        currentPerm = currentPerm->next;

      }


      /* version 3 output */
      //print to .cx on the fly
      if(globalArgs.v3 && globalArgs.out == 1) {
        if( allSets[i].pf <= 0) {
          fprintf( F_cx, "%% ");
        }

        fprintf( F_cx, "%d\t", lastCxId);

        for( j = 0; j <= nStrands - 1; j++) {
          fprintf( F_cx, "%d\t", allSets[i].code[j]); //strand composition
        }

        if( allSets[i].pf > 0) {
          if(!NUPACK_VALIDATE) {
            fprintf( F_cx, "%.8Le", -1*(kB*TEMP_K)*LOG_FUNC( allSets[i].pf));
          } else {
            fprintf( F_cx, "%.14Le", -1*(kB*TEMP_K)*LOG_FUNC( allSets[i].pf));
          }
        } else {
          fprintf( F_cx, "No legal secondary structures!");
        }

        fprintf( F_cx, "\n");
      }

      /* version 3 output */
      if(globalArgs.v3 && globalArgs.out == 1 && globalArgs.dopairs) {

        // Row of percent signs for separation
        fprintf(F_cxAvg,"\n%% ");
        for (l = 0; l < nPercent; l++) {
          fprintf(F_cxAvg,"%%");
        }
        fprintf(F_cxAvg," %%\n");

        // Print the record number
        fprintf(F_cxAvg,"%% complex%d\n",lastCxId);

        // Print N_distinct
        fprintf(F_cxAvg,"%d\n",totalStrandsLength);

      }

      if( globalArgs.out == 1 && globalArgs.permsOn && allSets[i].pf > 0.0) {
        printPerms( F_perm, lastCxId, nStrands, &(allSets[i]));
      }

      if( allSets[i].pf > 0.0) lastCxId++; //this will keep complex Ids consecutive

      if( globalArgs.dopairs && allSets[i].pf > 0.0) {
        for( q = 0; q < nStrands; q++) {
          for( r = 0; r < seqlength[q]; r++) {
            allSets[i].avgBp[q][r] /= allSets[i].pf;
            //divide by total pf to get average
          }
        }

        //print expected number of pairs of each type
        for( q = 0; q < totalStrandsLength*totalStrandsLength; q++) {
          cxAvgPairs[ q ] /= allSets[i].pf;
          row = (q / totalStrandsLength)+1; //1 to totalStrandsLength
          column = (q % totalStrandsLength)+1; //1 to totalStrandsLength
          /* version 3 output */
          if(globalArgs.v3) {
            if( cxAvgPairs[q] >= globalArgs.cutoff && row < column ) {
              if(!NUPACK_VALIDATE) {
                fprintf( F_cxAvg, "%d\t%d\t%.3Le\n", row, column, cxAvgPairs[q]);
              } else {
                fprintf( F_cxAvg, "%d\t%d\t%.14Le\n", row, column, cxAvgPairs[q]);
              }
            }
          }
        }

      }

      //print expected number of unpaired bases of each type
      if( globalArgs.dopairs &&  globalArgs.out == 1 && allSets[i].pf > 0.0) {
        // Avg Base Pairs
        for( q = 0; q < nStrands; q++) {
          for( r = 0; r < seqlength[q]; r++) {
            row = indexCx[q][r] + 1;
            column = totalStrandsLength + 1;
            expectedUnpaired = allSets[i].code[q] - allSets[i].avgBp[q][r];
            /* version 3 output */
            if(globalArgs.v3) {
              if( expectedUnpaired >= globalArgs.cutoff) {
                if(!NUPACK_VALIDATE) {
                  fprintf( F_cxAvg, "%d\t%d\t%.3Le\n", row, column, expectedUnpaired);
                } else {
                  fprintf( F_cxAvg, "%d\t%d\t%.14Le\n", row, column, expectedUnpaired);
                }
              }
            }
          }
        }
        /* version 3 output */
        // Row of percent signs for separation
        if(globalArgs.v3) {
          fprintf(F_cxAvg,"%% ");
          for (l = 0; l < nPercent; l++) {
            fprintf(F_cxAvg,"%%");
          }
          fprintf(F_cxAvg," %%\n");
        }
      }

      // Assess progress and write to a file
      if (globalArgs.progress) {
        prog += pow(allSets[i].totalLength,3);
        F_prog = fopen(progName,"w");
        if(!NUPACK_VALIDATE) {
          fprintf(F_prog,"%.3f\n\n",prog/progDenom);
        } else {
          fprintf(F_prog,"%.14f\n\n",prog/progDenom);
        }
        fclose(F_prog);
      }

    }
  }

  if( globalArgs.permsOn && !globalArgs.timeonly) {

    for( j = 0; j < nStrands; j++) { //free
      free( permPr[j]);
      permPr[j] = NULL;
    }
    free( permPr);
    permPr = NULL;
  }


  if( !(globalArgs.quiet)) {
    //debug by printing
    //print subsets, perms
    //time( &end);
    //printf("Final Time: %.0f\n", difftime( end,start));

    if (globalArgs.debug){
      for( i = setStart; i<=totalSets-1; i++) {
        for( j = 0; j <= nStrands - 1; j++) {
          printf("%d", allSets[i].code[j]);
        }
        if(!NUPACK_VALIDATE) {
          printf(" %.8Le\n", allSets[i].pf);
        } else {
          printf(" %.14Le\n", allSets[i].pf);
        }

        currentPerm = allSets[i].perms;
        while( currentPerm != NULL) {
          for( k = 0; k <= allSets[i].nSeqs - 1; k++) {
            printf("%d", (currentPerm->code)[k]);
          }
          if(!NUPACK_VALIDATE) {
            printf(" - %s %.8Le\n", currentPerm->seq,
                 currentPerm->pf);
          } else {
            printf(" - %s %.8Le\n", currentPerm->seq,
                 currentPerm->pf);
          }
          currentPerm = currentPerm->next;
        }
        printf("\n");
      }
    }
  }

  if( globalArgs.out == 1) {
    /* version 3 output */
    if(globalArgs.v3) {
      fclose( F_cx);
    }
    if( globalArgs.mfe) {
      fclose( F_permMfe);
    }
    if( globalArgs.permsOn) {
      fclose( F_perm);
      fclose( F_ocx);
    }
    if( globalArgs.dopairs) {
      /* version 3 output */
      if(globalArgs.v3) {
        fclose( F_cxAvg);
      }
      if( globalArgs.permsOn) {
        fclose( F_permPr);
        fclose( F_permAvg);
      }
    }
    if( globalArgs.dodefect) {
      fclose(F_defect);
    }
  }

  free( nicks); nicks = NULL;

  for( i = 0; i<=nStrands-1; i++) {
    free( seqs[i]); seqs[i] = NULL;

    if( globalArgs.dopairs) {
      free( indexCx[i]);
      indexCx[i] = NULL;
    }
  }

  free( seqlength); seqlength = NULL;
  free( seqs); seqs = NULL;

  if( globalArgs.dopairs) {
    free( indexCx); indexCx = NULL;
    free( cxAvgPairs); cxAvgPairs = NULL;
    free( permAvgPairs); permAvgPairs = NULL;
  }


  for( i = setStart; i <= totalSets-1; i++) {

    if( globalArgs.dopairs && !globalArgs.timeonly) {
      for(j = 0; j <= nStrands -1 ; j++) {
        free( allSets[i].avgBp[j]); allSets[i].avgBp[j] = NULL;
      }
      free( allSets[i].avgBp);
    }

    free( allSets[i].code); allSets[i].code = NULL;

    free( allSets[i].mfePerms);
    allSets[i].mfePerms = NULL;

    currentPerm = allSets[i].perms;
    while( currentPerm != NULL) {
      free( currentPerm->code); currentPerm->code = NULL;
      free( currentPerm->baseCode); currentPerm->code = NULL;
      free( currentPerm->seq); currentPerm->seq = NULL;
      free( currentPerm->strand_sums); currentPerm->strand_sums = NULL;
      tmpPerm = currentPerm;
      currentPerm = currentPerm->next;
    }
  }
  free( allPermutations); allPermutations = NULL;
  free( allSets); allSets = NULL;

  //free( seqlength); seqlength = NULL;
  free( pfSeq); pfSeq = NULL;

  if( !(globalArgs.quiet) && !(globalArgs.timeonly)) {
    printf("Total number of terms calculated: %d (%d)\n",
           totalOrders2, nTotalOrders);
    curtime = time(NULL); //current time
    loctime = localtime( &curtime);
    printf( "Calculation finished on: %s\n",
           asctime( loctime));
  }

  if( globalArgs.echo) {
    printf( "%s\n", filePrefix);
  }


  return 0;
}


