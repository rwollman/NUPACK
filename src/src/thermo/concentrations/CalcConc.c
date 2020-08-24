/*
  CalcConc.c is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.
  Coded by: Justin Bois 9/2006

  CALCCONC.C

  For use with Concentrations.c.

  Computes the equilibrium mole fractions of the products of
  aggregation reactions in dilute solution, given the identity of the
  aggregates (called "complexes") and their repsective free energies.

  This program solves the problem presented in Dirks, et al.,
  "Thermodynamic analysis of interacting nucleic acid strands", SIAM
  Review, in press.  All variable names are chosen to match those in
  that paper.

  This file contains CalcConc, the main function that computes the
  equilibrium mole fractions of the complexes, and the auxillary
  functions it calls.  Some of the functions it calls are standard
  utility functions, such as functions to sum entries in an array,
  etc.  These are included in utils.c.

  The trust region algorithm for solving the dual problem is that in
  Nocedal and Wright, Numerical Optimization, 1999, page 68, with the
  dogleg method on page 71.  There are some inherent precision issues.
  For most systems involving nucleic acids, this isn't a problem, but
  for some problematic cases, adjustments are made.  For some initial
  conditions, these precision issues cannot be overcome and a new
  initial condition must be generated.  This is done by randomly
  perturbing the standard initial condition (see comments in the
  function getInitialCondition, below), and re-running the trust
  region optimization.

  The inputs are as follows:
  A: A 2-D array; A[i][j] is the number of monomers of type i in complex j
  G: An array containing the corresponding complex free energies in units of kT.
     G[j] corresponds to the entries A[..][j].
  x0: Initial mole fractions of the unit-size complexes as mole fractions.
  numSS: The number of single-species (monomers) in the system.
  numTotal: The total number of complexes.
  maxIters: Maximum number of interations allowed in trust region method.
  tol: The tolerance for convergence.  The absolute tolerance is tol*(mininium 
       single-species initial mole fraction)
  deltaBar: The maximum step size allowed in the trust region method
  eta: The value for eta in the trust region method, 0 < eta < 1/4
  outputFile: The file to which the output is written.
  kT: kT in kcal/mol.
  MaxNoStep: The maximum number of iterations without a step being taken before
             the initial conditions are regenerated.
  MaxTrial: The maximum number of initial conditions to be tried.
  PerturbScale: The multiplier on the random perturbations to the initial conditions
                as new ones are generated.
  quiet: = 1 for no printing of messages (except error messages) to screen.  
  WriteLogFile: = 1 if log file is to be written
  logFile: file for printing information about the run.
  MolesWaterPerLiter: Number of moles of water per liter
*/


#include "CalcConc.h"
#include "constants.h"

/* ******************************************************************************** */

/* ******************************************************************************** */
/*                              BEGIN CALCCONC FUNCTION                             */
/* ******************************************************************************** */
int CalcConc(double *x, int **A, double *G, double *x0, int numSS, int numTotal, 
         int MaxIters, double tol, double deltaBar, double eta, double kT, 
         int MaxNoStep, int MaxTrial, double PerturbScale, int quiet, 
         int WriteLogFile, char *logFile, double MolesWaterPerLiter,
         unsigned long seed) {
  /*
    Computes the equilbrium mole fractions of species in dilute
    solution using a trust region algorithm on the dual problem.
    Discussion of the method is in Dirks, et al., Thermodynamic
    analysis of interacting nucleic acid strands, SIAM Review, (2006),
    in press.  The trust region algorithm for solving the dual problem
    is that in Nocedal and Wright, Numerical Optimization, 1999, page
    68, with the dogleg method on page 71.

    Returns 1 if converged and 0 otherwise.
  */

  int i,j; // Counters, i is over single-species and j is over all complexes
  int iters; // Number of iterations
  double *AbsTol; // The absolute tolerance on all values of gradient
  double rho; // Ratio of actual to predicted reduction in trust region method
  double delta; // Radius of trust region
  double *Grad; // The gradient of -g(lambda)
  double *lambda; // Lagrange multipliers (dual variables),x[j] = Q[j]*exp(lambda[j])
                  // for j \in \Psi^0
  double *p; // The step we take toward minimization
  double **Hes; // The Hessian
  double FreeEnergy; // The free energy of the solution
  unsigned long rand_seed = 0; // Random number seed
  int nNoStep; // Number of iterations without taking a step
  int nTrial; // Number of times we've perturbed lambda
  int **AT; // Transpose of A
  int RunStats[6]; // Statistics on results from getSearchDir (see comments below)
  FILE *fplog; // Log file

  // Initialize iters just so compiler doesn't give a warning when optimization is on
  iters = 0;

  // Allocate memory
  AT = (int **) malloc(numTotal * sizeof(int *));
  for (j = 0; j < numTotal; j++) {
    AT[j] = (int *) malloc(numSS * sizeof(int));
  }
  Hes = (double **) malloc(numSS * sizeof(double *));
  for (i = 0; i < numSS; i++) {
    Hes[i] = (double *) malloc(numSS * sizeof(double));
  }

  AbsTol = (double *) malloc(numSS * sizeof(double));
  Grad = (double *) malloc(numSS * sizeof(double));
  lambda = (double *) malloc(numSS * sizeof(double));
  p = (double *) malloc(numSS * sizeof(double));

  // The absolute tolerance is a percentage of the entries in x0
  for (i = 0; i < numSS; i++) {
    AbsTol[i] = tol * x0[i];
  }

  // Compute AT (transpose of A), useful to have around.
  IntTranspose(AT,A,numSS,numTotal);

  nTrial = 0;
  for (i = 0; i < numSS; i++) {
    Grad[i] = AbsTol[i] + 1.0; // Initialize just to get started.
  }
  while (CheckTol(Grad,AbsTol,numSS) == 0 && nTrial < MaxTrial) {

    if (nTrial == 1) {
      // Seed the random number generator if necessary
      rand_seed = GetRandSeed(seed);
      init_genrand(rand_seed);
    }

    // Set initial guess
    getInitialGuess(x0,lambda,G,AT,A,numSS,numTotal,PerturbScale,rand_seed);

    // Calculate the counts of the species based on lambda
    if (getx(x,lambda,G,AT,numSS,numTotal) == 0) { // Should be fine; checked prev.
      if (quiet == 0) {
    printf("Overflow error in calcution of mole fractions.\n\n");
    printf("Exiting....\n");
      }
      exit(ERR_OVERFLOW);
    }

    // Calculate the gradient
    getGrad(Grad,x0,x,A,numSS,numTotal);
    
    // Initialize delta to be just less than deltaBar
    delta = 0.99 * deltaBar;
    
    // Initializations
    iters = 0;
    nNoStep = 0;
    RunStats[0] = 0; // Number of pure Newton steps (didn't hit trust region boundary)
    RunStats[1] = 0; // Number of pure Cauchy steps (hit trust region boundary)
    RunStats[2] = 0; // Number of dogleg steps (part Newton and part Cauchy)
    RunStats[3] = 0; // Number of steps with Cholesky failure forcing Cauchy step
    RunStats[4] = 0; // Number of steps with irrelevant Cholesky failures
    RunStats[5] = 0; // Number of failed dogleg calculations
    
    // Run trust region with these initial conditions
    while (iters < MaxIters && CheckTol(Grad,AbsTol,numSS) == 0 
      && nNoStep < MaxNoStep) {
      

      // Compute the Hessian (symmetric, positive, positive definite)
      getHes(Hes,x,A,numSS,numTotal);

      // Solve for the search direction
      (RunStats[getSearchDir(p,Grad,Hes,delta,numSS) - 1])++;

      // Calculate rho, ratio of actual to predicted reduction
      rho = getRho(lambda,p,Grad,x,Hes,x0,G,AT,numSS,numTotal);
      
      // Adjust delta and make step based on rho
      if (rho < 0.25) {
        delta /= 4.0;
      }
      else if (rho > 0.75 && fabs(norm(p,numSS) - delta) < NUM_PRECISION) {
        delta = min2(2.0*delta,deltaBar);
      }
      if (rho > eta) {    
        for (i = 0; i < numSS; i++) {
          lambda[i] += p[i];
        }
        nNoStep = 0;    
      }
      else {
        nNoStep++;
      }
      
      // Calculate the mole fractions of the complexes based on lambda
      if (getx(x,lambda,G,AT,numSS,numTotal) == 0) {// Should be fine;checked prev.
        if (quiet == 0) {
          printf("Overflow error in calcution of mole fractions.\n\n");
          printf("Exiting....\n");
        }
        exit(ERR_OVERFLOW);
      }
      
      // Calculate the gradient
      getGrad(Grad,x0,x,A,numSS,numTotal);
      
      // Advance the iterations count
      iters++;
    }

    // Advance the number of perturbations we've tried
    nTrial++;
  }

  // Compute the free energy
  FreeEnergy = 0;
  // First the reference free energy
  for (i = 0; i < numSS; i++) {
    FreeEnergy += x0[i]*(1.0 - log(x0[i]));
  }
  // Now the free energy
  for (j = 0; j < numTotal; j++) {
    if (x[j] > 0) {
      FreeEnergy += x[j]*(log(x[j]) + G[j] - 1.0);
    }
  }
  // Convert to kcal/liter of solution
  FreeEnergy *= kT*MolesWaterPerLiter;

  /* **************** WRITE OUT RESULTS ********************************* */
  if ( nTrial == MaxTrial  && quiet == 0) {
    printf("\n\n   TRUST REGION METHOD DID NOT CONVERGE DUE TO PRECISION ISSUES\n\n");
  }

  // Report errors in conservation of mass to screen
  if (quiet == 0) {
    // Print out values of the gradient, which is the error in cons. of mass
    printf("Error in conservation of mass:\n");
    for (i = 0; i < numSS; i++) {
      printf("   %8.6e Molar\n",Grad[i]*MolesWaterPerLiter);
    }
    printf("\n");

    // Print out the free energy of the solution
    printf("Free energy = %8.6e kcal/litre of solution\n",FreeEnergy);
   }

  // Write out details of calculation to outfile
  if (WriteLogFile) {
    if ((fplog = fopen(logFile,"a")) == NULL) {
      if (quiet == 0) {
        printf("Error opening %s.\n\nExiting....\n",logFile);
      }
      exit(ERR_LOG);
    }
    if (nTrial == MaxTrial) {
      fprintf(fplog,"TRUST REGION DID NOT CONVERGE DUE TO PRECISION ISSUES\n\n");
    }
    fprintf(fplog,"   --Trust region results:\n");
    fprintf(fplog,"       No. of initial conditions tried: %d\n",nTrial);
    fprintf(fplog,"       Results from succesful trial:\n");
    fprintf(fplog,"         No. of iterations: %d\n",iters);
    fprintf(fplog,"         No. of Newton steps: %d\n",RunStats[0]);
    fprintf(fplog,"         No. of Cauchy steps: %d\n",RunStats[1]);
    fprintf(fplog,"         No. of dogleg steps: %d\n",RunStats[2]);
    fprintf(fplog,"         No. of Cholesky failures resulting in Cauchy steps: %d\n"
      ,RunStats[3]);
    fprintf(fplog,"         No. of inconsequential Cholesky failures: %d\n",
      RunStats[4]);
    fprintf(fplog,"         No. of dogleg failures: %d\n",RunStats[5]);
    fprintf(fplog,"         Error in conservation of mass (units of molarity):\n");
    fprintf(fplog,"              Error\tTolerance\n");
    for (i = 0; i < numSS; i++) {
      fprintf(fplog, "           %.14e\t%.14e\n",Grad[i]*MolesWaterPerLiter,
          AbsTol[i]*MolesWaterPerLiter);
    }
    fprintf(fplog,"   --Free energy of solution = %.14e kcal/litre of solution\n",
        FreeEnergy);
    fclose(fplog);
  }
  /* **************** END OF WRITING OUT RESULTS **************************** */

   // Free memory
  for (j = 0; j < numTotal; j++) {
    free(AT[j]);
  }
  for (i = 0; i < numSS; i++) {
    free(Hes[i]);
  }
  free(AbsTol);
  free(AT);
  free(Hes);
  free(Grad);
  free(p);
  free(lambda);

  // Return convergence
  if (nTrial == MaxTrial) {
    return 0;
  }
  else {
    return 1;
  }

}
/* ******************************************************************************** */
/*                              END CALCCONC FUNCTION                               */
/* ******************************************************************************** */


/* ******************************************************************************** */
void getInitialGuess(double *x0, double *lambda, double *G, int **AT, int **A, 
             int numSS, int numTotal, double PerturbScale, 
             unsigned long rand_seed) {
  /*
    Calculates an initial guess for lambda such that the maximum mole
    fraction calculated will not give an overflow error and the
    objective function $-g(\lambda)$ will be positive.  It is best to
    have a positive objective function because when the objective
    function is negative, it tends to be very close to zero and there
    are precision issues.

    We assume all the lambda's have the same value in the initial condition.
    We compute the maximal lambda such that all mole fractions of all complexes
    are below some maximum.
  */

  int i,j; // Counters
  double MaxLogx; // maximum log of the mole fraction allowed
  double LambdaVal; // Possible values of lambda s.t. conc is exp(MaxLogx).
  double NewLambdaVal; // Same as LambdaVal
  double tG;

  MaxLogx = 1.0;  // Maximum mole fraction is ~3

  LambdaVal = (MaxLogx + G[0]) / sumint(AT[0],numSS);
  for (j = 1; j < numTotal; j++) {
    NewLambdaVal = (MaxLogx + G[j]) / sumint(AT[j],numSS);
    if (NewLambdaVal < LambdaVal) {
      LambdaVal = NewLambdaVal;
    }
  }

  for (i = 0; i < numSS; i++) {
    lambda[i] = LambdaVal;
  }

  // Perturb Lambda if desired
  if (rand_seed != 0) {
    PerturbLambda(lambda,PerturbScale,G,AT,numSS,numTotal);
  }

  // If we already know concentration (ss species is inert), set lambda
  for (i = 0; i < numSS; i++) {
    if (sumint(A[i],numTotal) == 1) {
      tG = G[FindNonZero(A[i],numTotal)];
      lambda[i] = log(x0[i]) + tG;
    }
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int getx(double *x, double *lambda, double *G, int **AT, int numSS, int numTotal) {
  /* 
     Calculates the mole fractions of all species from lambda, G, and
     A.  Returns 1 if the calculation was ok and 0 if there will be an
     overflow error.
  */

  int j; // Counter
  double logx; // log of the mole fraction

  for (j = 0; j < numTotal; j++) {
    logx = -G[j] + didot(lambda,AT[j],numSS);
    if (logx > MAXLOGX) { // Will have an overflow error
      return 0;
    }
    x[j] = exp(logx);
  }

  // No overflow errors
  return 1;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void getGrad(double *Grad, double *x0, double *x, int **A, int numSS, int numTotal) {
  /*
    Calculates the gradient of -g(\lambda), the dual function for
    which we're trying to find the minimum.
  */

  int i; // Counter

  for (i = 0; i < numSS; i++) {
    Grad[i] = -x0[i] + didot(x,A[i],numTotal);
  }

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void getHes(double **Hes, double *x, int **A, int numSS, int numTotal) {
  /*
    Calculates the Hessian, Hes.  Must be preallocated and of size
    numSS by numSS.  We only need to calculate the upper triangle of
    Hes since it's symmetric.  We can then fill out the lower triangle.

    Avec is a vector consisting of elements of one row of A multiplied by 
    those of another.
  */


  int m,n,j; // Counters
  double *Avec;

  Avec = (double *) malloc(numTotal * sizeof(double));

  for (n = 0; n < numSS; n++) {
    for (m = 0; m <= n; m++) {
      for (j = 0; j < numTotal; j++) {
    Avec[j] = ((double) A[m][j]) * ((double) A[n][j]);
      }
      Hes[m][n] = dot(x,Avec,numTotal);
    }
  }

  // Fill out the lower entries in the Hessian (needed for matrix mults)
  for(m = 1; m < numSS; m++) {
    for (n = 0; n < m; n++) {
      Hes[m][n] = Hes[n][m];
    }
  }

  free(Avec);

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int getSearchDir(double *p, double *Grad, double **Hes, double delta, int numSS) {
  /*
    Computes the search direction using the dogleg method (Nocedal and Wright,
    page 71).  Notation is consistent with that in this reference.

    Due to the construction of the problem, the minimization routine
    to find tau can be solved exactly by solving a quadratic.  There
    can be precision issues with this, so this is checked.  If the
    argument of the square root in the quadratic formula is negative,
    there must be a precision error, as such a situation is not
    possible in the construction of the problem.

    Returns:
      1 if step was a pure Newton step (didn't hit trust region boundary)
      2 if the step was purely Cauchy in nature (hit trust region boundary)
      3 if the step was a dogleg step (part Newton and part Cauchy)
      4 if Cholesky decomposition failed and we had to take a Cauchy step
      5 if Cholesky decompostion failed but we would've taken Cauchy step anyways
      6 if the dogleg calculation failed (should never happen)
  */

  int i,j; // counters
  double a,b,c,sgnb; // Constants used in quadratic formula
  double q,x1,x2; // results from quadratic formula
  double tau; // Multipler in Newtonstep
  double *pB; // Unconstrained minimizer (the regular Newton step)
  double *pU; // Minimizer along steepest descent direction
  double pB2; // pj^2
  double pU2; // pu^2
  double pBpU; // pj . pc
  double *CholDiag; // Diagonal from Cholesky decomposition
  double **HesCopy; // Copy of the upper triangle of the Hessian (don't want to mess
                   // with actual Hessian).
  int CholSuccess; // Whether of not Cholesky decomposition is successful
  double *HGrad; // Hessian dotted with the gradient
  double pUcoeff;  // The coefficient on the Cauchy step
  double delta2; // delta^2
  double mag1;   // temporary variables for vector magnitudes
  double mag2;

  // Initialize pB2 just so compiler doesn't give a warning when optimization is on
  pB2 = 0.0;

  delta2 = pow(delta,2);

  /* ********** Compute the Newton step ************** */
  // Memory allocation for computation of  Newton step
  pB = (double *) malloc(numSS * sizeof(double));
  CholDiag = (double *) malloc(numSS * sizeof(double));
  HesCopy = (double **) malloc(numSS * sizeof(double *));
  for (j = 0; j < numSS; j++) {
    HesCopy[j] = (double *) malloc(numSS * sizeof(double));
  }
 
  // Make a copy of the Hessian because the Cholesky decomposition messes
  // with its entries.  We only have to copy the upper diagonal.
  for (j = 0; j < numSS; j++) {
    for (i = j; i < numSS; i++) {
      HesCopy[i][j] = Hes[i][j];
    }
  }

  CholSuccess = choleskyDecomposition(HesCopy,numSS);

  if (CholSuccess == 1) {
    choleskySolve(HesCopy,numSS,Grad,pB);

    // Free memory from Cholesky computation
    free(CholDiag);
    for (i = 0; i < numSS; i++) {
      free(HesCopy[i]);
    }
    free(HesCopy);

    // Newton step is -H^{-1} Grad
    for (i = 0; i < numSS; i++) {
      pB[i] *= -1.0;
    }

    // If Newton is in trust region, take it
    pB2 = dot(pB,pB,numSS);
    if (pB2 <= delta2) {
      for (i = 0; i < numSS; i++) {
    p[i] = pB[i];
      }
      free(pB);
      return 1; // Signifies we took a pure Newton step
    }
  }
  else { // Free memory from failed Newton step computation
    free(CholDiag);
    for (i = 0; i < numSS; i++) {
      free(HesCopy[i]);
    }
    free(HesCopy);
  }
  /* ************************************************* */


  /* ********** Compute the Cauchy step ************** */
  // Allocate necessary arrays
  HGrad = (double *) malloc(numSS * sizeof(double));
  pU = (double *) malloc(numSS * sizeof(double));

  // The direction of the Cauchy step
  for (i = 0; i < numSS; i++) {
    pU[i] =  -Grad[i];
  }

  // prefactor for the Cauchy step
  MatrixVectorMult(HGrad,Hes,Grad,numSS);
  // Should this be sqrt too?

  mag1 = dot(Grad,Grad,numSS);
  mag2 = dot(Grad,HGrad,numSS);
  pUcoeff = mag1 / mag2;

  for (i = 0; i < numSS; i++) {
    pU[i] = pUcoeff * pU[i];
  }
  free(HGrad); // Don't need this any more

  pU2 = dot(pU,pU,numSS);

  if (pU2 >= delta2) { // In this case we just take the Cauchy step, 0 < tau <= 1
    tau = sqrt(delta2/pU2);
    for (i = 0; i < numSS; i++) {
      p[i] = tau*pU[i];
    }
    free(pU);
    free(pB);
    if (CholSuccess != 1) {
      return 5; // Signifies Cholesky failure, but doesn't matter, would take Cauchy
    }           // regardless
    else {
      return 2; // Signifies that we just took the Cauchy step
    }
  }

  if (CholSuccess != 1) { // We failed computing Newton step and have to take Cauchy
    for (i = 0; i < numSS; i++) {
      p[i] = pU[i];
    }
    free(pU);
    free(pB);
    return 4; // Signifies Cholesky failure and we just took the Cauchy step
  }
  /* ************************************************* */


  /* ************ Take the dogleg step *************** */
  pBpU = dot(pB,pU,numSS); // Need this for dogleg calculation
  
  // Constants for quadratic formula for solving ||pU + (alpha)(pB-pU)||^2 = delta2
  a = pB2 + pU2 - 2.0*pBpU;
  b = 2*(pBpU - pU2);
  c = pU2 - delta2;
  sgnb = 1;

  if(b < 0) {
    sgnb = -1;
  }

  q = -0.5 * (b + sgnb * sqrt(b*b - 4*a*c));
  x1 = q / a;
  x2 = c / q;

  // x2 should be the positive root, x1 should be the negative root.
  if(x2 >= 0 && x2 <= 1.0) {
    for(i = 0; i < numSS; i++) {
      p[i] = pU[i] + x2 * (pB[i] - pU[i]);
    }
    free(pU);
    free(pB);
    return 3; // Signifies we took a dogleg step
  } else if(x1 >= 0 && x1 <= 1.0) {
    for(i = 0; i < numSS; i++) {
      p[i] = pU[i] + x1 * (pB[i] - pU[i]);
    }
    free(pU);
    free(pB);
    return 3;
  } else {
    for (i = 0; i < numSS; i++) {
      p[i] = pU[i];
    }
    free(pU);
    free(pB);
    return 6; // Signifies no root satisfies the dogleg step and we took a Cauchy step
  }
}
/* ******************************************************************************** */


/* ******************************************************************************** */
double getRho(double *lambda, double *p, double *Grad, double *x, double **Hes, 
          double *x0, double *G, int **AT, int numSS, int numTotal) {
  /*
    Calculates rho based on equations 4.4 and 4.1 of Nocedal and
    Wright.  This is the ratio of the actual correction based on the
    stepping a length delta along the search direction to the
    predicted correction of taking the same step.

    Function returns -1 if there is an overflow error in the
    calculation of rho.
  */

  int i; // Counter
  double rho; // That which we return
  double *newlambda; // y after the step
  double *newx; // x after the Newton step
  double *Hp; // the vector Hes*p
  double negh; // -h(lambda)
  double NewNegh; // New one
  double pHp; // p * H * p


  // Allocate memory for arrays used in the calculation.
  newlambda = (double *) malloc(numSS * sizeof(double));
  newx = (double *) malloc(numTotal * sizeof(double));
  Hp = (double *) malloc(numSS * sizeof(double));

  negh = sum(x,numTotal);
  negh -= dot(lambda,x0,numSS);

  // Calculate the new lambda
  for (i = 0; i < numSS; i++) {
    newlambda[i] = lambda[i] + p[i];
  }

  // Calculate the counts of the species based on new lambda
  if (getx(newx,newlambda,G,AT,numSS,numTotal)) {
    NewNegh = sum(newx,numTotal);
    NewNegh -= dot(newlambda,x0,numSS);

    MatrixVectorMult(Hp,Hes,p,numSS);
    pHp = dot(p,Hp,numSS);

    rho = (negh - NewNegh) / (-dot(Grad,p,numSS) - pHp/2.0);

  }
  else {
    // Since the denominator is always positive, an overflow error in the
    // calculation of the mole fractions for the "new lambda" implies
    // a negative denominator, and therefore the value of rho is negative.
  
    rho = -1.0;
  }

  free(newlambda);
  free(newx);
  free(Hp);

  return rho;

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void getCauchyPoint(double *CauchyPoint, double **Hes, double *Grad, double delta, 
            int numSS) {
  /*
   Computes the Cauchy point using the formulas 4.7 and 4.8 in Nocedal
   and Wright, Numerical Optimization (1999), page 70.

   Note that because the Hessian is symmetric that it is equal to its
   transpose.
  */

  int i; // Counters
  double coeff; // Coefficient on gradient for Cauchy point
  double tau; // Multiplier used in Cauchy point calculation
  double normGrad; // ||Grad||
  double *Hgrad; // Hessian dotted with the gradient
  double numerator;   // Used to hold temporary variables
  double denominator;
  
  Hgrad = (double *) malloc (numSS * sizeof(double));

  normGrad = norm(Grad,numSS);
  MatrixVectorMult(Hgrad,Hes,Grad,numSS);
  numerator = pow(normGrad,3);
  denominator = delta * dot(Grad,Hgrad,numSS);
  tau = min2(numerator/denominator , 1.0);
  coeff = -tau * delta / normGrad;

  for (i = 0; i < numSS; i++) {
    CauchyPoint[i] = coeff*Grad[i];
  }
   
  free(Hgrad);

}
/* ******************************************************************************** */


/* ******************************************************************************** */
void PerturbLambda(double *lambda, double PerturbScale, double *G, int **AT, 
           int numSS, int numTotal) {
  /*
    Perturbs the values of Lambda in case the trust region has shrunk to be 
    very small.

    Adds PerturbScale*random number to each entry in Lambda
  */

  int i; // Counter
  double *dummyx; // A dummy mole fraction vector for checking overflow
  double *newlambda; // The new perturbed lambda
  int xOK; // = 1 is there is no overflow error induced in the concentrations

  dummyx = (double *) malloc(numTotal * sizeof(double));
  newlambda = (double *) malloc(numSS * sizeof(double));

  xOK = 0;
  while (xOK == 0) {
    for (i = 0; i < numSS; i++) {
      newlambda[i] = lambda[i] + PerturbScale * 2.0*(genrand_real1() - 0.5);
    }
    xOK =  getx(dummyx,newlambda,G,AT,numSS,numTotal);
    PerturbScale /= 2.0; // Reduce scale in order to try not to have overflow problems
  }

  // Copy new perturbed lambda to the value of lambda
  for (i = 0; i < numSS; i++) {
    lambda[i] = newlambda[i];
  }
  
  free(dummyx);
  free(newlambda);

}
/* ******************************************************************************** */


/* ******************************************************************************** */
int CheckTol(double *Grad, double *AbsTol, int numSS) {
  /*
    Check the to see if the entries in the gradient satisfy the
    tolerance.  Returns 1 if they all do and 0 otherwise.
  */


  int i; // Counter

  for (i = 0; i < numSS; i++) {
    if (fabs(Grad[i]) > AbsTol[i]) {
      return 0;
    }
  }

  return 1;

}
/* ******************************************************************************** */


