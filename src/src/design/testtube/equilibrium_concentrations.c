/*
  equilibrium_concentrations.c

  Based on:
  CalcConc.c, part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.

  Computes the equilibrium mole fractions of the products of
  aggregation reactions in dilute solution, given the identity of the
  aggregates (called "compounds") and their repsective free energies.

  This file contains calculate_concentrations, the main function that
  computes the equilibrium mole fractions of the compounds, and the
  auxillary functions it calls.  Some of the functions it calls are
  standard utility functions, such as functions to sum entries in an
  array, etc.  These are included in utils.c.

  The trust region algorithm for solving the dual problem is that in
  Nocedal and Wright, Numerical Optimization, 1999, page 68, with the
  dogleg method on page 71.  There are some inherent precision issues.
  For most systems in biochemistry, this isn't a problem, but for some
  problematic cases, adjustments are made.  For some initial
  conditions, these precision issues cannot be overcome and a new
  initial condition must be generated.  This is done by randomly
  perturbing the standard initial condition (see comments in the
  function get_initial_condition, below), and re-running the trust
  region optimization.

  The inputs are as follows:
  A : A 2-D array; A[i][j] is the number of monomers of type i in complex j
  G : Array containing the corresponding compound free energies in units of kT.
     G[j] corresponds to the entries A[..][j].
  x0 : Initial mole fractions of the unit-size complexes as mole fractions.
  n_constraints : The number of single-species (monomers) in the system.
  n_compounds : The total number of complexes.
  max_iters : Maximum number of interations allowed in trust region method.
  tol : The tolerance for convergence.  The absolute tolerance is tol*(mininium 
       single-species initial mole fraction)
  delta_bar : The maximum step size allowed in the trust region method
  eta : The value for eta in the trust region method, 0 < eta < 1/4
  min_delta : The minimal trust region radius allowed before quitting.
  max_trial : The maximum number of initial conditions to be tried.
  perturb_scale : The multiplier on the random perturbations to the initial 
        conditions as new ones are generated.
  quiet : = 1 for no printing of messages (except error messages) to screen.  
  write_log_file : = 1 if log file is to be written
  log_file : file for printing information about the run.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // Must compile with -lm option
#include <float.h>
#include "pathway_debug.h"
#include "equilibrium_concentrations.h" // Header file for concentrations

/* ************************************************************************** */

/* ************************************************************************** */
/*                BEGIN CALC_CONC_FROM_FREE_ENERGIES                          */
/* ************************************************************************** */
int calc_conc_from_free_energies(double *x, double *A, double *G, double *x0, 
                                 int n_particles, int n_compounds, 
                                 int n_points, 
                                 int max_iters, double tol, double delta_bar, 
                                 double eta, double min_delta, 
                                 int max_trial, double perturb_scale, 
                                 int quiet, int write_log_file, char *log_file, 
                                 unsigned long seed, run_stats_t *run_stats) {
  /*
    Computes the equilbrium mole fractions of species in dilute
    solution using a trust region algorithm on the dual problem.
    Discussion of the method is in Dirks, et al., Thermodynamic
    analysis of interacting nucleic acid strands, SIAM Review, (2006),
    in press.  The trust region algorithm for solving the dual problem
    is that in Nocedal and Wright, Numerical Optimization, 1999, page
    68, with the dogleg method on page 71.

    Return codes:
    0: Convergence
    1: Failed to converge, too many iterations
    2: Overflow error in calculating the mole fractions
  */

  int i, j, i_tmp, j_tmp, k; // indices, i over particles and j over compounds
  int i_point; // the current point in the titration
  int cur_rank;
  int *active_particles = NULL; // Particles with nonzero concentration
  int *active_compounds = NULL; // Compounds that have nonzero concentration
  int *num_active = NULL;
  double *grad = NULL;
  double *abs_tol = NULL;
  double *x_tmp = NULL; // Temporary x for storing nonzero concentrations
  double *new_x0 = NULL; // x0 containing only particle counts
  double *x0_tmp = NULL; // Temporary x0 for storing nonzero entries in x0
  double *G_tmp = NULL; // Temporary G for storing free energies of active cmpds
  double *A_tmp = NULL; // Temporary x for storing info about active compounds
  int n_compounds_tmp; // Temporary number of compounds
  int n_particles_tmp; // Temporary number of particles
  int ret_val = ERR_OK;
  int changed = 1;
  run_stats_t * act_run_stats = NULL;
 
  /* ********** CUT OUT PARTICLES WITH ZERO CONC  ************************** */
  // Allocate memory for temporary variables
  active_particles = (int *) malloc(n_particles * sizeof(int));
  active_compounds = (int *) malloc(n_compounds * sizeof(int));
  num_active = (int *) malloc(n_particles * sizeof(int));

  // Build new problem, start with memory allocation
  A_tmp = (double *) malloc(n_particles * n_compounds * sizeof(double));
  x_tmp = (double *) malloc(n_compounds * sizeof(double));
  G_tmp = (double *) malloc(n_compounds * sizeof(double));
  x0_tmp = (double *) malloc(n_compounds * sizeof(double));
  new_x0 = (double *) malloc(n_particles * sizeof(double));
  grad = (double *) malloc(n_particles * sizeof(double));
  abs_tol = (double *) malloc(n_particles * sizeof(double));

  if(NULL == run_stats) {
    act_run_stats = (run_stats_t *) malloc(sizeof(run_stats_t));
  } 
  else {
    act_run_stats = run_stats;
  }


  if(NULL == active_particles ||
     NULL == active_compounds ||
     NULL == num_active ||
     NULL == A_tmp ||
     NULL == x_tmp || 
     NULL == G_tmp ||
     NULL == x0_tmp || 
     NULL == act_run_stats ||
     NULL == grad ||
     NULL == abs_tol) {
    ret_val = ERR_OOM;
  }

  if(ret_val == ERR_OK) {

    for (i_point = 0; i_point < n_points; i_point++) {
      // Find total particle mass for pruning
      for (i = 0; i < n_particles; i++) {
        new_x0[i] = 0.0;
        for (j = 0; j < n_compounds; j++) {
          new_x0[i] += x0[ij(i_point, j, n_compounds)] 
            * A[ij(i, j, n_compounds)];
        }
      }
      // Determine which particles have nonzero concentration
      for (i = 0; i < n_particles; i++) {
        if (new_x0[i] < DBL_MIN && -new_x0[i] < DBL_MIN) {
          active_particles[i] = 0;
        }
        else {
          active_particles[i] = 1;
        }
      }

      changed = 1;

      for (j = 0; j < n_compounds; j++) {
        active_compounds[j] = 1;
      }


      while (changed) {
        changed = 0;
        // Disable compounds if they contain an inactive particle
        for (i = 0; i < n_particles; i++) {
          if (active_particles[i] == 0) {
            for (j = 0; j < n_compounds; j++) {
              if (active_compounds[j] && 
                  !(A[ij(i, j, n_compounds)] < DBL_MIN && 
                    -A[ij(i, j, n_compounds)] < DBL_MIN)) {
                active_compounds[j] = 0;
                changed = 1;
              }
            }
          }
        }
        
        // Only enable particles if there is more than one active compound
        // that contains them
        for (i = 0; i < n_particles; i++) {
          num_active[i] = 0;
        }
        for (i = 0; i < n_particles; i++) {
          for (j = 0; j < n_compounds; j++) {
            if (active_compounds[j] &&
                (A[ij(i, j, n_compounds)] > DBL_MIN || 
                -A[ij(i, j, n_compounds)] < DBL_MIN)) {

              num_active[i] ++;
            }
          }
        }
        for (i = 0; i < n_particles; i++) {
          if (num_active[i] <= 1 && active_particles[i]) {
            changed = 1;
            active_particles[i] = 0;
          }
        }
      }


      for (k = 0; k < n_particles; k++) {
        if (active_particles[k]) {
          n_particles_tmp = sumint(active_particles, k+1);
          n_compounds_tmp = sumint(active_compounds, n_compounds);
          i_tmp = 0;
          for (i = 0; i < k+1; i++) {
            if (active_particles[i]) {
              j_tmp = 0;
              for (j = 0; j < n_compounds; j++) {
                if (active_compounds[j]) {
                  A_tmp[ij(i_tmp, j_tmp, n_compounds_tmp)]
                    = A[ij(i, j, n_compounds)];
                  j_tmp++;
                }
              }
              i_tmp++;
            }
          }
          cur_rank = get_approx_rank(A_tmp, n_particles_tmp, 
              n_compounds_tmp, NUM_PRECISION);
          if (cur_rank != n_particles_tmp) {
            active_particles[k] = 0;
          }
        }
      }
      
      n_particles_tmp = sumint(active_particles, n_particles);
      n_compounds_tmp = sumint(active_compounds, n_compounds);

      if (n_particles_tmp != 0 && n_compounds_tmp != 0) {
        // Build new A and x0, A_tmp and x0_tmp, respectively
        i_tmp = 0;
        for (i = 0; i < n_particles; i++) {
          if (active_particles[i]) {
            j_tmp = 0;
            for (j = 0; j < n_compounds; j++) {
              if (active_compounds[j]) {
                A_tmp[ij(i_tmp, j_tmp++, n_compounds_tmp)] 
                   = A[ij(i, j, n_compounds)];
              }
            }
            i_tmp++;
          }
        }
        
        // Build new G, x0
        j_tmp = 0;
        for (j = 0; j < n_compounds; j++) {
          if (active_compounds[j]) {
            x0_tmp[j_tmp] = x0[ij(i_point, j, n_compounds)];
            G_tmp[j_tmp++] = G[j];
          }
        }
      
        /* *************** END CUTTING OUT ZERO CONCS *********************** */
        ret_val = calc_conc_optimize(x_tmp, A_tmp, G_tmp, x0_tmp, 
              n_particles_tmp, n_compounds_tmp, max_iters, tol, delta_bar, eta, 
              min_delta, max_trial, perturb_scale, seed, act_run_stats, grad, 
              abs_tol);
      }

      /* *************** CONVERT BACK TO ORIGINAL PROBLEM ******************* */
      // Give the full gradient, zero where x0 was zero
      i_tmp = n_particles_tmp - 1;
      for (i = n_particles - 1; i >= 0; i--) {
        if (active_particles[i]) {
          abs_tol[i] = abs_tol[i_tmp];
          grad[i] = grad[i_tmp--];
        }
        else {
          abs_tol[i] = 0;
          grad[i] = 0.0;
        }
      }

      // Report concentrations
      j_tmp = 0;
      for (j = 0; j < n_compounds; j++) {
        if (active_compounds[j]) {
          x[ij(i_point, j, n_compounds)] = x_tmp[j_tmp++];
        }
        else {  // Not active. Use original concentration
          x[ij(i_point, j, n_compounds)] = x0[ij(i_point, j, n_compounds)];
        }
      }

      /* *************** DONE CONVERTING BACK TO ORIGINAL PROBLEM *********** */

      // Report errors in conservation of mass to screen
      if (quiet == 0) {
        // Print out values of the gradient, which is the error in cons. of mass
        printf("Error in conservation of mass (units of mole fraction):\n");
        for (i = 0; i < n_particles; i++) {
          printf("   %8.6e \n", grad[i]);
        }
        printf("\n");
      }

      if (write_log_file == 1) {
        write_concentration_log_file(act_run_stats, grad, abs_tol, log_file);
      }
    }
  }
  if(active_particles != NULL) free(active_particles);
  if(active_compounds != NULL) free(active_compounds);
  if(num_active != NULL) free(num_active);

  if(A_tmp != NULL) free(A_tmp);
  if(x_tmp != NULL) free(x_tmp); 
  if(G_tmp  != NULL) free(G_tmp); 
  if(x0_tmp  != NULL) free(x0_tmp); 
  if(new_x0 != NULL) free(new_x0);
  if(grad != NULL) free(grad);
  if(abs_tol != NULL) free(abs_tol);

  if(run_stats == NULL && act_run_stats != NULL) {
        free(act_run_stats);
  }

  return ret_val;
}
/* ************************************************************************** */
/*                 END CALC_CONC_FROM_FREE_ENERGIES                           */
/* ************************************************************************** */

/* ************************************************************************** */
/*                        BEGIN CALC_CONC_OPTIMIZE                            */
/* ************************************************************************** */
int calc_conc_optimize(double *x, double *A, double *G, double *x0,
          int n_constraints, int n_compounds, int max_iters, double tol,
          double delta_bar, double eta, double min_delta, int max_trial,
          double perturb_scale, unsigned long seed, run_stats_t * run_stats,
          double *final_gradient, double *absolute_tolerance) {

  int i, j ; // indices, i over particles and j over all compounds
  int iters = 0; // Number of iterations
  double *abs_tol; // The absolute tolerance on all values of gradient
  double rho; // Ratio of actual to predicted reduction in trust region method
  double delta; // Radius of trust region
  double *lambda; // Lagrange multipliers
  double *p; // The step we take toward minimization
  double *hes; // The Hessian
  unsigned long rand_seed = 0; // Random number seed
  double *new_x0; // Constraint values
  double *new_lambda; // Adjusted lambda if taking Newton step
  double new_log_scale_fact = 0.0; // Adjusted scale factor if taking Newt step
  double *new_x; // Adjusted x if taking Newton step
  double *new_grad; // Adjusted gradient if taking Newton step
  double *grad;
  double log_scale_fact = 0.0; // Log of scale factor dividing the x values

  double n_old; // Old norm of grad
  double n_new; // New norm of grad
  int n_trial; // Number of times we've perturbed lambda
  int newton_success; // Whether Newton's method was successful
  int ret_val = ERR_OK;
  
  // Allocate memory
  new_x0 = (double *) malloc(n_constraints * sizeof(double));
  hes = (double *) malloc(n_constraints * n_constraints * sizeof(double));
  abs_tol = (double *) malloc(n_constraints * sizeof(double));
  grad = (double *) malloc(n_constraints * sizeof(double));
  lambda = (double *) malloc(n_constraints * sizeof(double));
  p = (double *) malloc(n_constraints * sizeof(double));

  new_grad = (double *) malloc(n_constraints * sizeof(double));
  new_lambda = (double *) malloc(n_constraints * sizeof(double));
  new_x = (double *) malloc(n_compounds * sizeof(double));

  run_stats->max_n_trials = max_trial;
  run_stats->n_constraints = n_constraints;

  if(NULL == new_x0 || 
     NULL == hes || 
     NULL == abs_tol || 
     NULL == grad || 
     NULL == lambda || 
     NULL == p ||
     NULL == new_grad ||
     NULL == new_lambda || 
     NULL == new_x) {
    ret_val = ERR_OOM;
goto end_calc_conc_optimize;
  } 

  // Make an "x0" that is the constraints, i.e., dot(A, x0).
  for(i = 0; i < n_constraints ; i++) {
    new_x0[i] = 0.0;
    for(j = 0; j < n_compounds; j++) {
      new_x0[i] += A[ij(i,j,n_compounds)] * x0[j];
    }
  }

  n_trial = 0;

  // The absolute tolerance is a percentage of the entries in x0
  make_tol(abs_tol, A, x0, n_constraints, n_compounds, tol);

  for (i = 0; i < n_constraints; i++) {
    grad[i] = abs_tol[i] + 1.0; // Initialize just to get started.
  }

  // Run trust region to get solution
  while (check_tol(grad, abs_tol, n_constraints) == 0 && n_trial < max_trial) {

    if (n_trial == 1) {
      // Seed the random number generator if necessary
      rand_seed = get_random_seed(seed);
      init_genrand(rand_seed);
    }

    // Set initial guess
    if (get_initial_guess(new_x0, lambda, G, A, n_constraints, 
                          n_compounds, perturb_scale, rand_seed)) {
      ret_val = ERR_INITIAL;
goto end_calc_conc_optimize;
    }

    // Calculate the counts of the species based on lambda
    if (get_x(x, &log_scale_fact, lambda, G, A, n_constraints, n_compounds)){ 
      // Should be fine; checked prev.
      // If it overflows on the initial guess, we probably won't be able to
      // find a valid initial guess
      ret_val = ERR_OVERFLOW;
#ifdef DEBUG
      printf("Overflow error at location 1, trial = %d.\n", n_trial);
      printf("lambda = ");
      for (i = 0; i < n_constraints; i++) {
	printf("%.6e  ", lambda[i]);
      }
      printf("\nG = ");
      for (j = 0; j < n_compounds; j++) {
	printf("%.6e  ", G[j]);
      }
      printf("\nA = \n");
      for (i = 0; i < n_constraints; i++) {
	for (j = 0; j < n_compounds; j++) {
	  printf("%.8f  ", A[ij(i, j, n_compounds)]);
	}
	printf("\n");
      }
      printf("\n");
#endif

goto end_calc_conc_optimize;
    }

    // Calculate the gradient
    get_grad(grad, log_scale_fact, new_x0, x, A, n_constraints, n_compounds);
    
    // Compute the Hessian (symmetric, positive, positive definite)
    ret_val = get_hessian(hes, x, A, n_constraints, n_compounds);
    if(ret_val) goto end_calc_conc_optimize;
      
    // Initialize delta to be just less than delta_bar
    delta = 0.99 * delta_bar;
    
    // Initializations
    run_stats->n_newton_steps = 0; // Number of pure Newton steps
    run_stats->n_cauchy_steps = 0; // Number of pure Cauchy steps
                                   // (hit trust region boundary)
    run_stats->n_dogleg_steps = 0; // Number of dogleg steps 
                                   // (part Newton and part Cauchy)
    run_stats->n_chol_fail_cauchy_steps = 0; // Cholesky failure forcing 
                                             // Cauchy step
    run_stats->n_irrel_chol_fail = 0; // Number of steps with irrelovent 
                                      // Cholesky failures
    run_stats->n_dogleg_fail = 0; // Number of failed dogleg calculations

    make_tol(abs_tol, A, x, n_constraints, n_compounds, tol);


#ifdef DEBUG
    printf("log_scale_fact = %.6e\n", log_scale_fact);
    printf("lambda = ");
    for (i = 0; i < n_constraints; i++) {
      printf("%.6e  ", lambda[i]);
    }
    printf("\nx = ");
    for (j = 0; j < n_compounds; j++) {
      printf("%.6e  ", x[j]);
    }
    printf("\ngradient = ");
    for (i = 0; i < n_constraints; i++) {
      printf("%.6e  ", grad[i]);
    }
    printf("\n");
#endif


    // Run trust region with these initial conditions
    while (iters < max_iters && check_tol(grad, abs_tol, n_constraints) == 0 
           && delta > min_delta) {

      // Solve for the search direction
      ret_val = get_search_dir(p, grad, hes, delta, n_constraints, run_stats);
      if(ret_val) goto end_calc_conc_optimize;
      
      // Calculate rho, ratio of actual to predicted reduction
      rho = get_rho(lambda, p, grad, x, hes, new_x0, G, A, log_scale_fact,
                    n_constraints, n_compounds);
      
      // Adjust delta and make step based on rho
      if (rho < 0.25) {
        delta /= 4.0;
      }
      else if (rho > 0.75 
               && fabs(norm(p, n_constraints) - delta) < NUM_PRECISION) {
        delta = min2(2.0 * delta, delta_bar);
      }
      if (rho > eta) {
        for (i = 0; i < n_constraints; i++) {
          lambda[i] += p[i];
        }

        // Calculate the mole fractions of the complexes based on lambda
        if (get_x(x, &log_scale_fact, lambda, G, A, 
		  n_constraints, n_compounds)) {
          // Should be fine; checked prev.
          ret_val = ERR_OVERFLOW;
#ifdef DEBUG
	  printf("**********Overflow error at location 2.***********\n");
#endif
goto end_calc_conc_optimize;
        }
      
        // Calculate the gradient
        get_grad(grad, log_scale_fact, new_x0, x, A, n_constraints, 
		 n_compounds);

        // Compute the Hessian (symmetric, positive, positive definite)
        ret_val = get_hessian(hes, x, A, n_constraints, n_compounds);
        if(ret_val) goto end_calc_conc_optimize;
      }
      
      make_tol(abs_tol, A, x, n_constraints, n_compounds, tol);
      // Advance the iterations count
      iters++;
    }

    // If we shrank the trust region too much, try last Newton step procedure
    if (delta <= min_delta) {
      newton_success = 1;

      while (newton_success && check_tol(grad, abs_tol, n_constraints) == 0 
             && iters < max_iters) {
        // Attempt Newton step
        if (compute_newton_step(p, grad, hes, n_constraints) == 0) {
          // Compute new lambda to see if the gradient decreased
          for (i = 0; i < n_constraints; i++) {
            new_lambda[i] = lambda[i] + p[i];
          }

          // Calculate the mole fractions of the complexes based on new lambda
          ret_val = get_x(new_x, &new_log_scale_fact, new_lambda, G, A, 
			  n_constraints, n_compounds);
#ifdef DEBUG
	  if (ret_val == ERR_OVERFLOW) {
	    printf("**********Overflow error at location 3.**********\n");
	  }
#endif
          if(ret_val != ERR_OK) goto end_calc_conc_optimize;
          
          // Calculate the new gradient
          get_grad(new_grad, new_log_scale_fact, new_x0, new_x, A, 
		   n_constraints, n_compounds);
          
          // If we got better with the Newton step, accept it.
          n_new = norm(new_grad, n_constraints);
          n_new *= exp(-new_log_scale_fact);
          n_old = norm(grad, n_constraints);
          n_old *= exp(-log_scale_fact);

          if (n_new < n_old) {
            for (j = 0; j < n_compounds; j++) {
              x[j] = new_x[j];
            }
            for (i = 0; i < n_constraints; i++) {
              lambda[i] = new_lambda[i];
              grad[i] = new_grad[i];
	      log_scale_fact = new_log_scale_fact;
            }

            // Calculate the new Hessian
            ret_val = get_hessian(hes, x, A, n_constraints, n_compounds);
            if(ret_val) goto end_calc_conc_optimize;

            run_stats->n_newton_steps++;
          }
          else { // Gradient increased with Newton step; reject it.
            newton_success = 0;
          }
        }
        else { // Failed to invert matrix in Newton step.
          newton_success = 0;
        }
        make_tol(abs_tol, A, x, n_constraints, n_compounds, tol);
        iters++;
      }
    }
    // Advance the number of perturbations we've tried
    n_trial++;
  }
  run_stats->n_iterations = iters;
  run_stats->n_trials = n_trial;

  make_tol(abs_tol, A, x, n_constraints, n_compounds, tol);
  if (check_tol(grad, abs_tol, n_constraints) == 0) { // failed, too many iters
    ret_val = ERR_NOCONVERGE;
  }

  // Prepare final results
  for(i = 0; i < n_constraints; i++) {
    final_gradient[i] = grad[i] * exp(log_scale_fact);
    absolute_tolerance[i] = abs_tol[i] * exp(log_scale_fact);
  }
  for (j = 0; j < n_compounds; j++) {
    x[j] *= exp(log_scale_fact);
  }


#ifdef DEBUG
    printf("final log_scale_fact = %.6e\n", log_scale_fact);
    printf("final lambda = ");
    for (i = 0; i < n_constraints; i++) {
      printf("%.6e  ", lambda[i]);
    }
    printf("\nfinal gradient = ");
    for (i = 0; i < n_constraints; i++) {
      printf("%.6e  ", grad[i]);
    }
    printf("\nfinal tolerance = ");
    for (i = 0; i < n_constraints; i++) {
      printf("%.6e  ", abs_tol[i]);
    }
    printf("\nfinal x = ");
    for (j = 0; j < n_compounds; j++) {
      printf("%.6e  ", x[j]);
    }
    printf("\n");
    printf("n_iters = %d\n", iters);
    printf("\nA = \n");
    for (i = 0; i < n_constraints; i++) {
      for (j = 0; j < n_compounds; j++) {
	printf("%.6e   ", A[ij(i, j, n_compounds)]);
      }
      printf("\n");
    }
    printf("\nG = ");
    for (j = 0; j < n_compounds; j++) {
      printf("%.6e  ", G[j]);
    }
    printf("\nx0 = ");
    for (j = 0; j < n_compounds; j++) {
      printf("%.6e  ", x0[j]);
    }
    printf("\n\n");
#endif



end_calc_conc_optimize:
  // Free memory
  if(new_x0 != NULL) free(new_x0);
  if(hes != NULL) free(hes);
  if(abs_tol != NULL) free(abs_tol);
  if(grad != NULL) free(grad);
  if(lambda != NULL) free(lambda);
  if(p != NULL) free(p);
    
  if(new_grad != NULL) free(new_grad); 
  if(new_lambda != NULL) free(new_lambda ); 
  if(new_x != NULL) free(new_x ); 

  // Return convergence
  return ret_val;
}
/* ************************************************************************** */
/*                          END CALC_CONC_OPTIMIZE                            */
/* ************************************************************************** */


/* ************************************************************************** */
int get_G_from_K(double *G, double *minus_log_K, double *augmented_N, 
		 int n_reactions, int n_compounds) {
  /*
    Computes the free energy of the compounds given the stoichiometric
    matrix, N, and the minus log of the equilibrium constants for the chemical
    reactions, minus_log_K.

  */

  int i, j; // indices
  double *augmented_minus_log_K = NULL; // Augmented array of equilibrium constants
  int * augmented_N_perm = NULL; // Augmented N permutation matrix (for LUP decomposition)
  double *augmented_N_copy = NULL;
  int ret_val = ERR_OK;
  int n_constraints = n_compounds-n_reactions;

  // Memory allocation and error checking
  augmented_minus_log_K = (double *) 
                      malloc(n_compounds * sizeof(double));
  augmented_N_perm = (int *) 
                      malloc(n_compounds * sizeof(int));
  augmented_N_copy = (double *) 
                      malloc(n_compounds * n_compounds * sizeof(double));

  if(NULL == augmented_minus_log_K ||
     NULL == augmented_N_perm ||
     NULL == augmented_N_copy) {
    ret_val = ERR_OOM;
  }


  if(ret_val == ERR_OK) {
    // Build first n_constraints rows of augmented_N and augmented_minus_log_K
    for (i = 0; i < n_constraints; i++) {
      augmented_minus_log_K[i] = 0.0;
    }

    // Build rest of augmented_N and augmented_minus_log_K
    for (i = n_constraints; i < n_compounds; i++) {
      augmented_minus_log_K[i] = minus_log_K[i-n_constraints];
    }

    for (i = 0; i < n_compounds; i++) {
      for (j = 0; j < n_compounds; j++) {
        augmented_N_copy[ij(i, j, n_compounds)] = 
                  augmented_N[ij(i, j, n_compounds)];
      }
    }

    // Solve for the free energies by LUP decomposition
    lup_decomposition(augmented_N_copy, augmented_N_perm, n_compounds);

    ret_val = lup_solve(augmented_N_copy, augmented_N_perm, n_compounds, 
			augmented_minus_log_K, G);
  }

  // Free memory. Avoid freeing NULL pointer
  if(NULL != augmented_minus_log_K) free(augmented_minus_log_K);
  if(NULL != augmented_N_perm) free(augmented_N_perm);
  if(NULL != augmented_N_copy) free(augmented_N_copy);

  return ret_val;
}
/* ************************************************************************** */

/* ************************************************************************** */
int get_A(double *A, double *N, int n_reactions, int n_compounds) {
  int i,j;
  int rank = 0;
  int ret_val = ERR_OK;
  int n_constraints = n_compounds-n_reactions;
  ret_val = get_nullspace(A, &rank, N, n_reactions, n_compounds);
  if(rank != n_constraints) {
    ret_val = ERR_BADN;
  }
  if(ret_val == ERR_OK) {
    for(i = 0; i < n_reactions; i++) {
      for(j = 0; j < n_compounds; j++) {
        A[ij(i+n_constraints, j, n_compounds)] = N[ij(i, j, n_compounds)];
      }
    }
  }

  return ret_val;
}
/* ************************************************************************** */

/* ************************************************************************** */
int get_nullspace(double *N, int *N_rank, double *A, int m, int n) {
  /*
   * Calculates the null space of m x n matrix A, returns a basis of
   * the null space as row vectors
   *
   * This uses a householder reflections to perform QR decomposition on
   * A^T and reconstructs the null space by backwards accumulation of
   * Q, only accumulating the rows of Q which specify the null space.
   *
   * Note that the use of i and j can be confusing because I am 
   * decomposing A^T instead of A and i and j are referring to the
   * row and column of A^T
   *
   * arguments
   * =========
   * N : a preallocated matrix at least n x n in order to store all the
   *      rows necessary for a null matrix. Only N_rank*n rows will be
   *      written. A smaller matrix can be used if the rank of A is already
   *      known
   * 
   * N_rank : on return, the location pointed to by *N_rank will be set
   *      to the rank of N (which will be equal to the number of 
   *      rows in N)
   *
   * A : the matrix to get the null space of
   *
   * m : the number of rows in A
   * n : the number of columns in A
   */
  int i,j,k;
  double * v = NULL; // the householder vector
  double * w = NULL; // the transformed householder vector
  double * betas = NULL; // list of betas
  double beta, sigma, mu;// Householder constants
  double gamma;     // Temporary scaling variable
  int rank = m;        // rank of null space
  int ret_val = ERR_OK;
  double * A_copy = NULL;

  if(NULL == (A_copy = (double*)malloc(n*m*sizeof(double)))) {
    ret_val = ERR_OOM;
  }

  if(NULL == (v = (double*)malloc(n*sizeof(double)))) {
    ret_val = ERR_OOM;
  }
  if(NULL == (w = (double*)malloc(n*sizeof(double)))) {
    ret_val = ERR_OOM;
  }
  if(NULL == (betas = (double*)malloc(m*sizeof(double)))) {
    ret_val = ERR_OOM;
  }

  if(ret_val == ERR_OK) {
    for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
        A_copy[ij(i,j,n)] = A[ij(i,j,n)];
      }
    }
    // Compute QR decomposition
    for(j = 0; j < m; j++) {
      // Compute Householder vector
      sigma = 0;
      for(i = j; i < n; i++) {
        v[i] = A_copy[ij(j,i,n)];
      }
      for(i = j+1; i < n; i++) {
        sigma += v[i]*v[i];
      }
      v[j] = 1.0;

      if(sigma < NUM_PRECISION) {
        beta = 0;
      } else {
        gamma = A_copy[ij(j,j,n)];
        mu = sqrt(gamma*gamma + sigma);
        if(gamma <= 0) {
          gamma = gamma - mu;
        } else {
          gamma = -sigma / (gamma + mu);
        }
        beta = (2.0*gamma*gamma) ;
        beta /= (sigma+(gamma*gamma));
        for(i = j+1; i < n; i++) {
          v[i] = v[i] / gamma;
        }
      }
      betas[j] = beta;

      // Transform Householder vector
      for(i = j; i < m; i++) {
        w[i] = 0;
        for(k = j; k < n; k++) {
          w[i] += v[k]*A_copy[ij(i,k,n)];
        }
        w[i] *= beta;
      }
      // Perform Householder update
      for(i = j; i < m; i++) {
        for(k = j; k < n; k++) {
          A_copy[ij(i,k,n)] -= v[k]*w[i];
        }
      }
      // Save the essential part of the vector
      for(i = j+1; i < n; i++) {
        A_copy[ij(j,i,n)] = v[i];
      }
    }

    rank = n; 
    for(i = 0; i < m; i++) {
      if(A_copy[ij(i,i,n)] >= NUM_PRECISION) {
        rank --;
      }
    }

    // Select the last m rows of Q
    for(i = 0; i < rank; i++) {
      for(j = 0; j < n; j++) {
        N[ij(i,j,n)] = 0;
      }
      N[ij(i,n-rank+i,n)] = 1.0;
    }

    // Left multiply by A and back-accumulate
    for(j = m-1; j >= 0; j--) {
      v[j] = 1.0;
      for(k = j+1; k < n; k++) {
        v[k] = A_copy[ij(j,k,n)];
      }

      beta = betas[j];

      // Compute w
      for(i = 0; i < rank; i++) {
        w[i] = 0;
        for(k = j; k < n; k++) {
          w[i] += N[ij(i,k,n)]*v[k];
        }
        w[i] *= beta;
      }

      // Accumulate in N
      for(i = 0; i < rank; i++) {
        for(k = j; k < n; k++) {
          N[ij(i,k,n)] -= w[i]*v[k];
        }
      }
    }
  }
  *N_rank = rank;

  if(NULL != A_copy) free(A_copy);
  if(NULL != v) free(v);
  if(NULL != w) free(w);
  if(NULL != betas) free(betas);

  return ret_val;
}
/* ************************************************************************** */

/* ************************************************************************** */
int prune_N_K_x(double *N_new, double *minus_log_K_new, double * x0_new, 
            int *n_reactions_new, int *n_compounds_new, int *active_compounds, 
            double *N_old, double *minus_log_K_old, double *x0_old, 
            int n_reactions, int n_compounds) {
  /*
   * Prunes N and K based on Justin's integer mass-action integration
   * algorithm.
   *
   * Input Arguments
   * ---------------
   * N_old: pointer to the N reaction stoichiometry matrix to prune
   *        size = (n_reactions, n_compounds)
   * minus_log_K_old: pointer to the -log(K) reaction constant to 
   *        prune size=n_reactions
   * x0_old: pointer to the initial compound mole fractions
   *
   * Output Arguments
   * ----------------
   *  N_new: @start pointer to empty array
   *        size = (n_reactions, n_compounds)
   *         @end pointer to the pruned N reaction stoichiometry matrix
   *  minus_log_K_new: @start pointer to empty array
   *        size = (n_reactions)
   *         @end pointer to the pruned K reaction constant array
   *  x0_new: @start pointer to a free array of doubles
   *           size = (n_compounds)
   *         @end: pointer to the pruned initial compound mole fractions
   *  n_reactions_new: @start pointer to a free integer
   *         @end pointer to the number of reactions after pruning
   *  n_compounds_new: @start pointer to free integer
   *         @end pointer to the number of compounds after pruning
   *  active_compounds: @start pointer to empty array
   *        size = (n_compounds)
   *          @end pointer to an array of 0s and 1s. 
   *              If active_compounds[i] == 1, then compound i is included
   *                in the active compounds. Otherwise, it has been pruned
   *
   * Return
   * ------
   *  ERR_OK if everything worked out fine. An error code > 0 otherwise.
   */

  int done;
  int *active_reactions = NULL;

  int i,j,k,l;
  int f_rate;
  int b_rate;
  int ret_val = ERR_OK;
  int n_reacs_new;
  int n_comps_new;

  if(NULL == (active_reactions = (int*)malloc(n_reactions*sizeof(int)))) {
    ret_val = ERR_OOM;
  }

  if(ret_val == ERR_OK) {
    // Justin's algorithm
    for(j = 0; j < n_compounds; j++) {
      if(x0_old[j] > 0) {
        active_compounds[j] = 1;
      } else {
        active_compounds[j] = 0;
      }
    }
    for(i = 0; i < n_reactions; i++) {
      active_reactions[i] = 0;
    }
    
    done = 0;
    while(!done) {
      done = 1;
      for(i = 0; i < n_reactions; i++) {
        f_rate = 1;
        b_rate = 1;
        for(j = 0; j < n_compounds; j++) {
          if(N_old[ij(i,j,n_compounds)] < 0) {
            f_rate *= active_compounds[j];
          } else if(N_old[ij(i,j,n_compounds)] > 0) {
            b_rate *= active_compounds[j];
          }
        }
        if(f_rate > 0 && !active_reactions[i]) {
          active_reactions[i] = 1;
          done = 0;
          for(j = 0; j < n_compounds; j++) {
            if(N_old[ij(i,j,n_compounds)] > 0) {
              active_compounds[j] = 1;
            }
          }
        } else if(b_rate > 0 && !active_reactions[i]) {
          done = 0;
          for(j = 0; j < n_compounds; j++) {
            if(N_old[ij(i,j,n_compounds)] < 0) {
              active_compounds[j] = 1;
            }
          }
        }
      }
    }
    for(j = 0; j < n_compounds; j++) {
      active_compounds[j] = 0;
    }

    for(i = 0 ; i < n_reactions; i++) {
      if(active_reactions[i]) {
        for(j = 0; j < n_compounds; j++) {
          if(N_old[ij(i,j,n_compounds)] > 0 || 
              N_old[ij(i,j,n_compounds)] < 0 ) {
            active_compounds[j] = 1;
          }
        }
      }
    }
    n_reacs_new = 0;
    n_comps_new = 0;
    for(i = 0; i < n_reactions; i++) {
      n_reacs_new += active_reactions[i];
    }
    for(j = 0; j < n_compounds; j++) {
      n_comps_new += active_compounds[j];
    }
    *n_reactions_new = n_reacs_new;
    *n_compounds_new = n_comps_new;
    

    // Copy N, minus_log_K, x0 over
    // Could do all the copying without branching, I think simplicity
    // of memory allocation wins over simplicity of algorithm for this 
    // case.
    k = 0;
    for(i = 0; i < n_reactions; i++) {
      l = 0;
      if(active_reactions[i]) {
        for(j = 0; j < n_compounds; j++) {
          if(active_compounds[j]) {
            N_new[ij(k,l,n_comps_new)] = N_old[ij(i,j,n_compounds)];
            l++;
          }
        }
        minus_log_K_new[k] = minus_log_K_old[i];
        k++;
      }
    }
    l = 0;
    for(j = 0; j < n_compounds; j++) {
      if(active_compounds[j]) {
        x0_new[l] = x0_old[j];
        l++;
      }
    }
  }

  if(NULL != active_reactions) free(active_reactions);

  return ret_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
int calc_conc(double *x, double *N, double *minus_log_K,
              double *x0, int n_reactions, int n_compounds, 
              int n_points,
              int max_iters, double tol, double delta_bar, 
              double eta, double min_delta,
              int max_trial, double perturb_scale, 
              int quiet, int write_log_file, char *log_file, 
              unsigned long seed, run_stats_t *run_stats) {
  /*
    Calculates concentrations given N, log(K), and x0.
  */

  double *G = NULL; // The free energies
  double *A = NULL;
  int return_val = ERR_OK; // Keep track of error state
  int n_reactions_new, n_compounds_new; // pruned dimensions
  int i_point; // The index of the current point (in initial conditions)
  int *active_compounds = NULL; // [i] == 1 if compound i is active 0 o.w.
  double *N_new = NULL;   // Pruned stoichiometry matrix
  double *minus_log_K_new = NULL;   // Pruned equilibrium constants
  double *x0_tmp = NULL;  // Current initial concentration point
  double *x0_new = NULL;  // Pruned initial concentrations
  double *x_new = NULL;   // Pruned final concentrations
  double *final_gradient = NULL;
  double *absolute_tolerance = NULL;
  int *p = NULL; // p (pivots) for use in LUP decomposition
  double *log_x = NULL; // log of concentrations
  int j, l;

  run_stats_t * act_run_stats = NULL;;

  G = (double *) malloc(n_compounds*sizeof(double));
  A = (double *) malloc(n_compounds*n_compounds*sizeof(double));
  N_new = (double *) malloc(n_compounds*n_reactions*sizeof(double));
  minus_log_K_new = (double *) malloc(n_reactions*sizeof(double));
  active_compounds = (int *) malloc(n_compounds*sizeof(int));
  x0_tmp = (double *) malloc(n_compounds*sizeof(double)); 
  x0_new = (double *) malloc(n_compounds*sizeof(double));
  x_new = (double *) malloc(n_compounds*sizeof(double));
  final_gradient = (double *) malloc((n_compounds-n_reactions)*sizeof(double));
  absolute_tolerance = (double *) malloc((n_compounds-n_reactions)*sizeof(double));

  if(NULL == run_stats) {
    act_run_stats = (run_stats_t *) malloc(sizeof(run_stats_t));
    if(NULL != act_run_stats) {
    }
  } else {
    act_run_stats = run_stats;
  }

  if( NULL == G ||
      NULL == A ||
      NULL == N_new ||
      NULL == minus_log_K_new ||
      NULL == active_compounds ||
      NULL == x0_tmp ||
      NULL == x0_new ||
      NULL == x_new ||
      NULL == act_run_stats ||
      NULL == final_gradient ||
      NULL == absolute_tolerance) {
    return_val = ERR_OOM;
  }

  for (i_point = 0; i_point < n_points; i_point++) {
    for (j = 0; j < n_compounds; j++) {
      x0_tmp[j] = x0[ij(i_point, j, n_compounds)];
    }

    if (return_val == ERR_OK) {
      return_val = prune_N_K_x(N_new, minus_log_K_new, x0_new, 
			       &n_reactions_new, 
			       &n_compounds_new, active_compounds, 
			       N, minus_log_K, x0_tmp, n_reactions, 
			       n_compounds);
    }

    if (return_val == ERR_OK) {
      return_val = get_A(A, N_new, n_reactions_new, n_compounds_new);
    }
    if (return_val == ERR_OK && n_reactions_new > 0 && n_compounds_new > 0) {
      if (n_reactions_new == n_compounds_new) { // Fully constrained, just solve
	p = (int *) malloc(n_compounds_new*sizeof(int));
	log_x = (double *) malloc(n_compounds_new*sizeof(double));
	if (p == NULL || log_x == NULL) {
	  return_val = ERR_OOM;
	}
	else {
	  lup_decomposition(A, p, n_compounds_new);
	  return_val = lup_solve(A, p, n_compounds_new, minus_log_K_new, log_x);
	  if (return_val == 0){
	    for (j = 0; j < n_compounds_new; j++) {
	      if (log_x[j] > MAXLOGX) { // Will have an overflow error
		return_val = ERR_OVERFLOW;
#ifdef DEBUG
		printf("**********Overflow error at location 4**********.\n");
#endif
	      }
	      else {
		x_new[j] = exp(log_x[j]);
	      }
	    }
	  }
	}
      }
      else {
	return_val = get_G_from_K(G, minus_log_K_new, A, n_reactions_new, 
				  n_compounds_new);
      }
    } 

    if (return_val == ERR_OK && n_reactions_new > 0 && n_compounds_new > 0
	&& n_reactions_new != n_compounds_new) {

      // Solve for concentrations
      return_val = calc_conc_optimize(x_new, A, G, x0_new, 
				      n_compounds_new-n_reactions_new, 
				      n_compounds_new,  
				      max_iters, tol, 
				      delta_bar, eta, min_delta,
				      max_trial, perturb_scale, 
				      seed, act_run_stats,
				      final_gradient,
				      absolute_tolerance);
    }

    if(return_val == ERR_OK) {
      l = 0;
      for(j = 0; j < n_compounds; j++) {
        if(active_compounds[j]) {
          x[ij(i_point, j, n_compounds)] = x_new[l];
          l++;
        } else {
          x[ij(i_point, j, n_compounds)] = x0_tmp[j];
        }
      }
    }
  }

  if(quiet == 0 && return_val != ERR_OK) {
    write_error_message(return_val);
  }

  if(write_log_file && return_val != ERR_OOM) {
    write_concentration_log_file(act_run_stats, final_gradient, 
                                 absolute_tolerance, log_file);
  }

  // Free memory
  if(G != NULL) free(G);
  if(A != NULL) free(A);
  if(N_new != NULL) free(N_new);
  if(minus_log_K_new != NULL) free(minus_log_K_new);
  if(x0_new != NULL) free(x0_new);
  if(x0_tmp != NULL) free(x0_tmp);
  if(x_new != NULL) free(x_new);
  if(active_compounds != NULL) free(active_compounds);

  if(absolute_tolerance != NULL) free(absolute_tolerance);
  if(final_gradient != NULL) free(final_gradient);
  if(p != NULL) free(p);
  if(log_x != NULL) free(log_x);

  // If run_stats were NULL, free the temporary run stats
  if(run_stats == NULL) {
    if(act_run_stats != NULL) free(act_run_stats);
  }
  
  return return_val;
}
/* ************************************************************************** */

/* ************************************************************************** */
int calc_titration(double *x, double *N, double *minus_log_K,
                   double *x0, double *x0_titrated, double *vol_titrated,
                   int volume_titration, double initial_volume,  
                   double x_titrant, 
                   int n_titration_points, int titrated_species,
                   int n_reactions, int n_compounds, 
                   int max_iters, double tol, double delta_bar, 
                   double eta, double min_delta,
                   int max_trial, double perturb_scale, 
                   int quiet) {
  /*
    Calculates a titration curve.  Stores results in array x.

    If volume_titration is True, assumes we are doing a titration
    where volume is constantly added to the solution.  The units of
    initial_volume and vol_titrated are irrelovant; they just have to
    be the same.

    THIS CODE CAN BE MODIFIED SO THAT WE DON'T DO QR DECOMPOSITION
    ON EVERY ITERATION.
  */

  int i,j;
  int write_log_file = 0; // Do not write log file
  unsigned long seed = 0; // Seed off clock
  char * dummy_log_file = NULL; // Dummy log file
  run_stats_t * no_run_stats = NULL;

  double *n0 = NULL; // Total number of particles of species
  double *x0_new = NULL; // Input x0 for calc_conc
  double *x_new = NULL;  // Result of calc_conc
  int return_val = ERR_OK;

  if(NULL == (x0_new = (double *) malloc(n_compounds*sizeof(double)))) {
    return_val = ERR_OOM;
  } 
  else if(NULL == (x_new = (double *) malloc(n_compounds*sizeof(double)))) {
    return_val = ERR_OOM;
  }

  if(return_val == ERR_OK) {
    if (volume_titration) {
      if (NULL == (n0 = (double *) malloc(n_compounds*sizeof(double)))) {
        return_val = ERR_OOM;
      }

      // Compute the number of particles of each species
      for (j = 0; j < n_compounds; j++) {
        n0[j] = x0[j] * initial_volume;
      }

      // Loop through, adjusting for added volume each time
      for (i = 0; i < n_titration_points; i++) {
        for (j = 0; j < n_compounds; j++) {
          x0_new[j] = n0[j] / (initial_volume + vol_titrated[i]);
        }
        x0_new[titrated_species] = vol_titrated[i] * x_titrant 
          / (initial_volume + vol_titrated[i]);

        return_val = calc_conc(x_new, N, minus_log_K, x0_new, n_reactions, 
			       n_compounds, 
                               1, max_iters, tol, delta_bar, eta, min_delta, 
                               max_trial, perturb_scale, quiet, write_log_file, 
                               dummy_log_file, seed, no_run_stats);

        for (j = 0; j < n_compounds; j++) {
          x[ij(i,j,n_compounds)] = x_new[j];
        }
        
        if (return_val != ERR_OK) {
          break;
        }
      }
      
    }
    else {
      for (j = 0; j < n_compounds; j++) {
        x0_new[j] = x0[j];
      }
      
      // Loop through, putting titrated concs in x0
      for (i = 0; i < n_titration_points; i++) {
        x0_new[titrated_species] = x0_titrated[i];
        return_val = calc_conc(x_new, N, minus_log_K, x0_new, n_reactions, 
			       n_compounds, 
                               1, max_iters, tol, delta_bar, eta, min_delta, 
                               max_trial, perturb_scale, quiet, write_log_file, 
                               dummy_log_file, seed, no_run_stats);
        if (return_val != ERR_OK) {
          break;
        }
        for (j = 0; j < n_compounds; j++) {
          x[ij(i,j,n_compounds)] = x_new[j];
        }        
      }
    }
  }


  // Free memory
  if(n0 != NULL) free(n0);
  if(x0_new != NULL) free(x0_new);
  if(x_new != NULL) free(x_new);

  return return_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
int calc_titration_from_free_energies(
                       double *x, double *A, double *G,
                       double *x0, double *x0_titrated, double *vol_titrated,
                       int volume_titration, double initial_volume,  
                       double x_titrant, 
                       int n_titration_points, int titrated_species,
                       int n_constraints, int n_compounds,
                       int max_iters, double tol, double delta_bar, 
                       double eta, double min_delta,
                       int max_trial, double perturb_scale, 
                       int quiet) {
  /*
    Calculates a titration curve.  Stores results in array x.

    If volume_titration is True, assumes we are doing a titration
    where volume is constantly added to the solution.  The units of
    initial_volume and vol_titrated are irrelovant; they just have to
    be the same.
  */

  int i,j;
  int write_log_file = 0; // Do not write log file
  unsigned long seed = 0; // Seed off clock
  char * dummy_log_file = NULL; // Dummy log file
  run_stats_t * no_run_stats = NULL;

  double *n0 = NULL; // Total number of particles of species
  double *x0_new = NULL; // Input x0 for calc_conc_from_free_energies
  double *x_new = NULL;  // Result of calc_conc_from_free_energies
  int return_val = ERR_OK;

  if (NULL == (x0_new = (double *) malloc(n_compounds*sizeof(double)))) {
    return_val = ERR_OOM;
  } 
  if (NULL == (x_new = (double *) malloc(n_compounds*sizeof(double)))) {
    return_val = ERR_OOM;
  }

  if (return_val == ERR_OK) {
    if (volume_titration) {
      if (NULL == (n0 = (double *) malloc(n_compounds*sizeof(double)))) {
        return_val = ERR_OOM;
      }

      // Compute the number of particles of each species
      for (j = 0; j < n_compounds; j++) {
        n0[j] = x0[j] * initial_volume;
      }

      // Loop through, adjusting for added volume each time
      for (i = 0; i < n_titration_points; i++) {
        for (j = 0; j < n_compounds; j++) {
          x0_new[j] = n0[j] / (initial_volume + vol_titrated[i]);
        }
        x0_new[titrated_species] = vol_titrated[i] * x_titrant
          / (initial_volume + vol_titrated[i]);

        return_val = calc_conc_from_free_energies(x_new, A, G, x0_new, 
                                                  n_constraints, 
                                                  n_compounds, 1, 
                                                  max_iters, tol, 
                                                  delta_bar, eta, min_delta, 
                                                  max_trial, perturb_scale, 
                                                  quiet, write_log_file, 
                                                  dummy_log_file, 
                                                  seed, no_run_stats);
        if(return_val != ERR_OK) {
          break;
        }
        for(j = 0; j < n_compounds; j++) {
          x[ij(i,j,n_compounds)] = x_new[j];
        }
      }
    }
    else{

      // Copy inputting x0 to one that changes for titrated species
      for (j = 0; j < n_compounds; j++) {
        x0_new[j] = x0[j];
      }
      
      // Loop through, putting titrated concs in x0
      for (i = 0; i < n_titration_points; i++) {
        x0_new[titrated_species] = x0_titrated[i];
        return_val = calc_conc_from_free_energies(x_new, A, G, x0_new, 
                                                  n_constraints, 
                                                  n_compounds, 1,
                                                  max_iters, tol, 
                                                  delta_bar, eta, min_delta, 
                                                  max_trial, perturb_scale, 
                                                  quiet, write_log_file, 
                                                  dummy_log_file, 
                                                  seed, no_run_stats);
        if(return_val != ERR_OK) {
          break;
        }
        for(j = 0; j < n_compounds; j++) {
          x[ij(i,j,n_compounds)] = x_new[j];
        }
      }
    }
  }

  // Free memory
  if(x0_new != NULL) free(x0_new);
  if(x_new != NULL) free(x_new);

  return return_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
int get_initial_guess(double *x0, double *lambda, double *G, double *A, 
                      int n_constraints, int n_compounds, 
                      double perturb_scale, unsigned long rand_seed) {
  /*
    Pick initial lambda such that ln x = 1 for all x (x ~ 3).

    This is done by solving:
      A . A_transpose . lambda = A . (G + 1)

    lambda must be preallocated to have n_constraints entries.
  */

  int i, j, k; // indices
  double log_x = 1.0; // log of the mole fractions
  double *c = NULL; // Right hand side of eqtn to solve.
  int *p = NULL; // Permutation array for LUP decomposition
  double *G_plus_log_x = NULL;
  double *A_AT = NULL; // A times its transpose
  double sum_parts = 0;
  int n_nonzero = 0;
  int ret_val;  // Return value of LUP solve


  // Allocate memory
  G_plus_log_x = (double *) malloc(n_compounds * sizeof(double));
  c = (double *) malloc(n_constraints * sizeof(double));
  p = (int *) malloc(n_constraints * sizeof(int));
  A_AT = (double *) malloc(n_constraints * n_constraints * sizeof(double));
  if (NULL == A_AT ||
      NULL == G_plus_log_x ||
      NULL == c ||
      NULL == p) {
    return ERR_OOM;
  }

  // Compute G + log_x, a vector
  for (j = 0; j < n_compounds; j++) {
    G_plus_log_x[j] = G[j] + log_x;
  }

  // c is the right hand side of equation to solve.
  matrix_vector_mult(c, A, G_plus_log_x, n_constraints, n_compounds);

  // Compute matrix multiplying lambda (usually identity, but not assuming that)
  for (i = 0; i < n_constraints; i++) {
    for (j = 0; j < n_constraints; j++) {
      A_AT[ij(i, j, n_constraints)] = 0.0;
      for (k = 0; k < n_compounds; k++) {
        A_AT[ij(i, j, n_constraints)] 
          += A[ij(i, k, n_compounds)] * A[ij(j, k, n_compounds)];
      }
    }
  }

  // A_AT must be positive semi-definite, 
  // could use modified cholesky
  // Solve using LUP decomposition
  modified_cholesky(A_AT, p, n_constraints);  
  ret_val = modified_cholesky_solve(A_AT, p, n_constraints, c, lambda);
  
  // Perturb lambda if desired
  if (rand_seed != 0) {
    perturb_lambda(lambda, perturb_scale, G, A, n_constraints, n_compounds);
  }

  // If we already know concentration (particle is inert), set lambda
  for (i = 0; i < n_constraints; i++) {
    sum_parts = 0;
    n_nonzero = 0;
    for (j = 0; j < n_compounds; j++) {
      sum_parts += A[ij(i, j, n_compounds)];
      n_nonzero ++;
    }

    if (n_nonzero == 1) {
      // Find the first nonzero entry in the row of A
      j = 0;
      while (j < n_compounds && A[ij(i, j, n_compounds)] < DBL_MIN 
          && -A[ij(i, j, n_compounds)] < DBL_MIN) {
        j++;
      }
      
      lambda[i] = (log(x0[i]) + G[j]) / sum_parts;
    }
  }

  /*
  // THIS STUFF IS BRIAN'S HACK TO ADJUST LAMBDA IF WE HAVE OVERFLOW
  // THIS SHOULD NOT BE NECESSARY WITH OUR NEW SCALED CONCENTRATIONS METHOD
  double dot_prod;
  double magnitude;
  double diff;
  for (j = 0; j < n_compounds; j++) {
    dot_prod = 0.0;
    magnitude = 0.0;
    for (i = 0; i < n_constraints; i++) {
      dot_prod += lambda[i] * A[ij(i, j, n_compounds)];
      magnitude += fabs(A[ij(i, j, n_compounds)]);
    }
    log_x = -G[j] + dot_prod;
    if (log_x > MAXLOGX / 2) {
      diff = log_x - MAXLOGX / 2;
      diff /= magnitude;
      for (i = 0; i < n_constraints; i++) {
        lambda[i] -= diff * A[ij(i, j, n_compounds)];
      }
    }
  }
  */

  // Free memory. Avoid freeing NULL pointer
  if(NULL != G_plus_log_x) free(G_plus_log_x);
  if(NULL != A_AT) free(A_AT);
  if(NULL != p) free(p);
  if(NULL != c) free(c);

  return ret_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
int get_x(double *x, double *log_scale_fact, double *lambda, double *G, 
	  double *A, int n_constraints, int n_compounds) {
  /* 
     Calculates the mole fractions of all species from lambda, G, and
     A.  Returns 0 if the calculation was ok and ERR_OVERFLOW if there will
     be an overflow error.

    Gives a scaled x and log_scale_fact.  If x_real is the "real"
    concentrations, then x = x_real / exp(log_scale_fact).
  */

  int i, j; // indices
  int ret_val = 0; // Value to be returned
  double *log_x = NULL; // log of the mole fraction
  double log_scale_factor; // The log of the scaling factor for the concs
  double dot_prod; // Dot product
  double exp_arg; // Argument for exponential, log_x - log_scale_fact

  // Allocate memory
  if (NULL == (log_x = (double *) malloc(n_compounds * sizeof(double)))) {
    return ERR_OOM;
  }

  // Compute log_x for first compound
  dot_prod = 0.0;
  for (i = 0; i < n_constraints; i++) {
    dot_prod += lambda[i] * A[ij(i, 0, n_compounds)];
  }
  log_x[0] = -G[0] + dot_prod;
  log_scale_factor = log_x[0];

  // Compute log_x for the rest of the compounds
  for (j = 1; j < n_compounds; j++) {
    dot_prod = 0.0;
    for (i = 0; i < n_constraints; i++) {
      dot_prod += lambda[i] * A[ij(i, j, n_compounds)];
    }
    log_x[j] = -G[j] + dot_prod;
    if (log_x[j] > log_scale_factor) {
      log_scale_factor = log_x[j];
    }
  }

  // Compute scaled x's
  for (j = 0; j < n_compounds; j++) {
    exp_arg = log_x[j] - log_scale_factor;
    if (exp_arg > MAXLOGX) {
      ret_val = ERR_OVERFLOW;
goto end_get_x;
    }
    x[j] = exp(exp_arg);
  }

  *log_scale_fact = log_scale_factor;

end_get_x:
  // Free memory
  if(NULL != log_x) free(log_x);

  // No overflow errors
  return ret_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
void get_grad(double *grad, double log_scale_fact, double *x0, double *x, 
	      double *A, int n_constraints, int n_compounds) {
  /*
    Calculates the gradient of -h(\lambda), the dual function for
    which we're trying to find the minimum.

    The inputted x is the scaled x.  If x_real is the "real"
    concentrations, then x = x_real / exp(log_scale_fact).

    The outputted gradient is also scaled.
  */

  int i, j; // indices
  double dot_prod;

  for (i = 0; i < n_constraints; i++) {
    dot_prod = 0.0;
    for (j = 0; j < n_compounds; j++) {
      dot_prod += x[j] * A[ij(i, j, n_compounds)];
    }
    grad[i] = -x0[i] * exp(-log_scale_fact) + dot_prod;
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
int get_hessian(double *hes, double *x, double *A, int n_constraints, 
                 int n_compounds) {
  /*
    Calculates the Hessian, hes.  Must be preallocated and of size
    n_constraints by n_constraints.  We only need to calculate the
    upper triangle of hes since it's symmetric.  We can then fill out
    the lower triangle.

    The inputted x is the scaled x.  If x_real is the "real"
    concentrations, then x = x_real / exp(log_scale_fact).

    Avec is a vector consisting of elements of one row of A multiplied by 
    those of another.
  */


  int m, n, j; // index
  double *Avec;
  int ret_val = ERR_OK;

  Avec = (double *) malloc(n_compounds * sizeof(double));

  if(NULL == Avec) {
    ret_val = ERR_OOM;
  } else {
    for (n = 0; n < n_constraints; n++) {
      for (m = 0; m <= n; m++) {
        for (j = 0; j < n_compounds; j++) {
          Avec[j] = A[ij(m, j, n_compounds)] * A[ij(n, j, n_compounds)];
        }
        hes[ij(m, n, n_constraints)] = dot(x, Avec, n_compounds);
      }
    }

    // Fill out the lower entries in the Hessian (needed for matrix mults)
    for(m = 1; m < n_constraints; m++) {
      for (n = 0; n < m; n++) {
        hes[ij(m, n, n_constraints)] = hes[ij(n, m, n_constraints)];
      }
    }
  }

  if(Avec != NULL) free(Avec);
  return ret_val;
}
/* ************************************************************************** */


/* ************************************************************************** */
int get_search_dir(double *p, double *grad, double *hes, double delta, 
                   int n_constraints, run_stats_t * run_stats) {
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
      5 if Cholesky decompostion failed but we would've taken Cauchy step anyway
      6 if the dogleg calculation failed (should never happen)
  */

  int i; // index
  double a, b, c, q; // Constants used in quadratic formula
  double t1, t2, t3; // Temporary variables used in quadratic formula
  double tau; // Multipler in Newtonstep
  double beta; // = tau - 1
  double *pB; // Unconstrained minimizer (the regular Newton step)
  double *pU; // Minimizer along steepest descent direction
  double pB2; // pj^2
  double pU2; // pu^2
  double pBpU; // pj . pc
  double *hes_dot_grad; // Hessian dotted with the gradient
  double pUcoeff;  // The coefficient on the Cauchy step
  double delta2; // delta^2
  int newton_fail; // = 1 if Newton step failed
  int ret_val = ERR_OK;

  // Initialize pB2 so compiler doesn't give warning when optimization is on
  pB2 = 0.0;

  // Useful to have delta^2 around.
  delta2 = pow(delta, 2);

  // Allocate all memory arrays
  pB = (double *) malloc(n_constraints * sizeof(double));
  hes_dot_grad = (double *) malloc(n_constraints * sizeof(double));
  pU = (double *) malloc(n_constraints * sizeof(double));
  if(pB == NULL || hes_dot_grad == NULL || pU == NULL) {
    ret_val = ERR_OOM;
goto end_get_search_dir;
  }

  /* ********** Compute the Newton step ************** */

  newton_fail = compute_newton_step(pB, grad, hes, n_constraints);
  if (newton_fail == 0) {
    // If Newton step is within trust region, take it
    pB2 = dot(pB, pB, n_constraints);
    if (pB2 <= delta2) {
      for (i = 0; i < n_constraints; i++) {
        p[i] = pB[i];
      }
      run_stats->n_newton_steps ++;
goto end_get_search_dir;
    }
  }
  /* ************************************************* */

  /* ********** Compute the Cauchy step ************** */

  // The direction of the Cauchy step
  for (i = 0; i < n_constraints; i++) {
    pU[i] =  -grad[i];
  }

  // prefactor for the Cauchy step
  matrix_vector_mult(hes_dot_grad, hes, grad, n_constraints, n_constraints);
  pUcoeff = dot(grad, grad, n_constraints) ;
  pUcoeff /= dot(grad, hes_dot_grad, n_constraints);
  for (i = 0; i < n_constraints; i++) {
    pU[i] = pUcoeff * pU[i];
  }

  pU2 = dot(pU, pU, n_constraints);

  if (pU2 >= delta2) { // In this case we take the Cauchy step, 0 < tau <= 1
    tau = sqrt(delta2 / pU2);
    for (i = 0; i < n_constraints; i++) {
      p[i] = tau * pU[i];
    }
    if (newton_fail != 0) {
      // Cholesky failure, doesn't matter, would take Cauchy anyway 
      run_stats->n_irrel_chol_fail ++;
    } else {
      // Signifies that we just took the Cauchy step
      run_stats->n_cauchy_steps ++;
    }
goto end_get_search_dir; 
  }
  
  if (newton_fail != 0) { // Failed computing Newton step & must take Cauchy
    for (i = 0; i < n_constraints; i++) {
      p[i] = pU[i];
    }
    run_stats->n_chol_fail_cauchy_steps ++;
goto end_get_search_dir;
  }
  /* ************************************************* */

  /* ************ Take the dogleg step *************** */
  pBpU = dot(pB, pU, n_constraints); // Need this for dogleg calculation
  
  // Constants for quadratic formula to solve ||pU + beta (pB-pU)||^2 = delta2,
  // where beta = tau - 1.
  a = pB2 + pU2 - 2.0 * pBpU;  // a > 0
  b = 2.0 * (pBpU - pU2);      // b can be positive or negative
  c = pU2 - delta2;            // c <= 0 since pU2 <= delta2 to get here
  t1 = pow(b, 2);
  t2 = sqrt(t1 - 4.0 * a * c);
  t3 = sgn(b) * t2;
  q = -0.5 * (b + t3);
  // q = -0.5 * (b + sgn(b) * sqrt(pow(b, 2) - 4.0 * a * c));


  // Choose correct (positive) root (don't have to worry about a = 0 because
  // if pU \approx pB, we would have already taken Newton step
  if (fabs(b) < NUM_PRECISION) {
    beta = sqrt(-c / a);
  }
  else if (b < 0.0) {
    beta = q / a;
  }
  else {   // b > 0
    beta = c / q;
  }

  if (beta >= 0.0 && beta <= 1.0) { // This is ok and should always be the case
    for (i = 0; i < n_constraints; i++) {
      p[i] = pU[i] + beta * (pB[i] - pU[i]);
    }
    run_stats->n_dogleg_steps ++;
  }
  else { // Something is messed up, take Cauchy step (we should never get here)
    for (i = 0; i < n_constraints; i++) {
      p[i] = pU[i];
    }
    run_stats->n_dogleg_fail ++;
  }

end_get_search_dir:

  free(pB);
  free(pU);
  free(hes_dot_grad);
  return ret_val;
  /* ************************************************* */
}
/* ************************************************************************** */


/* ************************************************************************** */
int compute_newton_step(double *pB, double *grad, double *hes, 
                        int n_constraints) {
  /*
    Computes a Newton step, given the Hessian and the gradient
    (assuming the Hessian is positive definite).  Returns 0 is the
    Hessian is not positive definite (Cholesky decomposition fails)
    and 1 if the Newton step is successfully computed.  Stores the
    Newton step in pB.
  */

  int i, j; // indices
  double *hes_copy; // Copy of the upper triangle of the Hessian 
                    // (don't want to mess with actual Hessian).
  int *p;
  int chol_success; // Whether of not Cholesky decomposition is successful

  // Make a copy of the Hessian because the Cholesky decomposition messes
  // with its entries.  We only have to copy the upper diagonal.
  hes_copy = (double *) malloc(n_constraints * n_constraints * sizeof(double));
  p = (int *) malloc(n_constraints * sizeof(int));
 
  for (j = 0; j < n_constraints; j++) {
    for (i = j; i < n_constraints; i++) {
      hes_copy[ij(i, j, n_constraints)] = hes[ij(i, j, n_constraints)];
    }
  }

  chol_success = modified_cholesky(hes_copy, p, n_constraints);

  if (chol_success == 0) {
    modified_cholesky_solve(hes_copy, p, n_constraints, grad, pB);

    // Free memory from the Hessian
    free(hes_copy);
    free(p);

    // Newton step is -H^{-1} grad
    for (i = 0; i < n_constraints; i++) {
      pB[i] *= -1.0;
    }

    return 0;
  }
  else { // Free memory from failed Newton step computation and return 0
    for(i = 0 ; i < n_constraints; i++) {
      for(j = 0; j < n_constraints; j++) {
        hes_copy[ij(i, j, n_constraints)] = hes[ij(i, j, n_constraints)];
      }
    }
    lup_decomposition(hes_copy, p, n_constraints);
    if(0 == lup_solve(hes_copy, p, n_constraints, grad, pB)) {
      for(i = 0; i < n_constraints; i++) {
        pB[i] *= -1.0;
      }
      free(hes_copy);
      free(p);
      return 0;
    }
    free(hes_copy);
    free(p);
    return 1;
  }
}
/* ************************************************************************** */


/* ************************************************************************** */
double get_rho(double *lambda, double *p, double *grad, double *x, double *hes, 
               double *x0, double *G, double *A, double log_scale_fact, 
	       int n_constraints, int n_compounds) {
  /*
    Calculates rho based on equations 4.4 and 4.1 of Nocedal and
    Wright.  This is the ratio of the actual correction based on the
    stepping a length delta along the search direction to the
    predicted correction of taking the same step.

    Function returns -1 if there is an overflow error in the
    calculation of rho.
  */

  int i; // index
  double rho; // That which we return
  double *new_lambda; // y after the step
  double *new_x; // x after the Newton step
  double *hes_dot_p; // the vector hes*p
  double obj_fun; // -h(lambda)
  double new_obj_fun; // New one
  double p_dot_hes_dot_p; // p * H * p
  double new_log_scale_fact; // New scale factor with new x
  double l_dot_x0;


  // Allocate memory for arrays used in the calculation.
  new_lambda = (double *) malloc(n_constraints * sizeof(double));
  new_x = (double *) malloc(n_compounds * sizeof(double));
  hes_dot_p = (double *) malloc(n_constraints * sizeof(double));

  l_dot_x0 = dot(lambda, x0, n_constraints);
  l_dot_x0 *= exp(-log_scale_fact);

  obj_fun = sum(x, n_compounds);
  obj_fun -= l_dot_x0;

  // Calculate the new lambda
  for (i = 0; i < n_constraints; i++) {
    new_lambda[i] = lambda[i] + p[i];
  }

  // Calculate the counts of the species based on new lambda
  if (get_x(new_x, &new_log_scale_fact, new_lambda, G, A, n_constraints, 
	    n_compounds) == 0) {

    new_obj_fun = sum(new_x, n_compounds) ;
    new_obj_fun *= exp(new_log_scale_fact - log_scale_fact);

    l_dot_x0 = dot(new_lambda, x0, n_constraints) ;
    l_dot_x0 *= exp(-log_scale_fact);
    new_obj_fun -= l_dot_x0;

    matrix_vector_mult(hes_dot_p, hes, p, n_constraints, n_constraints);
    p_dot_hes_dot_p = dot(p, hes_dot_p, n_constraints);

    rho = (obj_fun - new_obj_fun) 
      / fabs(-dot(grad, p, n_constraints) - p_dot_hes_dot_p / 2.0);
  }
  else {
    // Since the denominator is always positive, an overflow error in the
    // calculation of the mole fractions for the "new lambda" implies
    // a negative denominator, and therefore the value of rho is negative.

#ifdef DEBUG
    printf("Overflow error at location 5.\n");
#endif

    rho = -1.0;
  }

  free(new_lambda);
  free(new_x);
  free(hes_dot_p);

  return rho;
}
/* ************************************************************************** */


/* ************************************************************************** */
void get_cauchy_point(double *cauchy_point, double *hes, double *grad, 
                      double delta, int n_constraints) {
  /*
   Computes the Cauchy point using the formulas 4.7 and 4.8 in Nocedal
   and Wright, Numerical Optimization (1999), page 70.

   Note that because the Hessian is symmetric that it is equal to its
   transpose.
  */

  int i; // index
  double coeff; // Coefficient on gradient for Cauchy point
  double tau; // Multiplier used in Cauchy point calculation
  double norm_grad; // ||Grad||
  double *hes_dot_grad; // Hessian dotted with the gradient
  
  hes_dot_grad = (double *) malloc (n_constraints * sizeof(double));

  norm_grad = norm(grad, n_constraints);
  matrix_vector_mult(hes_dot_grad, hes, grad, n_constraints, n_constraints);
  tau = pow(norm_grad, 3);
  tau /= delta * dot(grad, hes_dot_grad, n_constraints);
  tau = min2(tau, 1.0);
  coeff = -tau * delta / norm_grad;

  for (i = 0; i < n_constraints; i++) {
    cauchy_point[i] = coeff * grad[i];
  }
   
  free(hes_dot_grad);
}
/* ************************************************************************** */


/* ************************************************************************** */
int perturb_lambda(double *lambda, double perturb_scale, double *G, double *A, 
                int n_constraints, int n_compounds) {
  /*
    Perturbs the values of lambda in case the trust region has shrunk to be 
    very small.

    Adds perturb_scale*random number to each entry in lambda
  */

  int i; // index
  double *dummy_x; // A dummy mole fraction vector for checking overflow
  double dummy_log_scale_fact; // Dummy variable for checking overflow
  double *new_lambda; // The new perturbed lambda
  int xOK; // = 0 is there is no overflow error induced in the concentrations
  int ret_val = ERR_OK;
  int counter = 0;
  int max_counter = 10000;

  dummy_x = (double *) malloc(n_compounds * sizeof(double));
  new_lambda = (double *) malloc(n_constraints * sizeof(double));

  if(NULL == dummy_x ||
     NULL == new_lambda) {
    ret_val = ERR_OOM;
  }

  xOK = 1;
  while (xOK == 1 && ret_val == ERR_OK && counter < max_counter) {
    for (i = 0; i < n_constraints; i++) {
      new_lambda[i] = lambda[i] + perturb_scale * 2.0 * (genrand_real1() - 0.5);
    }
    xOK = get_x(dummy_x, &dummy_log_scale_fact, new_lambda, G, A, 
		n_constraints, n_compounds);
    perturb_scale /= 2.0; // Reduce scale to try not to have overflow problems
  }

  if(counter == max_counter) {
    ret_val = ERR_OVERFLOW;
#ifdef DEBUG
    printf("**********Overflow error at location 6.**********\n");
#endif
  }

  // Copy new perturbed lambda to the value of lambda
  for (i = 0; i < n_constraints; i++) {
    lambda[i] = new_lambda[i];
  }
  
  free(dummy_x);
  free(new_lambda);
  return ret_val;
}
/* ************************************************************************** */

/* ************************************************************************** */
void make_tol(double *abs_tol, double *A, double *x, int n_constraints, 
        int n_compounds, double tol) {
  /*
    Determines what the absolute tolerance is for each of the constraints.
  */

  int i, j;

  for (i = 0; i < n_constraints; i++) {
    abs_tol[i] = 0.0;
    for (j = 0; j < n_compounds; j++) {
      abs_tol[i] += tol * fabs(A[ij(i, j, n_compounds)]) * x[j];
    }
  }
}
/* ************************************************************************** */

/* ************************************************************************** */
int check_tol(double *grad, double *abs_tol, int n_constraints) {
  /*
    Check the to see if the entries in the gradient satisfy the
    tolerance.  Returns 1 if they all do and 0 otherwise.
  */

  int i; // index

  for (i = 0; i < n_constraints; i++) {
    if (fabs(grad[i]) > abs_tol[i]) {
      return 0;
    }
  }


  return 1;
}
/* ************************************************************************** */

/* ************************************************************************** */
void write_error_message(int error_code) {
  /*
   * Write the error message corresponding to error_code to stderr
   */
  switch(error_code) {
    case ERR_NOCONVERGE:
      fprintf(stderr,"\nOptimization failed to converge\n");
      break;
    case ERR_OVERFLOW:
      fprintf(stderr,"\nNumerical overflow encountered\n");
      break;
    case ERR_BADN:
      fprintf(stderr,"\nReactions are linearly dependent\n");
      break;
    case ERR_BADK:
      fprintf(stderr,"\nInvalid values for equilibrium constants K provided\n");
      break;
    case ERR_NOINPUT:
      fprintf(stderr,"\nNo valid input files specified\n");
      break;
    case ERR_NOCMP:
      fprintf(stderr,"\nNo composition .cmp file found\n");
    case ERR_SOLVDENS:
      fprintf(stderr,"\nError in solvent density\n");
      break;
    case ERR_BADFREEENERGY:
      fprintf(stderr,"\nError in free energy input\n");
      break;
    case ERR_CON:
      fprintf(stderr,"\nError in concentration input specification\n");
    case ERR_LOG:
      fprintf(stderr,"\nError writing to the log file\n");
      break;
    case ERR_EQ:
      fprintf(stderr,"\nError opening output file .eq\n");
      break;
    case ERR_INITIAL:
      fprintf(stderr, "\nError setting initial guess.\n");
      break;
    case ERR_OOM:
      fprintf(stderr,"\nError allocating memory\n");
    default:
      break;
  }
}
/* ************************************************************************** */

/* ************************************************************************** */
int write_concentration_log_file(run_stats_t *run_stats, double * grad, 
                              double * abs_tol, char *log_file) {
  /*
    Writes the calculation statistics from an equilibrium
    concentrations calculation to the log file.
  */

  int i; // Counter
  FILE *fplog;  // file handle for log file

  if ((fplog = fopen(log_file, "a")) == NULL) {
    fprintf(stderr,"Error opening %s.\n\nExiting....\n", log_file);
    // Return failure to open log file
    return ERR_LOG;
  }

  if (check_tol(grad,abs_tol,run_stats->n_constraints)) {
    fprintf(fplog, "TRUST REGION DID NOT CONVERGE DUE TO PRECISION ISSUES\n\n");
  }

  fprintf(fplog, "   --Trust region results:\n");
  fprintf(fplog, "       No. of initial conditions tried: %d\n", run_stats->n_trials);
  fprintf(fplog, "       Results from final trial:\n");
  fprintf(fplog, "         No. of iterations: %d\n", run_stats->n_iterations);
  fprintf(fplog, "         No. of Newton steps: %d\n", 
                                                    run_stats->n_newton_steps);
  fprintf(fplog, "         No. of Cauchy steps: %d\n", 
                                                    run_stats->n_cauchy_steps);
  fprintf(fplog, "         No. of dogleg steps: %d\n", 
                                                    run_stats->n_dogleg_steps);
  fprintf(fplog, 
          "         No. of Cholesky failures resulting in Cauchy steps: %d\n",
          run_stats->n_chol_fail_cauchy_steps);
  fprintf(fplog, "         No. of inconsequential Cholesky failures: %d\n",
          run_stats->n_irrel_chol_fail);
  fprintf(fplog, "         No. of dogleg failures: %d\n", 
                                                    run_stats->n_dogleg_fail);

  fprintf(fplog,
          "         Error in conservation of mass (units of molarity):\n");
  fprintf(fplog, "              Error\tTolerance (units of mole fraction)\n");
  for (i = 0; i < run_stats->n_constraints; i++) {
    fprintf(fplog, "           %.14e\t%.14e\n", 
            grad[i], abs_tol[i]);
  }
  fclose(fplog);

  return 0;
}
/* ************************************************************************** */
