#pragma once

/*
  concentrations.h is part of the NUPACK software suite
  Copyright (c) 2007 Caltech. All rights reserved.

  Header file with global variables and function prototypes for use
  with concentrations.c Contains trust region parameters, among
  others.
*/

#include "utils.h"
  
#include <getopt.h>// Takes options from the command line
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif
// Version of EQCON
#define EQCONVERSION "0.05"

// Constants used in trust region
#define TRUST_REGION_DELTABAR 1000.0 // Maximal size of trust region
#define TRUST_REGION_ETA 0.125 // Decision criterion for trust region

#define MAXLOGX 250 // Maximum logarithm of a concentration (prevents overflow)

/* *************************** Struct Definitions *************************** */
typedef struct _RUN_STATS_T_ {
  // Run statistics
  int max_n_trials;
  int n_trials;
  int n_constraints;
  
  // Trial statistics
  int n_iterations;
  int n_newton_steps;
  int n_cauchy_steps;
  int n_dogleg_steps;
  int n_chol_fail_cauchy_steps;
  int n_irrel_chol_fail;
  int n_dogleg_fail;
} run_stats_t;
/* ************************************************************************** */


// Function prototypes
/* ***************** IN EQUILIBRIUM_CONCENTRATIONS.C  *********************** */
int calc_conc_from_free_energies(double *x, double *A, double *G, double *x0,
                                 int n_particles, int n_compounds, 
                                 int n_points, int max_iters, double tol, 
                                 double delta_bar, double eta, 
                                 double min_delta, int max_trial, 
                                 double perturb_scale, int quiet, 
                                 int write_log_file, char *log_file, 
				 unsigned long seed, run_stats_t * run_stats);
int calc_conc(double *x, double *N, double *minus_log_K,
	      double *x0, int n_reactions, int n_compounds, 
              int n_points, int max_iters, double tol, 
              double delta_bar, double eta, double min_delta,
	      int max_trial, double perturb_scale, 
	      int quiet, int write_log_file, char *log_file, 
	      unsigned long seed, run_stats_t * run_stats);
int calc_conc_optimize(double *x, double *A, double *G, double *x0,
              int n_particles, int n_compounds, int max_iters,
              double tol, double delta_bar, double eta,
              double min_delta, int max_trial, double perturb_scale,
              unsigned long seed, run_stats_t * run_stats,
              double * gradient, double * abs_tol);
int get_G_from_K(double *G, double *minus_log_K, double *N, int n_compounds, 
		 int n_particles);
int prune_N_K_x(double *N_new, double *minus_log_K_new, double * x0_new, 
            int *n_reactions_new, int *n_compounds_new, int *active_compounds, 
            double *N_old, double *minus_log_K_old, double *x0_old, 
            int n_reactions, int n_compounds);
int get_N_prime(double *N_prime, double *N, 
                int n_reactions, int n_compounds);
int get_nullspace(double *N, int *N_rank, double *A, int m, int n);
int get_A(double *A, double *N, int n_reactions, int n_compounds);
int calc_titration(double *x, double *N, double *minus_log_K,
                   double *x0, double *x0_titrated, double *vol_titrated,
		   int volume_titration, double initial_volume,  
		   double x_titrant, 
                   int n_titration_points, int titrated_species,
                   int n_reactions, int n_compounds, 
                   int max_iters, double tol, double delta_bar, 
                   double eta, double min_delta,
                   int max_trial, double perturb_scale, 
                   int quiet);
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
                       int quiet);
int get_initial_guess(double *x0, double *lambda, double *G, double *A, 
		      int n_particles, int n_compounds, double perturb_scale, 
		      unsigned long rand_seed);
int get_x(double *x, double *log_scale_fact, double *lambda, double *G, 
	  double *A, int n_particles, int n_compounds);
void get_grad(double *grad, double log_scale_fact, double *x0, double *x, 
	      double *A, int n_particles, int n_compounds);
int get_hessian(double *hes, double *x, double *A, int n_particles, 
		 int n_compounds);
int get_search_dir(double *p, double *grad, double *hes, double delta, 
		   int n_particles, run_stats_t * run_stats);
int compute_newton_step(double *pB, double *grad, double *hes, 
			int n_particles);
double get_rho(double *lambda, double *p, double *grad, double *x, double *hes, 
	       double *x0, double *G, double *A, double log_scale_fact, 
	       int n_particles, int n_compounds);
void get_cauchy_point(double *cauchy_point, double *hes, double *grad, 
		      double delta, int n_particles);
int perturb_lambda(double *lambda, double perturb_scale, double *G, double *A, 
		int n_particles, int n_compounds);
int check_tol(double *grad, double *abs_tol, int n_particles);
void make_tol(double *abs_tol, double *A, double *x, int n_constraints, 
        int n_compounds, double tol);
int write_concentration_log_file(run_stats_t * run_stats, double * grad,
                            double * abs_tol, char *log_file);
void write_error_message(int key);
/* ************************************************************************** */

/* ************************ IN READ_COMMAND_LINE.C ************************** */
void read_command_line(int nargs, char **args, char *spec_file, char *con_file, 
		       char *log_file, char *eq_file, int *sort_output,
		       int *max_iters, double *tol, int *free_energy_input,
		       double *solvent_number_density, double *kT, 
		       int *dimensionless, double *min_delta, int *max_trial, 
		       double *perturb_scale, int *quiet, int *write_log_file, 
		       unsigned long *seed);
void display_help(void);
/* ************************************************************************** */

/* ************************* IN READ_INPUT_FILES.C ************************** */
void get_size(int *n_particles, int *n_compounds, char *spec_file);
void read_spec_file(int n_particles, int n_compounds, char *spec_file, 
		    double *A, double *N, double *K, double *G, 
		    double solvent_number_density, double kT, 
		    int free_energy_input, int dimensionless, 
		    int write_log_file, char *eq_file, char *log_file);
void read_con_file(int n_compounds, char *con_file, double *x0, 
                   double solvent_number_density, int dimensionless);
/* ************************************************************************** */

/* ************************** IN MT19937AR.C ******************************** */
// Random number generation by Mersenne Twister
void init_genrand(unsigned long s);
double genrand_real1(void);
/* ************************************************************************** */

#ifdef __cplusplus
}
#endif
