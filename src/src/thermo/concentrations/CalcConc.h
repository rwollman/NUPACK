#ifndef NUPACK_THERMO_CONCENTRATIONS_CALCCONC_H__
#define NUPACK_THERMO_CONCENTRATIONS_CALCCONC_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


void getInitialGuess(double *x0, double *lambda, double *G, int **AT, int **A,
         int numSS, int numTotal, double PerturbScale, 
         unsigned long rand_seed);

int getx(double *x, double *lambda, double *G, int **AT, int numSS, int numTotal);

void getGrad(double *Grad, double *x0, double *x, int **A, int numSS, int numTotal);

void getHes(double **Hes, double *x, int **A, int numSS, int numTotal);

int getSearchDir(double *p, double *Grad, double **Hes, double delta, int numSS);

double getRho(double *lambda, double *p, double *Grad, double *x, double **Hes, 
        double *x0, double *G, int **AT, int numSS, int numTotal);

void getCauchyPoint(double *CauchyPoint, double **Hes, double *Grad, double delta, 
        int numSS);

void PerturbLambda(double *lambda, double PerturbScale, double *G, int **AT, 
       int numSS, int numTotal);

int CheckTol(double *Grad, double *AbsTol, int numSS);

int CalcConc(double *x, int **A, double *G, double *x0, int numSS, int numTotal, 
       int MaxIters, double tol, double deltaBar, double eta, double kT, 
       int MaxNoStep, int MaxTrial, double PerturbScale, int quiet, 
       int WriteLogFile, char *logFile, double MolesWaterPerLiter, 
       unsigned long seed);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUPACK_THERMO_CONCENTRATIONS_CALCCONC_H__ */
