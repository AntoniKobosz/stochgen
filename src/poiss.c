#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "poiss.h"

void poiss_jump_times(double* t,double lambda,size_t n,gsl_rng* rng){
    // Generuje momenty n skoków procesu Poissona Nt
    // Zapisuje wynik do wektora t

    // Odstępy czasowe między skokami mają rozkład Exp(lambda)
    double tau = 0;
    double mu = 1/lambda;

    for(size_t i = 0; i<n; i++){
        tau += gsl_ran_exponential(rng,mu);
        t[i] = tau;
    }
    gsl_rng_free(rng);
}
