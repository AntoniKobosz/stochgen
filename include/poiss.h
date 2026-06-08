#ifndef POISS
#define POISS
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


void poiss_jump_times(double* t,double lambda,size_t n,gsl_rng* rng);


#endif
