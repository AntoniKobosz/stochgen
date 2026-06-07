#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include "SEXP_converter.h"
#include "poiss.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng* create_rng()
{
    // Creates gsl_rng* object consistant with
    // current R seed
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

    GetRNGstate();
    unsigned long seed = (unsigned long) (unif_rand() * 0x100000000UL); // 2^32
    PutRNGstate();

    gsl_rng_set(rng,seed);
    return rng;
}

SEXP poiss_jump_times_r(SEXP lamba_,SEXP n_)
{
    double lambda = as_double(lamba_);
    size_t n = as_size_t(n_);
    gsl_rng* rng = create_rng();
    // Stwórz wynikowy SEXP
    SEXP T;
    PROTECT(T = Rf_allocVector(REALSXP, n));
    double* t = REAL(T);
    poiss_jump_times(t, lambda, n, rng);

    UNPROTECT(1);
    return T;
}

SEXP brownian_vector_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_, SEXP n_)
{
    double x0 = as_double(x0_);
    double mu = as_double(mu_);
    double sigma = as_positive_double(sigma_);
    double dt = as_positive_double(dt_);
    size_t n = as_size_t(n_);
    gsl_rng* rng = create_rng();

    SEXP X = PROTECT(Rf_allocVector(REALSXP,n));
    double* x = REAL(X);
    brownian_vector(x,x0,mu,sigma,dt,n,rng);
    gsl_rng_free(rng);
    UNPROTECT(1);
    return X;
}