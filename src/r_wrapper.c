#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "SEXP_converter.h"
#include "poiss.h"
#include "brownian.h"


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

SEXP brownian_matrix_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_,
    SEXP n_, SEXP m_)
{
    double x0 = as_double(x0_);
    double mu = as_double(mu_);
    double sigma = as_positive_double(sigma_);
    double dt = as_positive_double(dt_);
    size_t n = as_size_t(n_);
    size_t m = as_size_t(m_);
    gsl_rng* rng = create_rng();

    SEXP M = PROTECT(Rf_allocMatrix(REALSXP,n,m));
    double* A = REAL(M);
    brownian_matrix(A,x0,mu,sigma,dt,n,m,rng);
    gsl_rng_free(rng);
    UNPROTECT(1);
    return M;
}

SEXP brownian_end_info_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_,
    SEXP lower_thresh_, SEXP upper_thresh_, SEXP n_, SEXP m_)
{
    // Convert types to C
    double x0 = as_double(x0_);
    double mu = as_double(mu_);
    double sigma = as_positive_double(sigma_);
    double dt = as_positive_double(dt_);
    double lower_thresh = as_double(lower_thresh_);
    double upper_thresh = as_double(upper_thresh_);
    size_t n = as_size_t(n_);
    size_t m = as_size_t(m_);
    gsl_rng* rng = create_rng();

    // Create list with results
    SEXP out = PROTECT(Rf_allocVector(VECSXP,3));
    SEXP xT_ = PROTECT(Rf_allocVector(REALSXP,n));
    SEXP crossed_lower_ = PROTECT(Rf_allocVector(LGLSXP,n));
    SEXP crossed_upper_ = PROTECT(Rf_allocVector(LGLSXP,n));
    // Extract arrays
    double* xT = REAL(xT_);
    int* crossed_lower = LOGICAL(crossed_lower_);
    int* crossed_upper = LOGICAL(crossed_upper_);
    // Fill vectors with result
    brownian_end_info(xT,crossed_lower,crossed_upper,x0,mu,sigma,dt,
    lower_thresh,upper_thresh,n,m,rng);
    // Attach vectors to the list
    SET_VECTOR_ELT(out,0,xT_);
    SET_VECTOR_ELT(out,1,crossed_lower_);
    SET_VECTOR_ELT(out,2,crossed_upper_);
    // Names
    SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
    SET_VECTOR_ELT(names,0,Rf_mkChar("xT"));
    SET_VECTOR_ELT(names,0,Rf_mkChar("crossed_lower"));
    SET_VECTOR_ELT(names,0,Rf_mkChar("crossed_upper"));
    Rf_setAttrib(out,R_NamesSymbol,names);
    
    UNPROTECT(5);
    gsl_rng_free(rng);
    return out;
}