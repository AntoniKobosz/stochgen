#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "SEXP_converter.h"
#include "poiss.h"
#include "brownian.h"
#include "geom_brownian.h"
#include "stochastic_generator.h"
#include <omp.h>

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
gsl_rng** create_rng_array(int n)
{
    // Creates an array of gsl_rng* objects for multithreading
    gsl_rng** rng_array = malloc(n * sizeof(gsl_rng*));
    for (size_t i = 0; i < n; i++)
    {
        rng_array[i] = create_rng();
    }
    return rng_array;
}
// Frees the array and each gsl_rng element
void free_rng_array(gsl_rng** rng_array,int n)
{
    for (size_t i = 0; i < n; i++)
        gsl_rng_free(rng_array[i]);
    free(rng_array);    
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

    gsl_rng_free(rng);
    UNPROTECT(1);
    return T;
}

SEXP stochastic_vector_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP theta_, SEXP dt_, SEXP n_, SEXP process_type)
{
    double x0 = as_double(x0_);
    process_params p = {
        .mu = as_double(mu_),
        .sigma = as_positive_double(sigma_),
        .theta = as_positive_double(theta_)
    };
    double dt = as_positive_double(dt_);
    size_t n = as_size_t(n_);
    int type = (int) as_size_t(process_type);
    dx_fn dx = select_process(type);


    gsl_rng* rng = create_rng();

    SEXP X = PROTECT(Rf_allocVector(REALSXP,n));
    double* x = REAL(X);
    stochastic_vector(x,x0,p,dt,n,dx,rng);
    gsl_rng_free(rng);
    UNPROTECT(1);
    return X;
}

SEXP stochastic_matrix_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP theta_, SEXP dt_,
    SEXP n_, SEXP m_, SEXP process_type)
{
    double x0 = as_double(x0_);
    process_params p = {
        .mu = as_double(mu_),
        .sigma = as_positive_double(sigma_),
        .theta = as_positive_double(theta_)
    };
    double dt = as_positive_double(dt_);
    size_t n = as_size_t(n_);
    size_t m = as_size_t(m_);
    int type = (int) as_size_t(process_type);
    dx_fn dx = select_process(type);

    int thread_count = omp_get_max_threads();
    gsl_rng** rng_array = create_rng_array(thread_count);

    SEXP M = PROTECT(Rf_allocMatrix(REALSXP,n,m));
    double* A = REAL(M);
    stochastic_matrix(A,x0,p,dt,n,m,dx,rng_array);
    free_rng_array(rng_array,thread_count);
    UNPROTECT(1);
    return M;
}

SEXP stochastic_end_info_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP theta_, SEXP dt_,
    SEXP lower_thresh_, SEXP upper_thresh_, SEXP n_, SEXP m_, SEXP process_type)
{
    // Convert types to C
    double x0 = as_double(x0_);
    process_params p = {
        .mu = as_double(mu_),
        .sigma = as_positive_double(sigma_),
        .theta = as_positive_double(theta_)
    };
    double dt = as_positive_double(dt_);
    double lower_thresh = as_double(lower_thresh_);
    double upper_thresh = as_double(upper_thresh_);
    size_t n = as_size_t(n_);
    size_t m = as_size_t(m_);
    int type = (int) as_size_t(process_type);
    dx_fn dx = select_process(type);

    int thread_count = omp_get_max_threads();
    gsl_rng** rng_array = create_rng_array(thread_count);

    // Create list with results
    SEXP out = PROTECT(Rf_allocVector(VECSXP,3));
    SEXP xT_ = PROTECT(Rf_allocVector(REALSXP,m));
    SEXP crossed_lower_ = PROTECT(Rf_allocVector(LGLSXP,m));
    SEXP crossed_upper_ = PROTECT(Rf_allocVector(LGLSXP,m));
    // Extract arrays
    double* xT = REAL(xT_);
    int* crossed_lower = LOGICAL(crossed_lower_);
    int* crossed_upper = LOGICAL(crossed_upper_);
    // Fill vectors with result
    stochastic_end_info(xT, crossed_lower, crossed_upper, x0, p, dt,
    lower_thresh, upper_thresh, n, m, dx, rng_array);
    // Attach vectors to the list
    SET_VECTOR_ELT(out,0,xT_);
    SET_VECTOR_ELT(out,1,crossed_lower_);
    SET_VECTOR_ELT(out,2,crossed_upper_);
    // Names
    SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
    SET_STRING_ELT(names,0,Rf_mkChar("xT"));
    SET_STRING_ELT(names,1,Rf_mkChar("crossed_lower"));
    SET_STRING_ELT(names,2,Rf_mkChar("crossed_upper"));
    Rf_setAttrib(out,R_NamesSymbol,names);
    
    UNPROTECT(5);
    free_rng_array(rng_array,thread_count);
    return out;
}

























// SEXP brownian_vector_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_, SEXP n_)
// {
//     double x0 = as_double(x0_);
//     double mu = as_double(mu_);
//     double sigma = as_positive_double(sigma_);
//     double dt = as_positive_double(dt_);
//     size_t n = as_size_t(n_);
//     gsl_rng* rng = create_rng();

//     SEXP X = PROTECT(Rf_allocVector(REALSXP,n));
//     double* x = REAL(X);
//     brownian_vector(x,x0,mu,sigma,dt,n,rng);
//     gsl_rng_free(rng);
//     UNPROTECT(1);
//     return X;
// }

// SEXP brownian_matrix_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_,
//     SEXP n_steps_, SEXP n_traj_)
// {
//     double x0 = as_double(x0_);
//     double mu = as_double(mu_);
//     double sigma = as_positive_double(sigma_);
//     double dt = as_positive_double(dt_);
//     size_t n_steps = as_size_t(n_steps_);
//     size_t n_traj = as_size_t(n_traj_);
//     gsl_rng* rng = create_rng();

//     SEXP M = PROTECT(Rf_allocMatrix(REALSXP,n_steps,n_traj));
//     double* A = REAL(M);
//     brownian_matrix(A,x0,mu,sigma,dt,n_steps,n_traj,rng);
//     gsl_rng_free(rng);
//     UNPROTECT(1);
//     return M;
// }

// SEXP brownian_end_info_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_,
//     SEXP lower_thresh_, SEXP upper_thresh_, SEXP n_, SEXP m_)
// {
//     // Convert types to C
//     double x0 = as_double(x0_);
//     double mu = as_double(mu_);
//     double sigma = as_positive_double(sigma_);
//     double dt = as_positive_double(dt_);
//     double lower_thresh = as_double(lower_thresh_);
//     double upper_thresh = as_double(upper_thresh_);
//     size_t n = as_size_t(n_);
//     size_t m = as_size_t(m_);
//     gsl_rng* rng = create_rng();

//     // Create list with results
//     SEXP out = PROTECT(Rf_allocVector(VECSXP,3));
//     SEXP xT_ = PROTECT(Rf_allocVector(REALSXP,m));
//     SEXP crossed_lower_ = PROTECT(Rf_allocVector(LGLSXP,m));
//     SEXP crossed_upper_ = PROTECT(Rf_allocVector(LGLSXP,m));
//     // Extract arrays
//     double* xT = REAL(xT_);
//     int* crossed_lower = LOGICAL(crossed_lower_);
//     int* crossed_upper = LOGICAL(crossed_upper_);
//     // Fill vectors with result
//     brownian_end_info(xT,crossed_lower,crossed_upper,x0,mu,sigma,dt,
//     lower_thresh,upper_thresh,n,m,rng);
//     // Attach vectors to the list
//     SET_VECTOR_ELT(out,0,xT_);
//     SET_VECTOR_ELT(out,1,crossed_lower_);
//     SET_VECTOR_ELT(out,2,crossed_upper_);
//     // Names
//     SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
//     SET_STRING_ELT(names,0,Rf_mkChar("xT"));
//     SET_STRING_ELT(names,1,Rf_mkChar("crossed_lower"));
//     SET_STRING_ELT(names,2,Rf_mkChar("crossed_upper"));
//     Rf_setAttrib(out,R_NamesSymbol,names);
    
//     UNPROTECT(5);
//     gsl_rng_free(rng);
//     return out;
// }

// SEXP geom_brownian_vector_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_, SEXP n_)
// {
//     double x0 = as_double(x0_);
//     double mu = as_double(mu_);
//     double sigma = as_positive_double(sigma_);
//     double dt = as_positive_double(dt_);
//     size_t n = as_size_t(n_);
//     gsl_rng* rng = create_rng();

//     SEXP X = PROTECT(Rf_allocVector(REALSXP,n));
//     double* x = REAL(X);
//     geom_brownian_vector(x,x0,mu,sigma,dt,n,rng);
//     gsl_rng_free(rng);
//     UNPROTECT(1);
//     return X;
// }

// SEXP geom_brownian_matrix_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_,
//     SEXP n_, SEXP m_)
// {
//     double x0 = as_double(x0_);
//     double mu = as_double(mu_);
//     double sigma = as_positive_double(sigma_);
//     double dt = as_positive_double(dt_);
//     size_t n = as_size_t(n_);
//     size_t m = as_size_t(m_);
//     gsl_rng* rng = create_rng();

//     SEXP M = PROTECT(Rf_allocMatrix(REALSXP,n,m));
//     double* A = REAL(M);
//     geom_brownian_matrix(A,x0,mu,sigma,dt,n,m,rng);
//     gsl_rng_free(rng);
//     UNPROTECT(1);
//     return M;
// }

// SEXP geom_brownian_end_info_r(SEXP x0_, SEXP mu_, SEXP sigma_, SEXP dt_,
//     SEXP lower_thresh_, SEXP upper_thresh_, SEXP n_, SEXP m_)
// {
//     // Convert types to C
//     double x0 = as_double(x0_);
//     double mu = as_double(mu_);
//     double sigma = as_positive_double(sigma_);
//     double dt = as_positive_double(dt_);
//     double lower_thresh = as_double(lower_thresh_);
//     double upper_thresh = as_double(upper_thresh_);
//     size_t n = as_size_t(n_);
//     size_t m = as_size_t(m_);
//     gsl_rng* rng = create_rng();

//     // Create list with results
//     SEXP out = PROTECT(Rf_allocVector(VECSXP,3));
//     SEXP xT_ = PROTECT(Rf_allocVector(REALSXP,m));
//     SEXP crossed_lower_ = PROTECT(Rf_allocVector(LGLSXP,m));
//     SEXP crossed_upper_ = PROTECT(Rf_allocVector(LGLSXP,m));
//     // Extract arrays
//     double* xT = REAL(xT_);
//     int* crossed_lower = LOGICAL(crossed_lower_);
//     int* crossed_upper = LOGICAL(crossed_upper_);
//     // Fill vectors with result
//     geom_brownian_end_info(xT,crossed_lower,crossed_upper,x0,mu,sigma,dt,
//     lower_thresh,upper_thresh,n,m,rng);
//     // Attach vectors to the list
//     SET_VECTOR_ELT(out,0,xT_);
//     SET_VECTOR_ELT(out,1,crossed_lower_);
//     SET_VECTOR_ELT(out,2,crossed_upper_);
//     // Names
//     SEXP names = PROTECT(Rf_allocVector(STRSXP,3));
//     SET_STRING_ELT(names,0,Rf_mkChar("xT"));
//     SET_STRING_ELT(names,1,Rf_mkChar("crossed_lower"));
//     SET_STRING_ELT(names,2,Rf_mkChar("crossed_upper"));
//     Rf_setAttrib(out,R_NamesSymbol,names);
    
//     UNPROTECT(5);
//     gsl_rng_free(rng);
//     return out;
// }