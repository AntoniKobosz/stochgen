#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include "stochastic_generator.h"
#include <R.h>
#include <Rinternals.h>
#include <omp.h>


// Brownian increment function
inline double dx_brownian(double x, process_params p, double dt, double dW)
{
    return p.mu * dt + p.sigma * dW;
}
// Geometric Brownian increment function
inline double dx_geom_brownian(double x, process_params p, double dt, double dW)
{
    return x * (p.mu * dt + p.sigma * dW);
}
// Ornstein Ohlenbeck increment function
inline double dx_ornstein_uhlenbeck(double x, process_params p, double dt, double dW)
{
    return -p.theta*(x-p.mu)*dt + p.sigma * dW;
}

// Returns one of the above increment functions according to the process type
dx_fn select_process(process_type type)
{
    switch (type)
    {
    case PROCESS_BROWNIAN:
        return dx_brownian;

    case PROCESS_GEOM_BROWNIAN:
        return dx_geom_brownian;

    case PROCESS_ORNSTEIN_UHLENBECK:
        return dx_ornstein_uhlenbeck;

    default:
        Rf_error("No such process exists :(");
    }
}


void stochastic_vector(double* x,double x0,process_params p, double dt,
     size_t n,dx_fn dx,gsl_rng* rng)
{
    // Generates a stochastic process with initial value x0, 
    // by iteratively adding dx to the previous value.
    // dW is the increment of a standard Weiner process over a time step dt

    // resultes are saved into the array x of length n
    x[0] = x0;
    double sqrt_dt = sqrt(dt);

    for(size_t i = 1; i < n; ++i){
        double dW = gsl_ran_gaussian(rng,sqrt_dt);
        x[i] = x[i-1] + dx(x[i-1],p,dt,dW);
    }
}

void stochastic_matrix(double* A,double x0, process_params p, double dt, 
    size_t n,size_t m,dx_fn dx,gsl_rng** rng_array)
{
    // Generates a matrix of m syochastic vectors, each of length n.
    // writes result into an nxm matrix in Fortran order
    double sqrt_dt = sqrt(dt);
    
    #pragma omp parallel for
    for(size_t j = 0; j < m; ++j){
        gsl_rng* rng = rng_array[omp_get_thread_num()];
        A[0+n*j] = x0;
        for(size_t i = 1; i < n; ++i){
                double dW = gsl_ran_gaussian(rng,sqrt_dt);
                double x = A[i-1+n*j];
                A[i+n*j] = x + dx(x,p,dt,dW);
        }
    }
}

void stochastic_end_info(double* xT,int* crossed_lower,int* crossed_upper, 
    double x0, process_params p, double dt, double lower_thresh,
    double upper_thresh, size_t n, size_t m, dx_fn dx, gsl_rng** rng_array)
{
    // Simulates m independent stochastic processes, 
    // returnes their values at the final time T = n*dt and 
    // whether they went below/above the provided thresholds during their runtime
    double sqrt_dt = sqrt(dt);

    #pragma omp parallel for
    for(size_t j=0; j < m; j++){
        gsl_rng* rng = rng_array[omp_get_thread_num()];
        double x = x0;
        int low_flag = 0, upper_flag = 0;
        for(size_t i = 1; i < n; i++){
            double dW = gsl_ran_gaussian(rng,sqrt_dt);
            x += dx(x,p,dt,dW);
            if(x >= upper_thresh) upper_flag = 1;
            if(x <= lower_thresh) low_flag = 1;
        }
        xT[j] = x;
        crossed_lower[j] = low_flag;
        crossed_upper[j] = upper_flag;
    }
}