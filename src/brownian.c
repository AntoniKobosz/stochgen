#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>

void brownian_vector(double* x,double x0, double mu, double sigma, double dt,
     size_t n,gsl_rng* rng)
{
    // Generates a stochastic process given by 
    // dx = mu*dt + sigma*dW
    // With initial value x0, 
    // where dW is the increment of a standard Weiner process over a time step dt

    // resultes are saved into the array x of length n
    x[0] = x0;
    double sqrt_dt = sqrt(dt);
    for(size_t i = 1; i < n; ++i){
        double dW = gsl_ran_gaussian(rng,sqrt_dt);
        x[i] = x[i-1] + mu*dt + sigma*dW;
    }
}
void brownian_matrix(double* A,double x0, double mu, double sigma, double dt, 
    size_t n,size_t m,gsl_rng* rng)
{
    // Generates a matrix of m brownian vectors, each of length n.
    // writes result into an nxm matrix in Fortran order
    double sqrt_dt = sqrt(dt);
    double drift = mu*dt;
    
    for(size_t j = 0; j < m; ++j){
        A[0+n*j] = x0;
        for(size_t i = 1; i < n; ++i){
                double dW = gsl_ran_gaussian(rng,sqrt_dt);
                A[i+n*j] = A[i-1+n*j] + drift + sigma*dW;
        }
    }
}
void brownian_end_info(double* xT,int* crossed_lower,int* crossed_upper, 
    double x0, double mu, double sigma, double dt, 
    double lower_thresh, double upper_thresh, size_t n, size_t m, gsl_rng* rng)
{
    // Simulates m independent brownian processes, 
    // returnes their values at the final time T = n*dt and 
    // whether they went below/above the provided thresholds during their runtime
    double sqrt_dt = sqrt(dt);
    double drift = mu*dt;

    for(size_t j=0; j < m; j++){
        double x = x0;
        int low_flag = 0, upper_flag = 0;
        for(size_t i = 1; i < n; i++){
            double dW = gsl_ran_gaussian(rng,sqrt_dt);
            x += drift + sigma * dW;
            if(x >= upper_thresh) upper_flag = 1;
            if(x <= lower_thresh) low_flag = 1;
        }
        xT[j] = x;
        crossed_lower[j] = low_flag;
        crossed_upper[j] = upper_flag;
    }
}