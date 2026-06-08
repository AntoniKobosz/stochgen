#ifndef GEOMBROWNIAN
#define GEOMBROWNIAN

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void geom_brownian_vector(double* x,double x0, double mu, double sigma, double dt,size_t n,gsl_rng* rng);


void geom_brownian_matrix(double* A,double A0, double mu, double sigma, double dt, 
    size_t n,size_t m,gsl_rng* rng);

void geom_brownian_end_info(double* xT,int* crossed_lower,int* crossed_upper, 
    double x0, double mu, double sigma, double dt, 
    double lower_thresh, double upper_thresh, size_t n, size_t m, gsl_rng* rng);


#endif