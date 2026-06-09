#ifndef STOCHASTIC_GENERATOR
#define STOCHASTIC_GENERATOR


typedef enum {
    PROCESS_BROWNIAN = 0,
    PROCESS_GEOM_BROWNIAN = 1,
    PROCESS_ORNSTEIN_UHLENBECK = 2
} process_type;

typedef struct
{
    double mu;
    double sigma;
    double theta; // Used only by Ornstein–Uhlenbeck process
} process_params;

// template for a single-step increment function
typedef double (*dx_fn) (double x, process_params p, double dt, double dW);

double dx_brownian(double x, process_params p, double dt, double dW);

double dx_geom_brownian(double x, process_params p, double dt, double dW);

double dx_ornstein_uhlenbeck(double x, process_params p, double dt, double dW);

dx_fn select_process(process_type type);

void stochastic_vector(double* x,double x0,process_params p, double dt,
     size_t n,dx_fn dx,gsl_rng* rng);

void stochastic_matrix(double* A,double x0, process_params p, double dt, 
    size_t n,size_t m,dx_fn dx,gsl_rng** rng_array);

void stochastic_end_info(double* xT,int* crossed_lower,int* crossed_upper, 
    double x0, process_params p, double dt, double lower_thresh,
    double upper_thresh, size_t n, size_t m, dx_fn dx, gsl_rng** rng_array);

#endif