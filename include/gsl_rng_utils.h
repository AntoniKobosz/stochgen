#ifndef GSLRNGUTILS
#define GSLRNGUTILS

gsl_rng* create_rng();

gsl_rng** create_rng_array(int n);

void free_rng_array(gsl_rng** rng_array,int n);

#endif