#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



gsl_rng* create_rng()
{
    // Creates gsl_rng* object consistant with
    // current R seed
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

    GetRNGstate();
    unsigned long seed = (unsigned long) (unif_rand() * 0x100000000UL); // 2^32
    PutRNGstate(); // synchronise with R

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