#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include "poiss.h"


typedef unsigned long long ull;


SEXP poiss_jump_times_r(SEXP lamba_,SEXP n_,SEXP seed_)
{
    // Tymczasowo ręcznie zamień typy
    double lambda = Rf_asReal(lamba_);
    int n = Rf_asInteger(n);
    ull seed =(ull) Rf_asReal(seed_);

    // Stwórz wynikowy SEXP
    SEXP T;
    PROTECT(T = Rf_allocVector(REALSXP,n));
    double* t = REAL(T);
    poiss_jump_times(t,lambda,n,seed);

    UNPROTECT(1);
}   