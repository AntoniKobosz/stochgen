#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

size_t as_size_t(SEXP n)
{
    if(!Rf_isInteger(n))
        Rf_error("Expected an integer");

    if(Rf_length(n) != 1)
        Rf_error("Expected scalar, got vector of length %d",Rf_length(n));
    
    int val = Rf_asInteger(n);

    if(val == NA_INTEGER)
        Rf_error("Na not allowed");

    if(val < 0)
        Rf_error("expected non negative integer");
    
    return (size_t) val;
}


double as_double(SEXP x)
{
    if(!Rf_isNumeric(x))
        Rf_error("Expected numeric");

    if(Rf_length(x) != 1)
        Rf_error("Expected scalar, got vector of length %d",Rf_length(x));
   
    double val = Rf_asReal(x);

    if(ISNA(val))
        Rf_error("Na not allowed");

    return val;
}

double as_positive_double(SEXP x)
{
    double val = as_double(x);
    if(val <= 0)
        Rf_error("value must be positive");
    return val;
}