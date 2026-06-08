#ifndef SEXP_CONVERTER
#define SEXP_CONVERTER

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

size_t as_size_t(SEXP n);

double as_double(SEXP x);

double as_positive_double(SEXP x);

#endif