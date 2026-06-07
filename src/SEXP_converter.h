#ifndef SEXP_CONVERTER.H
#define SEXP_CONVERTER.H

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

size_t as_size_t(SEXP n);

double as_double(SEXP x);

double as_positive_double(SEXP x);

#endif