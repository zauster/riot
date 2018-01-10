#ifndef HELPER
#define HELPER

#include "RcppArmadillo.h"

using namespace Rcpp;

arma::colvec sinvc(arma::colvec x);

arma::rowvec sinvr(arma::rowvec x);

#endif
