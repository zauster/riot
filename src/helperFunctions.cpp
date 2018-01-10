#include "stdio.h"
#include "RcppArmadillo.h"

using namespace Rcpp;

//' Invert a column vector
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec sinvc(arma::colvec x) {
  x.transform( [](double val) { return 1/val; });
  x.elem(find_nonfinite(x)).ones();
  return x;
}

//' Invert a row vector
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::rowvec sinvr(arma::rowvec x) {
  x.transform( [](double val) { return 1/val; });
  x.elem(find_nonfinite(x)).ones();
  return x;
}

