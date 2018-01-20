#include "stdio.h"
#include "RcppArmadillo.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Invert a column vector
//'
//' @param x a vector
// [[Rcpp::export]]
arma::colvec sinvc(arma::colvec x) {
  x.transform( [](double val) { return 1/val; });
  x.elem(find_nonfinite(x)).ones();
  return x;
}

//' Invert a row vector
//'
//' @param x a vector
// [[Rcpp::export]]
arma::rowvec sinvr(arma::rowvec x) {
  x.transform( [](double val) { return 1/val; });
  x.elem(find_nonfinite(x)).ones();
  return x;
}
