#include "stdio.h"
#include "RcppArmadillo.h"

using namespace Rcpp;

//' Update the inverse by column
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void ShermanUpdateCol(arma::mat &A, arma::colvec u, int i) {
  A -= ((A * u) * A.row(i)) / arma::as_scalar(1 + A.row(i) * u);
}


//' Update the inverse by row
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void ShermanUpdateRow(arma::mat &A, arma::rowvec v, int i) {
  A -= (A.col(i) * (v * A)) / arma::as_scalar(1 + v * A.col(i));
}

//' Return a vector of indices except index i
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec index_noti(int i, int nm1) {
  arma::uvec index = arma::uvec(nm1);
  int pos = 0;
  for(int j = 0; j < nm1 + 1; j++) {
    if(j != i) {
      index(pos) = j;
      pos++;
    }
  }
  return index;
}

//' Calculate random walk centrality
//'
//' Calculate random walk centrality of an input-output matrix
//' @param A the input-output matrix
//' @param verbose should some information of the iterations be
//'     displayed? Default is FALSE
//' @return a vector of centralities
//' @author Oliver Reiter
//' @references Blöchl F, Theis FJ, Vega-Redondo F, and Fisher E: Vertex Centralities in Input-Output Networks Reveal the Structure of Modern Economies, Physical Review E, 83(4):046127, 2011
//' @keywords input-output analysis
//' @examples
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::rowvec rwCentrality(arma::mat A,
                          bool verbose = false) {

  int n = A.n_rows;
  int nm1 = A.n_rows - 1; // nm1 = n minus 1
  arma::colvec u = arma::colvec(nm1);
  arma::rowvec v = arma::rowvec(nm1);
  arma::mat IMinv = arma::mat(nm1, nm1);
  arma::mat H = arma::mat(A).zeros();

  arma::uvec index_i = arma::uvec(1);
  arma::uvec index_ip1 = arma::uvec(1); // ip1 = i plus 1

  // calculate the transformation matrix M
  arma::colvec tmpSum = sum(A, 1);
  A.each_col() /= tmpSum;
  A = eye(size(A)) - A;

  // calculate the first inverse
  IMinv = A.submat(1, 1, nm1, nm1).i();

  for(int i = 0; i < n; i++) {

    index_i(0) = i;
    index_ip1(0) = i + 1;

    // calculate matrix H (as row sums of the inverse(IM))
    H.submat(index_noti(i, nm1), index_i) = sum(IMinv, 1);
    // Rcout << "H:\n" << H << std::endl;

    if(i < nm1) {
      // column update
      u = A.submat(index_noti(i + 1, nm1), index_i) -
        A.submat(index_noti(i, nm1), index_ip1);
      u(i) /= 2;
      // Rcout << "u:\n" << u << std::endl;
      ShermanUpdateCol(IMinv, u, i);

      // row update
      v = A.submat(index_i, index_noti(i + 1, nm1)) -
        A.submat(index_ip1, index_noti(i, nm1));
      v(i) /= 2;
      // Rcout << "v:\n" << v << std::endl;
      ShermanUpdateRow(IMinv, v, i);

      // if some elements are not finite after the updates,
      // calculate a new inverse from scratch
      if(!IMinv.is_finite()) {
        Rcout << "Calc new IMinv\n" << IMinv << std::endl;
        IMinv = A.submat(index_noti(i + 1, nm1), index_noti(i + 1, nm1)).i();
      }
      // Rcout << "IMinv upd:\n" << IMinv << std::endl;
    }
  }
  arma::rowvec res = sum(H, 0);
  return n / res;
}
