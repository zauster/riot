#include "stdio.h"
#include "RcppArmadillo.h"

using namespace std;

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

//' Do the GRAS algorithm for updating matrizes
//'
//' This function calculates an updated matrix X, based on a matrix A,
//' which meets the given row and columns totals. Note: A can contain
//' negative elements.
//' @param A the "base" matrix
//' @param u vector with the row totals
//' @param v vector with the column totals
//' @param epsilon the error tolerance level, default is 1e-10.
//' @param max.iter maximum number of iterations, default is 10000.
//' @param verbose should some information of the iterations be
//'     displayed? Default is FALSE
//' @return the updated matrix X
//' @author Oliver Reiter
//' @references Junius T. and J. Oosterhaven (2003), The solution of
//'     updating or regionalizing a matrix with both positive and
//'     negative entries, Economic Systems Research, 15, pp. 87-96.
//'
//'     Lenzen M., R. Wood and B. Gallego (2007), Some comments on the
//'     GRAS method, Economic Systems Research, 19, pp. 461-465.
//'
//'     Temurshoev, U., R.E. Miller and M.C. Bouwmeester (2013), A note
//'     on the GRAS method, Economic Systems Research, 25, pp. 361-367.
//' @keywords gras, matrix updating
//' @examples
//' ## example from the papers
//' A <- matrix(c(7, 3, 5, -3, 2, 9, 8, 0, -2, 0, 2, 0),
//'             ncol = 4, nrow = 3, byrow = TRUE)
//' u <- c(15, 26, -1)
//' v <- c(9, 16, 17, -2)
//' Rcpp_doGRAS(A, u, v)
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Rcpp_doGRAS(arma::mat A,
                      arma::colvec u,
                      arma::rowvec v,
                      double epsilon = 1e-10,
                      int maxiter = 1000,
                      bool verbose = true) {
  int nrowsA = A.n_rows;
  int ncolsA = A.n_cols;
  int nrowsu = u.n_rows;
  int ncolsv = v.n_cols;

  if(nrowsA != nrowsu) {
    // throw std::range_error("Unequal number of rows of A and u!");
    ::Rf_error("Unequal number of rows of A and u!");
  }
  if(ncolsA != ncolsv) {
    // throw std::range_error("Unequal number of cols of A and v!");
    ::Rf_error("Unequal number of cols of A and v!");
  }

  // separate A into a positive (P) and a negative (N) matrix
  arma::mat P = arma::mat(A);
  P.transform( [](double val) { return val >= 0? val : 0; });
  arma::mat N = arma::mat(A);
  N.transform( [](double val) { return val < 0? -val : 0; });

  // initialize the multiplicators
  arma::rowvec s = arma::rowvec(ncolsA).ones();
  arma::rowvec sold, salt, pjr, njr = arma::rowvec(ncolsA);

  arma::colvec r = arma::colvec(nrowsA).ones();
  arma::colvec ralt, pis, nis = arma::colvec(nrowsA);

  double error = 1;
  int iter = 1;
  int errorpos;

  do {
    sold = s;

    // row multiplicator s
    pjr = r.t() * P;
    njr = sinvc(r).t() * N;
    s = sinvr(2 * pjr) % (v + arma::sqrt(pow(v, 2) + 4 * pjr % njr));
    salt = -1 * sinvr(v) % njr;
    s(find(s == 0)) = salt(find(s == 0));
    // cout << "s: " << s; // << endl;

    // col multiplicator r
    pis = P * s.t();
    nis = N * sinvr(s).t();
    r = sinvc(2 * pis) % (u + arma::sqrt(pow(u, 2) + 4 * pis % nis));
    ralt = -1 * sinvc(u) % nis;
    r(find(r == 0)) = ralt(find(r == 0));
    // cout << "r': " << r.t(); // << endl;

    error = abs(sold - s).max();
    errorpos = abs(sold - s).index_max();

    if(verbose) {
      cout << "iter: " << iter << " -> error: " << error
           << " at " << errorpos << endl;
    }

    iter++;

  } while((error > epsilon) & (iter <= maxiter));

  if(iter >= maxiter) {
    cout << "Warning: Reached maxiter!" << endl;
  }

  if(verbose) {
    cout << " => Iterations needed: " << iter;
    cout << " => Resulting error: " << error << " at " << errorpos << endl;
  }

  P.each_row() %= s;
  P.each_col() %= r;
  N.each_row() %= sinvr(s);
  N.each_col() %= sinvc(r);

  return P - N;
}
