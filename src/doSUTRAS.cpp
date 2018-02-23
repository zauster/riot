#include "stdio.h"
#include "RcppArmadillo.h"
#include "helperFunctions.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Do the SUTRAS algorithm for updating supply and use tables
//'
//' This function calculates an updated supply and (domestic and import) use
//' table (SUTs), based on a given SUT and national accounts data.
//' @param V the supply (make) matrix
//' @param Ud the domestic use matrix
//' @param Um the import use matrix
//' @param m0 the import vector
//' @param u_bar vector with the new intermediate input sums
//' @param x_bar vector with the new gross output sums
//' @param M scalar, new total imports
//' @param c vector with the negative entries of the trade margins, given in absolute values
//' @param epsilon the error tolerance level, default is 1e-10.
//' @param maxiter maximum number of iterations, default is 10000.
//' @param verbose should some information of the iterations be
//'     displayed? Default is FALSE
//' @return the updated SUTs, a list
//' @author Oliver Reiter
//' @references Temurshoev, U., & Timmer, M. P. (2011). Joint estimation of supply and use tables. Papers in Regional Science, 90(4), 863-882.
//' @keywords SUTRAS, input-output updating
//' @examples
//' data(TestSUTs_AUT2010)
//' data(TestNA_AUT2011)
//' c <- rep(0, length(m0))
//' res <- doSUTRAS(V = V, Ud = Ud, Um = Um, m0 = m0, u_bar = ubar,
//'                 x_bar = xbar, M = M, c = c, epsilon = 1e-10,
//'                 maxiter = 100, verbose = FALSE)
//' @export
// [[Rcpp::export]]
Rcpp::List doSUTRAS(arma::mat V, arma::mat Ud,
                    arma::mat Um, arma::vec m0,
                    arma::rowvec u_bar,
                    arma::rowvec x_bar,
                    double M, arma::colvec c,
                    double epsilon = 1e-10,
                    int maxiter = 1000,
                    bool verbose = false) {

  // separate supply matrix into positive and negative part
  arma::mat Pv0 = arma::mat(V);
  Pv0.transform( [](double val) { return val >= 0? val : 0; });
  arma::mat Nv0 = arma::mat(V);
  Nv0.transform( [](double val) { return val < 0? -val : 0; });

  // and the domestic use matrix too
  arma::mat Pd0 = arma::mat(Ud);
  Pd0.transform( [](double val) { return val >= 0? val : 0; });
  arma::mat Nd0 = arma::mat(Ud);
  Nd0.transform( [](double val) { return val < 0? -val : 0; });

  // and the import use matrix too
  arma::mat Pm0 = arma::mat(Um);
  Pm0.transform( [](double val) { return val >= 0? val : 0; });
  arma::mat Nm0 = arma::mat(Um);
  Nm0.transform( [](double val) { return val < 0? -val : 0; });

  // initialize the multiplicators
  double r = 1;
  arma::colvec rm = arma::colvec(Pm0.n_rows).ones();
  arma::colvec rm_m1 = arma::colvec(Pm0.n_rows);
  arma::colvec rd = arma::colvec(Pd0.n_rows).ones();
  arma::colvec rd_m1 = arma::colvec(Pd0.n_rows);
  arma::rowvec su = arma::rowvec(Pd0.n_cols).ones();
  arma::rowvec rv = arma::rowvec(Pv0.n_cols).ones();
  arma::rowvec rv_tmp = arma::rowvec(Pv0.n_cols);

  arma::colvec pd = arma::colvec(Pd0.n_rows);
  arma::colvec nd = arma::colvec(Pd0.n_rows);
  arma::rowvec ps = arma::rowvec(Pd0.n_cols);
  arma::rowvec ns = arma::rowvec(Pd0.n_cols);

  arma::colvec rd_diff = arma::colvec(Pd0.n_rows);
  arma::colvec rm_diff = arma::colvec(Pm0.n_rows);

  bool rd_test, rm_test;

  int iter = 1;

  do {
    rm_m1 = rm;
    rd_m1 = rd;

    // rd
    pd = Pd0 * su.t() + sum(Nv0.each_row() % sinvr(rv), 1);
    nd = sum(Nd0.each_row() % sinvr(su), 1) + Pv0 * rv.t();
    rd = (0.5 * sinvc(pd)) % (-1 * c + arma::sqrt((c % c) + 4 * (pd % nd)));

    // rm
    rm = arma::sqrt(sinvc(Pm0 * su.t()) % (sum(Nm0.each_row() % sinvr(su), 1) + r * m0));

    // rv
    rv_tmp = sum(Pv0.each_col() % sinvc(rd), 0);
    rv = x_bar + arma::sqrt(x_bar % x_bar + 4 * (rv % (Nv0.t() * rd).t()));
    rv = 0.5 * sinvr(rv_tmp) % rv;

    // su
    ps = (Pd0.t() * rd + Pm0.t() * rm).t();
    ns = sum(Nd0.each_col() % sinvc(rd), 0) + sum(Nm0.each_col() % sinvc(rm), 0);
    ns = u_bar + arma::sqrt(u_bar % u_bar + 4 * (ps % ns));
    su = 0.5 * sinvr(ps) % ns;

    // r
    r = M / sum(m0 % sinvc(rm));

    rd_test = any(abs(rd_m1 - rd) > epsilon);
    rm_test = any(abs(rm_m1 - rm) > epsilon);

    if(verbose) {
      rd_diff = abs(rd_m1 - rd);
      rm_diff = abs(rm_m1 - rm);
      Rcout << "Iter " << iter;
      Rcout << "  rd diff: " << sum(rd_diff)
            << "  rm diff: " << sum(rm_diff) << std::endl;
    } //else if(iter %% 50 == ) {}

    iter++;
  } while ((rd_test | rm_test) & (iter < maxiter));


  // apply the multiplicators, update the matrizes. We save the results in
  // the positive parts (eg. Pv0) of the respective matrizes
  // V, supply matrix:
  Pv0.each_col() %= sinvc(rd);
  Pv0.each_row() %= rv;
  Nv0.each_col() %= rd;
  Nv0.each_row() %= sinvr(rv);
  Pv0 = Pv0 - Nv0;


  // Use matrizes
  // Domestic Use
  Pd0.each_row() %= su;
  Pd0.each_col() %= rd;
  Nd0.each_row() %= sinvr(su);
  Nd0.each_col() %= sinvc(rd);
  Pd0 = Pd0 - Nd0;

  // Import Use
  Pm0.each_row() %= su;
  Pm0.each_col() %= rm;
  Nm0.each_row() %= sinvr(su);
  Nm0.each_col() %= sinvc(rm);
  Pm0 = Pm0 - Nm0;

  // Import vector
  m0 = r * (m0 % sinvc(rm));

  // test if the column sums are equal to the specified values (ie
  // ubar, xbar and M)
  rv = abs(sum(Pv0, 0) - x_bar);
  su = abs(sum(Pd0, 0) + sum(Pm0, 0) - u_bar);
  r = sum(m0) - M;
  if(any(rv > epsilon) | any(su > epsilon) | (r > epsilon)) {
    Rcout << "  M: " << r << std::endl;
    Rcout << "  xbar diff: " << rv << std::endl;
    Rcout << "  ubar diff: " << su << std::endl;
    Rcpp::warning("Some column sums are not completely equalized!");
  }

  // test if the product sums are equal to each other
  rm = sum(Pm0, 1) + sum(Pd0, 1);
  rd = sum(Pv0, 1) + m0;

  // product/row sums are never that exact as the column sums, that's why we
  // adjust the epsilon below.

  rm = abs(rm - rd);
  if(verbose) {
    Rcout << "Biggest diff in product sums: "
          << max(rm)
          << std::endl;
  }
  // double epsilon1000 = epsilon * 5000;
  // if(any(rm > epsilon1000)) {
  //   Rcout << "diff: " << rm.t() << std::endl;
  //   Rcpp::warning("Some product/row sums are not completely equalized!");
  // }

  return Rcpp::List::create(Rcpp::Named("iter") = iter,
                            Rcpp::Named("V")    = Pv0,
                            Rcpp::Named("m")    = m0,
                            Rcpp::Named("Ud")   = Pd0,
                            Rcpp::Named("Um")   = Pm0);
}
