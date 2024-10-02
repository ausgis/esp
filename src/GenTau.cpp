#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericVector GenTau(double n) {
  // Generate an Armadillo vector with n points between 0 and 1
  arma::vec partition = arma::linspace(1.0/n, 1.0, n);

  // Convert arma::vec to Rcpp::NumericVector and return
  return Rcpp::NumericVector(partition.begin(), partition.end());
}
