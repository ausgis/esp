#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericVector QuantileDisc(const arma::vec& x, int n) {
  // Ensure n > 0
  if (n <= 0) {
    Rcpp::stop("n must be greater than 0");
  }

  // Get the size of the vector
  int size = x.n_elem;

  // Create a NumericVector to store the results
  Rcpp::NumericVector result(size);

  // Calculate quantile cut points
  arma::vec quantiles = arma::quantile(x, arma::linspace<arma::vec>(0.0, 1.0, n + 1));

  // Iterate over the original vector and classify based on quantiles
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < n; ++j) {
      if (x[i] <= quantiles[j + 1]) {
        result[i] = j + 1;  // Classification starts from 1
        break;
      }
    }
  }

  return result;  // Return classification results
}
