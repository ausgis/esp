#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function to find unique elements in a column
// [[Rcpp::export]]
arma::ivec ArmaRunique(const arma::ivec& x) {
  std::vector<int> seen;
  for (int i = 0; i < x.n_elem; ++i) {
    if (std::find(seen.begin(), seen.end(), x[i]) == seen.end()) {
      seen.push_back(x[i]);
    }
  }
  return arma::conv_to<arma::ivec>::from(seen);
}

// Function to generate dummy variables
// [[Rcpp::export]]
arma::mat ArmaDummyMat(const arma::imat& mat) {
  int n = mat.n_rows;   // Number of rows
  int p = mat.n_cols;   // Number of columns

  // Vector to hold dummy matrices
  arma::mat combined_dummies;

  // Create dummy variables for each column
  for (int col_idx = 0; col_idx < p; ++col_idx) {
    arma::ivec x = mat.col(col_idx);
    arma::ivec levels = ArmaRunique(x);  // Get unique levels

    // n-1 dummy variables for each column
    arma::mat dummies(n, levels.n_elem - 1, fill::zeros);

    for (int level_idx = 0; level_idx < levels.n_elem - 1; ++level_idx) {
      dummies.col(level_idx) = arma::conv_to<arma::vec>::from(x == levels(level_idx));
    }

    // Concatenate horizontally with the previously generated dummy columns
    if (combined_dummies.n_elem == 0) {
      combined_dummies = dummies;  // Initialize with the first set of dummies
    } else {
      combined_dummies = arma::join_horiz(combined_dummies, dummies);  // Append new dummies
    }
  }

  return combined_dummies;  // Return combined dummy matrix
}
