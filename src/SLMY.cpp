#include <RcppArmadillo.h>
#include "Dummy.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// transform each column of the input matrix into dummy variables, and then,
// for each column, retain only the corresponding dummy variables while setting
// all other dummy variables to zero.

// [[Rcpp::export]]
arma::mat PredictDummyY(const arma::imat& mat, const arma::vec& vec){
  int p = mat.n_cols; // Number of columns in the original matrix
  int n = mat.n_rows; // Number of rows

  // Generate the full dummy matrix
  arma::mat dummy_matrix = ArmaDummyMat(mat);

  // Initialize the result matrix
  arma::mat result(n, p, fill::zeros);

  // Number of dummy variables per column
  std::vector<int> levels_count(p);
  int total_dummy_cols = 0;

  // Calculate the number of dummy variables per column
  for (int col_idx = 0; col_idx < p; ++col_idx) {
    arma::ivec levels = ArmaRunique(mat.col(col_idx));
    levels_count[col_idx] = levels.n_elem - 1;
    total_dummy_cols += levels.n_elem - 1;
  }

  // Iterate over each column of the original matrix
  int col_offset = 0;  // Offset to keep track of dummy variable positions
  for (int col_idx = 0; col_idx < p; ++col_idx) {
    // Copy the dummy matrix
    arma::mat modified_dummy = arma::zeros(n, total_dummy_cols);

    // Keep only the dummy variables for the current column, zero out others
    modified_dummy.cols(col_offset, col_offset + levels_count[col_idx] - 1) = dummy_matrix.cols(col_offset, col_offset + levels_count[col_idx] - 1);

    // Perform the multiplication with vec(1:p) and add vec(0) to the result
    result.col(col_idx) = modified_dummy * vec.subvec(1, p) + vec(0);

    // Update the offset for the next column
    col_offset += levels_count[col_idx];
  }

  return result; // Return the final result matrix
}
