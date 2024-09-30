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

// calculating equivalent q-values under the framework of linear regression

// [[Rcpp::export]]
Rcpp::NumericVector CalculateQ(const arma::mat& y_pred,
                               const arma::imat& discmat,
                               const arma::vec& y) {
  int p = y_pred.n_cols;  // Number of columns (matches discmat)
  // int n = y_pred.n_rows;  // Number of rows (matches discmat and y)

  arma::vec q(p, fill::zeros);  // Initialize result vector

  // Iterate over each column of y_pred and discmat
  for (int col_idx = 0; col_idx < p; ++col_idx) {
    arma::ivec disc_col = discmat.col(col_idx); // Current column of discmat
    arma::ivec levels = ArmaRunique(disc_col);  // Get unique levels in this column

    arma::vec total_diff = square(y_pred.col(col_idx) - y);  // Squared differences for y_pred and y
    double total_sum = sum(total_diff);  // Total sum of squared differences

    arma::vec group_sums(levels.n_elem, fill::zeros);  // Store group-wise sum of squared differences

    // Compute the sum of squared differences for each level
    int levelidx_length = levels.n_elem;
    for (int level_idx = 0; level_idx < levelidx_length; ++level_idx) {
      int current_level = levels(level_idx);

      // Mask for the current level
      arma::uvec mask = find(disc_col == current_level);

      // Sum the squared differences for this level
      group_sums(level_idx) = sum(total_diff.elem(mask));
    }

    // Compute q for this column as 1 - (group_sum / total_sum)
    double sum_of_ratios = sum(group_sums / total_sum);
    q(col_idx) = 1.0 - sum_of_ratios;
  }

  // Convert q to Rcpp::NumericVector and return
  return Rcpp::wrap(q);
}

// [[Rcpp::export]]
Rcpp::NumericVector SLMQ(const arma::imat& levelmat,
                         const arma::vec& coefs,
                         const arma::vec& y){
  arma::mat y_pred = PredictDummyY(levelmat,coefs);
  Rcpp::NumericVector qv = CalculateQ(y_pred,levelmat,y);
  return qv;
}
