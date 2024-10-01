#include <RcppArmadillo.h>
#include "Dummy.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double ComputeR2(const arma::vec& y, const arma::vec& y_pred) {
  // Calculate the mean of y
  double mean_y = arma::mean(y);

  // Calculate total variance (total sum of squares)
  double total_variance = arma::sum(arma::pow(y - mean_y, 2));

  // Calculate residual sum of squares (RSS)
  double residual_sum_of_squares = arma::sum(arma::pow(y - y_pred, 2));

  // Calculate R-squared
  double r_squared = 1 - (residual_sum_of_squares / total_variance);

  return r_squared;
}

// [[Rcpp::export]]
Rcpp::NumericVector SLMQ(const arma::mat& FitY, const arma::vec& Y) {
  // Get the number of columns in FitY
  int n_cols = FitY.n_cols;

  // Initialize a NumericVector to store R-squared values
  Rcpp::NumericVector r_squared_values(n_cols);

  // Loop through each column of FitY
  for (int i = 0; i < n_cols; ++i) {
    // Extract the i-th column of FitY
    arma::vec y_pred = FitY.col(i);

    // Compute R-squared for the current column
    double r2 = ComputeR2(Y, y_pred);

    // Store the R-squared value in the result vector
    r_squared_values[i] = r2;
  }

  // Return the vector of R-squared values
  return r_squared_values;
}

// Function to compute R-squared for each group defined by Zones
// [[Rcpp::export]]
Rcpp::List SLMLocalQ(const arma::mat& FitY,
                     const arma::vec& Y,
                     const arma::ivec& Zones) {
  // Get unique levels in Zones using ArmaRunique
  arma::ivec unique_zones = ArmaRunique(Zones);
  int n_zones = unique_zones.n_elem; // Number of unique zones
  int n_cols = FitY.n_cols;          // Number of columns in FitY

  // Create a list to store the R-squared vectors for each zone
  Rcpp::List results(n_zones);

  // Loop through each unique zone
  for (int i = 0; i < n_zones; ++i) {
    int zone = unique_zones[i];  // Current zone value

    // Get the indices of the rows in the current zone
    arma::uvec indices = arma::find(Zones == zone);

    // Subset Y and FitY by the current zone
    arma::vec Y_zone = Y.elem(indices);          // Subset of Y in current zone
    arma::mat FitY_zone = FitY.rows(indices);    // Subset of FitY in current zone

    // Initialize a NumericVector to store R-squared for each column in this zone
    Rcpp::NumericVector r_squared_values(n_cols);

    // Compute R-squared for each column in the current zone
    for (int j = 0; j < n_cols; ++j) {
      arma::vec y_pred = FitY_zone.col(j);       // Predicted values in current zone
      double r2 = ComputeR2(Y_zone, y_pred);     // Compute R-squared for the current column
      r_squared_values[j] = r2;                  // Store the R-squared value
    }

    // Store the R-squared vector for this zone
    results[i] = r_squared_values;
  }

  // Return the list of R-squared values for each zone
  return results;
}
