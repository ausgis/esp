#include <Rcpp.h>
#include <iomanip>
using namespace Rcpp;

// [[Rcpp::export]]
std::string SLMUsed(std::string model, bool durbin) {
  if (model == "ols") {
    return "Ordinary Least Square";
  } else if (model == "lag") {
    if (durbin) {
      return "Spatial Durbin Model";
    } else {
      return "Spatial Lag Model";
    }
  } else if (model == "error") {
    if (durbin) {
      return "Spatial Durbin Error Model";
    } else {
      return "Spatial Error Model";
    }
  } else if (model == "gwr") {
    return "Geographically Weighted Regression";
  } else {
    return "Unknown Model";
  }
}

// [[Rcpp::export]]
void PrintGloalQ(DataFrame df) {
  int n = df.nrows();
  CharacterVector colNames = df.names();

  // Print column names
  for (int j = 0; j < df.size(); ++j) {
    Rcpp::Rcout << std::setw(15) << as<std::string>(colNames[j]) << " ";
  }
  Rcpp::Rcout << std::endl;

  // Print data
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < df.size(); ++j) {
      if (Rcpp::is<CharacterVector>(df[j])) {
        // Non-numeric column
        CharacterVector col = df[j];
        Rcpp::Rcout << std::setw(15) << as<std::string>(col[i]) << " ";
      } else if (Rcpp::is<NumericVector>(df[j])) {
        // Numeric column
        NumericVector col = df[j];
        Rcpp::Rcout << std::setw(15) << std::fixed << std::setprecision(3) << col[i] << " ";
      }
    }
    Rcpp::Rcout << std::endl;
  }
}
