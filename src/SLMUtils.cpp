#include <Rcpp.h>
#include <iomanip> // for setting decimal precision and alignment
#include <algorithm> // for max_element
#include <sstream> // for stringstream
using namespace Rcpp;

// Helper function to calculate the maximum width of a column (for character vectors)
int maxColumnWidth(CharacterVector col) {
  int maxLen = 0;
  for (int i = 0; i < col.size(); ++i) {
    maxLen = std::max(maxLen, (int)as<std::string>(col[i]).size());
  }
  return maxLen;
}

// Helper function to calculate the maximum width of a column (for numeric vectors)
int maxColumnWidth(NumericVector col, int precision) {
  int maxLen = 0;
  for (int i = 0; i < col.size(); ++i) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << col[i];
    maxLen = std::max(maxLen, (int)oss.str().size());
  }
  return maxLen;
}

// [[Rcpp::export]]
void PrintGloalQ(DataFrame df) {
  int n = df.nrows();
  CharacterVector colNames = df.names();
  std::vector<int> colWidths(df.size());
  int precision = 3; // Set precision for numeric columns
  int spacing = 1; // Space between columns

  // Calculate the width of each column
  for (int j = 0; j < df.size(); ++j) {
    if (Rcpp::is<CharacterVector>(df[j])) {
      CharacterVector col = df[j];
      colWidths[j] = maxColumnWidth(col);
    } else if (Rcpp::is<NumericVector>(df[j])) {
      NumericVector col = df[j];
      colWidths[j] = maxColumnWidth(col, precision);
    }
    // Ensure minimum width based on the column name length
    colWidths[j] = std::max(colWidths[j], (int)as<std::string>(colNames[j]).size());
  }

  // Print column names with appropriate spacing and left alignment
  for (int j = 0; j < df.size(); ++j) {
    Rcpp::Rcout << std::left << std::setw(colWidths[j] + spacing) << as<std::string>(colNames[j]);
  }
  Rcpp::Rcout << std::endl;

  // Print data rows with left alignment and adaptive column width
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < df.size(); ++j) {
      if (Rcpp::is<CharacterVector>(df[j])) {
        // Non-numeric column
        CharacterVector col = df[j];
        Rcpp::Rcout << std::left << std::setw(colWidths[j] + spacing) << as<std::string>(col[i]);
      } else if (Rcpp::is<NumericVector>(df[j])) {
        // Numeric column
        NumericVector col = df[j];
        Rcpp::Rcout << std::left << std::setw(colWidths[j] + spacing) << std::fixed << std::setprecision(precision) << col[i];
      }
    }
    Rcpp::Rcout << std::endl;
  }
}

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
