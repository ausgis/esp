#include <Rcpp.h>
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
