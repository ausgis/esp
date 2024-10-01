#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string InteractionType(double qv12, double qv1, double qv2) {

  // Step 1: Initialize the interaction string
  std::string interaction;

  // Step 2: Apply conditional logic based on the provided R code
  if (qv12 < std::min(qv1, qv2)) {
    interaction = "Weaken, nonlinear";
  } else if (qv12 >= std::min(qv1, qv2) && qv12 <= std::max(qv1, qv2)) {
    interaction = "Weaken, uni-";
  } else if (qv12 > std::max(qv1, qv2) && qv12 < (qv1 + qv2)) {
    interaction = "Enhance, bi-";
  } else if (qv12 == (qv1 + qv2)) {
    interaction = "Independent";
  } else {
    interaction = "Enhance, nonlinear";
  }

  // Step 3: Return the interaction result
  return interaction;
}
