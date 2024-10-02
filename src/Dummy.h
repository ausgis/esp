#ifndef Dummy_H
#define Dummy_H

#include <RcppArmadillo.h>

arma::ivec ArmaRunique(const arma::ivec& x);
arma::mat ArmaDummyVec(const arma::vec& vec);
arma::mat ArmaDummyMat(const arma::imat& mat);

#endif // Dummy_H
