// [[Rcpp::depends(RcppEigen)]]

#include "Data.hpp"

// [[Rcpp::export]]
void mcmc (const Eigen::MatrixXd& distances){
    Data data(distances);
    data.print();
    // Additional MCMC logic would go here
}