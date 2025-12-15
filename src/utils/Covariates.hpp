#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

struct Covariates {
  
    /** @brief Adjacency matrix defining spatial neighborhood structure */
    Eigen::MatrixXi W;

    /** @brief Coefficient controlling the strength of spatial dependency */
    double spatial_coefficient;

    /** @brief Ages list  */
    Eigen::VectorXi ages;

    Covariates(Eigen::MatrixXi W = Eigen::MatrixXi(), 
               double spatial_coefficient = 1,
               Eigen::VectorXi ages = Eigen::VectorXi()) 
               : W(W), spatial_coefficient(spatial_coefficient), ages(ages) {}

};

// Expose Covariates to R
RCPP_EXPOSED_CLASS(Covariates);