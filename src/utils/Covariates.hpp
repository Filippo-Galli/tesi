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

    double B; // prior variance of the means
    double m; // prior mean of the means
    double v; // prior variance of the covariates effect

    Covariates(Eigen::MatrixXi W = Eigen::MatrixXi(), 
               double spatial_coefficient = 1,
               Eigen::VectorXi ages = Eigen::VectorXi(), 
               double B = 1.0, 
               double m = 0.0, 
               double v = 1.0) 
               : W(W), spatial_coefficient(spatial_coefficient), ages(ages), 
               B(B), m(m), v(v) {}

};

// Expose Covariates to R
RCPP_EXPOSED_CLASS(Covariates);