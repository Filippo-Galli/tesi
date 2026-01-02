/**
 * @file Covariates.hpp
 * @brief Covariates data structure for clustering processes.
 */

#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

/**
 * @struct Covariates
 * @brief Data structure to hold covariate information for clustering.
 */

struct Covariates {
    
    /**
     * @name Spatial Adjacency Matrix and Coefficient
     * @{
     */
    /** @brief Adjacency matrix defining spatial neighborhood structure */
    Eigen::MatrixXi W;

    /** @brief Coefficient controlling the strength of spatial dependency */
    double spatial_coefficient;

    /** @} */

    /**
     * @name Continuous Covariates
     * @{
     */

    bool fixed_v; ///< Flag indicating if variance v is fixed
    double B;     ///< prior variance of the means or variance multiplier (higher B means diffuse prior)
    double m;     ///< prior mean of the means
    double v;     ///< prior variance of the covariates effect
    double nu;    ///< degrees of freedom for the variance prior
    double S0;    ///< scale parameter for the variance prior

    /** @brief Ages list  */
    Eigen::VectorXd ages;

    /** @} */

    Covariates(Eigen::MatrixXi W = Eigen::MatrixXi(), double spatial_coefficient = 1,
               Eigen::VectorXd ages = Eigen::VectorXd(), double B = 1.0, double m = 0.0, double v = 1.0,
               bool fixed_v = false, double nu = 1.0, double S0 = 1.0)
        : W(W), spatial_coefficient(spatial_coefficient), ages(ages), B(B), m(m), v(v), fixed_v(fixed_v), nu(nu),
          S0(S0) {}
};
