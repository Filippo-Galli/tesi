#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

struct Params {
    // Distribution parameters
    double delta1; // parameter for the first gamma
    double alpha; // parameter for the lambda_k gamma
    double beta; // parameter for the lambda_k gamma

    double delta2; // parameter for the second gamma
    double gamma; // parameter for the theta_kt gamma
    double zeta; // parameter for the theta_kt gamma

    int BI; // number of burn-in iterations
    int NI; // number of iterations after burn-in
    
    double a; // total mass
    double sigma; // second parameter of the NGGP
    double tau; // third parameter of the NGGP

    double coefficient; // coefficient for the spatial dependency
    Eigen::MatrixXi W; // adjacency matrix for the points


    Params(double delta1 = 0.5, double alpha = 2, double beta = 2,
           double delta2 = 2, double gamma = 2, double zeta = 2,
           int BI = 1000, int NI = 10000, double a = 1.0, double sigma = 1.0, double tau = 1.0, double coefficient = 1, Eigen::MatrixXi W = Eigen::MatrixXi())
        : delta1(delta1), alpha(alpha), beta(beta),
          delta2(delta2), gamma(gamma), zeta(zeta),
          BI(BI), NI(NI), a(a), sigma(sigma), tau(tau), coefficient(coefficient), W(W) {}
};

// Expose Params to R
RCPP_EXPOSED_CLASS(Params);