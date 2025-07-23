// [[Rcpp::depends(RcppEigen)]]

#include "Data.hpp"
#include "Params.hpp"
#include "Likelihood.hpp"

// Expose Params class to R using Rcpp modules
RCPP_MODULE(params_module) {
    Rcpp::class_<Params>("Params")
        .constructor<double, double, double, double, double, double, int, int, double, double, double>("Constructor with all parameters")
        .constructor("Default constructor")
        .field("delta1", &Params::delta1, "Parameter for the first gamma")
        .field("alpha", &Params::alpha, "Parameter for the lambda_k gamma")
        .field("beta", &Params::beta, "Parameter for the lambda_k gamma")
        .field("delta2", &Params::delta2, "Parameter for the second gamma")
        .field("gamma", &Params::gamma, "Parameter for the theta_kt gamma")
        .field("zeta", &Params::zeta, "Parameter for the theta_kt gamma")
        .field("BI", &Params::BI, "Number of burn-in iterations")
        .field("NI", &Params::NI, "Number of iterations after burn-in")
        .field("a", &Params::a, "Total mass")
        .field("sigma", &Params::sigma, "Second parameter of the NGGP")
        .field("tau", &Params::tau, "Third parameter of the NGGP")
    ;
}

// [[Rcpp::export]]
void mcmc (const Eigen::MatrixXd& distances, const Params& param, const std::string& first_allocation = "all-in-one") {
    Data data(distances, first_allocation);
    Likelihood likelihood(data, param);

    // test likelihood
    double lk = likelihood.cluster_loglikelihood(0);
    std::cout << "Likelihood cluster 0: " << lk << std::endl;

    data.set_allocation(5, 1);
    data.set_allocation(4, 1);

    lk = likelihood.cluster_loglikelihood(0);
    std::cout << "Likelihood cluster 0: " << lk << std::endl;

    data.set_allocation(0, 0);
    data.set_allocation(1, 0);
    data.set_allocation(2, 1);
    data.set_allocation(3, 1);
    data.set_allocation(4, 2);
    data.set_allocation(5, 2);


    lk = likelihood.cluster_loglikelihood(0);
    std::cout << "Likelihood cluster 0: " << lk << std::endl;
    



}