// [[Rcpp::depends(RcppEigen)]]

#include "DP_neal2.hpp"
#include "Data.hpp"
#include "Likelihood.hpp"
#include "Params.hpp"
#include "Rcpp/vector/instantiation.h"

// Expose Params class to R using Rcpp modules
RCPP_MODULE(params_module) {
  Rcpp::class_<Params>("Params")
      .constructor<double, double, double, double, double, double, int, int,
                   double, double, double>("Constructor with all parameters")
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
      .field("tau", &Params::tau, "Third parameter of the NGGP");
}

// [[Rcpp::export]]
Rcpp::List
mcmc(const Eigen::MatrixXd &distances, Params &param,
     const Rcpp::IntegerVector &initial_allocations_r = Rcpp::IntegerVector()) {
  // Convert R IntegerVector to Eigen::VectorXi
  Eigen::VectorXi initial_allocations;
  if (initial_allocations_r.size() > 0) {
    initial_allocations = Rcpp::as<Eigen::VectorXi>(initial_allocations_r);
  }

  Data data(distances, initial_allocations);
  Likelihood likelihood(data, param);
  DPNeal2 sampler(data, param, likelihood);

  Rcpp::List results = Rcpp::List::create(
      Rcpp::Named("allocations") = Rcpp::List(param.NI),
      Rcpp::Named("K") = Rcpp::List(param.NI),
      Rcpp::Named("loglikelihood") = Rcpp::NumericVector(param.NI, 0.0));

  for (int i = 0; i < param.NI + param.BI; ++i) {
    for (int j = 0; j < data.get_n(); ++j) {
      sampler.step(j);
    }

    // Save intermediate results
    if (i >= param.BI) {
      Rcpp::as<Rcpp::List>(results["allocations"])[i - param.BI] =
          Rcpp::wrap(data.get_allocations());
      Rcpp::as<Rcpp::List>(results["K"])[i - param.BI] = data.get_K();
      for (int k = 0; k < data.get_K(); ++k) {
        Rcpp::as<Rcpp::NumericVector>(results["loglikelihood"])[i - param.BI] +=
            likelihood.cluster_loglikelihood(k);
      }
    }

    // print intermediate results
    if (i % 100 == 0) {
      std::cout << "Iteration " << i << ": ";
      std::cout << "Number of clusters: " << data.get_K() << std::endl;
    }
  }

  std::cout << "MCMC completed." << std::endl;

  return results;
}