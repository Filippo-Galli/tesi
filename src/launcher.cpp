// [[Rcpp::depends(RcppEigen)]]

#include "neal.hpp"

#include "DP.hpp"

#include "Data.hpp"
#include "Likelihood.hpp"
#include "Params.hpp"
#include "Process.hpp"
#include "Rcpp/vector/instantiation.h"
#include <chrono>

// Expose Params class to R using Rcpp modules
RCPP_MODULE(params_module) {
  Rcpp::class_<Params>("Params")
      .constructor<double, double, double, double, double, double, int, int,
                   double, double, double, double, Eigen::MatrixXi>("Constructor with W matrix")
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
      .field("coefficient", &Params::coefficient, "Coefficient for the spatial dependency")
      .field("W", &Params::W, "Adjacency matrix for the points");
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
  
  DP process(data, param);

  Neal3 sampler(data, param, likelihood, process);

  Rcpp::List results = Rcpp::List::create(
      Rcpp::Named("allocations") = Rcpp::List(param.NI + param.BI),
      Rcpp::Named("K") = Rcpp::List(param.NI + param.BI),
      Rcpp::Named("loglikelihood") = Rcpp::NumericVector(param.NI + param.BI, 0.0),
      Rcpp::Named("U") = Rcpp::NumericVector(param.NI + param.BI, 0.0)
    );


  std::cout << "Starting MCMC with " << param.NI << " iterations after "
            << param.BI << " burn-in iterations." << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  int printing_interval = int((param.NI + param.BI)*0.05);
  double iter_s = 0.0;

  for (int i = 0; i < param.NI + param.BI; ++i) {
    
    sampler.step();

    // if(i % 20 == 0)
    //   gibbs.step();

    // Save intermediate results
    Rcpp::as<Rcpp::List>(results["allocations"])[i] =
        Rcpp::wrap(data.get_allocations());
    Rcpp::as<Rcpp::List>(results["K"])[i] = data.get_K();
    //Rcpp::as<Rcpp::NumericVector>(results["U"])[i] = sampler.get_U();

    // print intermediate results
    if ((i + 1) % printing_interval == 0) {
      std::cout << "Iteration " << i + 1 << ": ";
      auto elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count() / 1000.0;
      if (elapsed_seconds > 0) {
        iter_s = static_cast<double>(i + 1) / elapsed_seconds;
      } else {
        iter_s = 0.0;
      }
      std::cout << "Number of clusters: " << data.get_K() << " iter/s: " << iter_s << std::endl;
    }
  }

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "MCMC completed in : " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " seconds." << std::endl;

  return results;
}