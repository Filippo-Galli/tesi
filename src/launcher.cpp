/**
 * @file launcher.cpp
 * @brief Main launcher file for MCMC simulation with various Bayesian
 * non-parametric processes
 * @author Filippo Galli
 * @date 2025
 *
 * This file contains the main MCMC implementation and R interface for running
 * Bayesian non-parametric clustering using different processes (DP, DPW, NGGP,
 * NGGPW) with various sampling strategies (Neal3, Split-Merge, Split-Merge
 * SAMS).
 */

// [[Rcpp::depends(RcppEigen)]]

#include "samplers/neal.hpp"
#include "samplers/neal_ZDNAM.hpp"
#include "samplers/splitmerge.hpp"
#include "samplers/splitmerge_SAMS.hpp"
// #include "splitmerge_SDDS.hpp"

#include "processes/DP.hpp"
#include "processes/DPW.hpp"
#include "processes/NGGP.hpp"
#include "processes/NGGPW.hpp"

#include "samplers/U_sampler/RWMH.hpp"
#include "samplers/U_sampler/MALA.hpp"

#include "utils/Data.hpp"
#include "utils/Likelihood.hpp"
#include "utils/Params.hpp"

#include <chrono>

/**
 * @brief Rcpp module to expose the Params class to R
 *
 * This module creates an R interface for the Params class, allowing R users
 * to create and manipulate parameter objects for the MCMC simulation.
 * The module exposes constructors and all parameter fields with their
 * descriptions.
 */
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
      .field("coefficient", &Params::coefficient,
             "Coefficient for the spatial dependency")
      .field("W", &Params::W, "Adjacency matrix for the points");
}

/**
 * @brief Main MCMC function for Bayesian non-parametric clustering
 *
 * This function performs Markov Chain Monte Carlo (MCMC) sampling for
 * clustering analysis using various Bayesian non-parametric processes. It
 * supports different process types (DP, DPW, NGGP, NGGPW) and sampling
 * algorithms (Neal3, Split-Merge).
 *
 * @param distances A matrix of distances between data points (Eigen::MatrixXd)
 * @param param Reference to a Params object containing all MCMC parameters and
 * hyperparameters
 * @param initial_allocations_r Optional initial cluster allocations
 * (Rcpp::IntegerVector). If empty, random initialization will be used.
 *
 * @return Rcpp::List containing:
 *         - allocations: List of cluster allocations for each iteration
 *         - K: List of number of clusters for each iteration
 *         - loglikelihood: Vector of log-likelihood values for each iteration
 *         - U: Vector of process-specific parameter U for each iteration
 *
 * @details The function:
 *          1. Initializes data structures and likelihood computation
 *          2. Sets up the chosen Bayesian non-parametric process (currently
 * NGGPW)
 *          3. Configures the sampling algorithm (currently Split-Merge)
 *          4. Runs MCMC for burn-in + sampling iterations
 *          5. Updates process parameters and saves results at each iteration
 *          6. Provides progress output and timing information
 *
 * @note The process type and sampler can be changed by uncommenting/commenting
 *       the appropriate lines in the function body.
 */
// [[Rcpp::export]]
Rcpp::List
mcmc(const Eigen::MatrixXd &distances, Params &param,
     const Rcpp::IntegerVector &initial_allocations_r = Rcpp::IntegerVector()) {
  // Convert R IntegerVector to Eigen::VectorXi for compatibility with
  // Eigen-based classes
  Eigen::VectorXi initial_allocations;
  if (initial_allocations_r.size() > 0) {
    initial_allocations = Rcpp::as<Eigen::VectorXi>(initial_allocations_r);
  }

  // Initialize data structure with distances and optional initial allocations
  Data data(distances, initial_allocations);

  // Initialize likelihood computation object
  Likelihood likelihood(data, param);

  // Sampler for U parameter if needed by the process
  RWMH U_sampler(param, data, true, 2.0, true);
  //MALA U_sampler(param, data, true, 2, true);

  // Initialize the Bayesian non-parametric process
  // Uncomment the desired process type:
  //DP process(data, param);      // Dirichlet Process
  //DPW process(data, param);     // Dirichlet Process with Weights
  //NGGP process(data, param, U_sampler);    // Normalized Generalized Gamma Process
  NGGPW process(data,param, U_sampler); // Normalized Generalized Gamma Process with Weights

  // Initialize sampling algorithms
  Neal3 gibbs(data, param, likelihood,process); // Gibbs sampler (Neal Algorithm 3)
  //Neal3ZDNAM gibbs(data, param, likelihood, process); // Gibbs sampler with ZDNAM

  // Choose the main sampling strategy:
  SplitMerge sampler(data,param, likelihood, process, true); // Split-Merge sampler
  //SplitMerge_SAMS sampler(data, param, likelihood, process, true);    // Split-Merge with SAMS 
  // SplitMerge_SDDS sampler(data, param, likelihood, process, true);    // Split-Merge with SDDS

  // Initialize results container to store MCMC output
  Rcpp::List results = Rcpp::List::create(
      Rcpp::Named("allocations") = Rcpp::List(param.NI + param.BI),                // Cluster assignments
      Rcpp::Named("K") = Rcpp::List(param.NI + param.BI), // Number of clusters
      Rcpp::Named("loglikelihood") = Rcpp::NumericVector(param.NI + param.BI, 0.0), // Log-likelihood values
      Rcpp::Named("U") = Rcpp::NumericVector(param.NI + param.BI, 0.0) // Process parameter U
  );

  // Print MCMC configuration and start timing
  std::cout << "Starting MCMC with " << param.NI << " iterations after " << param.BI << " burn-in iterations." << std::endl;

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Calculate progress reporting interval (5% of total iterations)
  int printing_interval = int((param.NI + param.BI) * 0.05);
  //int printing_interval = 1;
  double iter_s = 0.0; // Iterations per second

  // Main MCMC loop
  for (int i = 0; i < param.NI + param.BI; ++i) {

    // Update process-specific parameters
    process.update_params();

    // Perform one MCMC step using the chosen sampler
    sampler.step();

    // Optional: Perform Gibbs step 
    if(i % 25 == 0)
      gibbs.step();

    // Store results for current iteration
    Rcpp::as<Rcpp::List>(results["allocations"])[i] = Rcpp::wrap(data.get_allocations());
    Rcpp::as<Rcpp::List>(results["K"])[i] = data.get_K();
    Rcpp::as<Rcpp::NumericVector>(results["U"])[i] = U_sampler.get_U();

    // Print progress information at regular intervals
    if ((i + 1) % printing_interval == 0) {
      std::cout << "Iteration " << i + 1 << ": ";
      auto elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count() / 1000.0;
      if (elapsed_seconds > 0) {
        iter_s = static_cast<double>(i + 1) / elapsed_seconds;
      } else {
        iter_s = 0.0;
      }
      std::cout << "Clusters: " << data.get_K() << " - iter/s: " << iter_s 
      << " | time to complete: " << (param.BI + param.NI - i)/iter_s<< " s" << std::endl;
    }
  }

  // Calculate and display final timing and acceptance statistics
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "MCMC completed in : " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " seconds." << std::endl;
  std::cout << "U acceptance rate: " << U_sampler.get_acceptance_rate() * 100 << " %." << std::endl;

  // std::cout << "Accepted split ratio: " << sampler.get_accepted_split() * 100 * 2 / (param.NI + param.BI) << " %." << std::endl;
  // std::cout << "Accepted merge ratio: " << sampler.get_accepted_merge() * 100 * 2 / (param.NI + param.BI) << " %." << std::endl;  
  // std::cout << "Accepted shuffle ratio: " << sampler.get_accepted_shuffle() * 100 / (param.NI + param.BI) << " %." << std::endl;

  return results;
}
