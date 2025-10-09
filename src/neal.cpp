/**
 * @file neal.cpp
 * @brief Implementation of Neal's Algorithm 3 for collapsed Gibbs sampling
 *
 * This file contains the complete implementation of the Neal3 class, providing
 * efficient collapsed Gibbs sampling for Bayesian nonparametric mixture models.
 * The implementation includes sequential observation updates with automatic
 * cluster management and proper probability computations.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "neal.hpp"
#include "Process.hpp"
#include <Rcpp.h>
#include <algorithm>
#include <random>
using namespace Rcpp;

void Neal3::step_1_observation(int index) {
  /**
   * @brief Performs a step in the DPNeal2 sampling process.
   * @details This method is responsible for updating the allocations of the
   * data points based on the current state of the model.
   * @param index The index of the data point to update.
   */

  // Set unallocated the index
  data.set_allocation(index, -1);

  // for each cluster, compute the log likelihood of the point being in that
  // cluster
  std::vector<double> log_likelihoods(data.get_K() + 1, 0.0);
  for (int k = 0; k < data.get_K(); ++k)
    log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k);

  // Compute the log likelihood of the point being in a new cluster
  log_likelihoods[data.get_K()] =
      likelihood.point_loglikelihood_cond(index, data.get_K());

  // multiply by the prior probability of the cluster
  for (int k = 0; k < data.get_K(); ++k) {
    log_likelihoods[k] += process.gibbs_prior_existing_cluster(k, index);
  }
  log_likelihoods[data.get_K()] += process.gibbs_prior_new_cluster();

  // Normalize the log likelihoods
  double max_loglik =
      *std::max_element(log_likelihoods.begin(), log_likelihoods.end());
  std::vector<double> probs(log_likelihoods);

  for (double &prob : probs) {
    prob = exp(prob - max_loglik);
  }
  double sum_probs = std::accumulate(probs.begin(), probs.end(), 0.0);
  for (double &prob : probs) {
    prob /= sum_probs;
  }

  // DEBUG
  // Rcpp::Rcout << "Log likelihoods and probabilities for data point " << index
  // << ":" << std::endl; for (int i = 0; i < probs.size(); ++i) {
  //     Rcpp::Rcout << "\tCluster " << i << ": l = " << log_likelihoods[i] << "
  //     p = "<< probs[i] << std::endl;
  // }

  // Sample a cluster based on the probabilities
  std::discrete_distribution<int> dist(probs.begin(), probs.end());
  int sampled_cluster = dist(gen);

  // Set the allocation for the data point
  data.set_allocation(index, sampled_cluster);
}

void Neal3::step() {
  /**
   * @brief Performs a single step of the DP Neal 2 algorithm for all the
   * dataset.
   */

  for (int j = 0; j < data.get_n(); ++j) {
    step_1_observation(j);
  }

  // Rcpp::Rcout << std::endl << "Element per cluster: " << std::endl;
  // for (int k = 0; k < data.get_K(); ++k) {
  //     Rcpp::Rcout << "\tCluster " << k << ": " << data.get_cluster_size(k) <<
  //     std::endl;
  // }
  // Rcpp::Rcout <<
  // "-----------------------------------------------------------------------------"
  // << std::endl;
}