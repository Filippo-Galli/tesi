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

using namespace Rcpp;

// In neal.cpp, could be a private helper in Neal3 or a free function
int Neal3::sample_from_log_probs(const std::vector<double>& log_probs) {
  // 1. Log-Sum-Exp trick for numerical stability
  double max_loglik = *std::max_element(log_probs.begin(), log_probs.end());
  
  std::vector<double> weights(log_probs.size());
  double sum_weights = 0.0;
  for (size_t i = 0; i < log_probs.size(); ++i) {
      weights[i] = exp(log_probs[i] - max_loglik);
      sum_weights += weights[i];
  }

  // 2. Roulette wheel selection
  std::uniform_real_distribution<double> unif(0.0, sum_weights);
  double u = unif(gen);

  int sampled_idx = -1;
  for (size_t i = 0; i < weights.size(); ++i) {
      u -= weights[i];
      if (u < 0.0) {
          sampled_idx = i;
          break;
      }
  }
  // Fallback for floating point inaccuracies
  return (sampled_idx != -1) ? sampled_idx : weights.size() - 1;
}

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

  #pragma omp parallel for
  for (int k = 0; k < data.get_K(); ++k)
    log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k);


  // Compute the log likelihood of the point being in a new cluster
  log_likelihoods[data.get_K()] = likelihood.point_loglikelihood_cond(index, data.get_K());

  // multiply by the prior probability of the cluster
  auto priors = process.gibbs_prior_existing_clusters(index);
  for (int k = 0; k < data.get_K(); ++k) {
    log_likelihoods[k] += priors(k);
  }
  log_likelihoods[data.get_K()] += process.gibbs_prior_new_cluster();

  // Sample a cluster based on the probabilities
  int sampled_cluster = sample_from_log_probs(log_likelihoods);

  // Set the allocation for the data point
  data.set_allocation(index, sampled_cluster);
}

void Neal3::step() {
  /**
   * @brief Performs a single step of the DP Neal 2 algorithm for all the
   * dataset.
   */

  // // Create and shuffle a vector of indices
  // std::vector<int> indices(data.get_n());
  // std::iota(indices.begin(), indices.end(), 0);
  // std::shuffle(indices.begin(), indices.end(), gen);

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