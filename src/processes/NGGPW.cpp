/**
 * @file NGGPW.cpp
 * @brief Implementation of Weighted Normalized Generalized Gamma Process
 *
 * This file contains the implementation of the NGGPW class, which combines
 * the flexibility of NGGP with spatial dependency modeling. This provides
 * the most sophisticated clustering approach with both enhanced cluster
 * control and spatial awareness.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "NGGPW.hpp"

Eigen::VectorXd NGGPW::gibbs_prior_existing_clusters(int obs_idx) const {
  /**
   * @brief Computes the log prior probabilities of assigning a data point to
   * all existing clusters.
   *
   * This method incorporates spatial information by considering the number of
   * neighbors in each target cluster when computing the prior probabilities.
   * @param obs_idx The index of the observation to assign.
   * @return A vector of log prior probabilities for assigning the data point to
   * each existing cluster.
   */

  Eigen::VectorXd priors = Eigen::VectorXd::Zero(data.get_K());
  Eigen::VectorXi neighbors = get_neighbors_obs(obs_idx);

  // Compute prior for each existing cluster
  for (int k = 0; k < data.get_K(); ++k) {
    const int cluster_size = data.get_cluster_size(k);
    double prior = cluster_size  - params.sigma > 0 ? log(cluster_size - params.sigma)
                                    : std::numeric_limits<double>::lowest();
    priors(k) = prior + params.coefficient * neighbors(k);
  }

  return priors;
}

double NGGPW::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {
  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   *
   * For NGGPW, this combines the NGGP discount parameter (n_k - sigma) with
   * spatial information from the adjacency matrix W.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation to assign.
   * @return The log prior probability of assigning the data point to the
   * existing cluster.
   */

  const int cluster_size = data.get_cluster_size(cls_idx);
  double prior = params.coefficient * get_neighbors_obs(obs_idx, cls_idx);
  prior = cluster_size - params.sigma > 0 ? prior + log(cluster_size - params.sigma) : std::numeric_limits<double>::lowest();
  return prior;
}

double NGGPW::gibbs_prior_new_cluster() const {
  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   *
   * For NGGPW, this follows the NGGP formulation and is proportional to
   * alpha * sigma * (tau + U)^sigma. New clusters have zero spatial neighbors.
   * @return The log prior probability of assigning the data point to a new
   * cluster.
   */

  return log_a + params.sigma * log(params.tau + U_sampler_method.get_U());
}

double NGGPW::prior_ratio_split(int ci, int cj) const {
  /**
   * @brief Computes the prior ratio for a split operation in a spatially-aware
   * NGGP-based split-merge MCMC algorithm.
   *
   * This method accounts for both the generalized gamma process prior and
   * spatial dependencies when computing the acceptance ratio for splitting
   * clusters.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */

  const int n_ci = data.get_cluster_size(ci);
  const int n_cj = data.get_cluster_size(cj);

  // NGGP prior part
  double log_acceptance_ratio = 0.0;
  log_acceptance_ratio += log_a;
  log_acceptance_ratio += params.sigma * log(params.tau + U_sampler_method.get_U());
  log_acceptance_ratio -= n_ci + n_cj - params.sigma > 0 ? lgamma(n_ci + n_cj - params.sigma) : 0;
  log_acceptance_ratio += n_ci - params.sigma > 0 ? lgamma(n_ci - params.sigma) : 0;
  log_acceptance_ratio += n_cj - params.sigma > 0 ? lgamma(n_cj - params.sigma) : 0;

  // Spatial part: add new clusters, subtract old merged cluster
  log_acceptance_ratio += params.coefficient * get_neighbors_cls(ci);
  log_acceptance_ratio += params.coefficient * get_neighbors_cls(cj);
  log_acceptance_ratio -= params.coefficient * get_neighbors_cls(ci, true);

  return log_acceptance_ratio;
}

double NGGPW::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
  /**
   * @brief Computes the prior ratio for a merge operation in a spatially-aware
   * NGGP-based split-merge MCMC algorithm.
   *
   * This method accounts for both the generalized gamma process prior and
   * spatial dependencies when computing the acceptance ratio for merging
   * clusters.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */

  const int size_merge = size_old_ci + size_old_cj;

  // NGGP prior part
  double log_acceptance_ratio = -log_a;
  log_acceptance_ratio -= params.sigma * log(params.tau + U_sampler_method.get_U());
  log_acceptance_ratio += size_merge - params.sigma > 0 ? lgamma(size_merge - params.sigma) : 0;
  log_acceptance_ratio -= size_old_ci - params.sigma > 0 ? lgamma(size_old_ci - params.sigma) : 0;
  log_acceptance_ratio -= size_old_cj - params.sigma > 0 ? lgamma(size_old_cj - params.sigma) : 0;

  // Spatial part
  const int old_ci = old_allocations[idx_i];
  const int old_cj = old_allocations[idx_j];
  log_acceptance_ratio -= params.coefficient * get_neighbors_cls(old_ci, true);
  log_acceptance_ratio -= params.coefficient * get_neighbors_cls(old_cj, true);

  const int new_ci = data.get_allocations()[idx_i];
  log_acceptance_ratio += params.coefficient * get_neighbors_cls(new_ci);

  return log_acceptance_ratio;
}

double NGGPW::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci,
                                  int cj) const {
  /**
   * @brief Computes the prior ratio for a shuffle operation in a
   * spatially-aware NGGP-based split-merge MCMC algorithm.
   *
   * This method accounts for both the generalized gamma process prior and
   * spatial dependencies when computing the acceptance ratio for shuffling
   * observations between clusters.
   * @param size_old_ci The size of the first cluster before the shuffle.
   * @param size_old_cj The size of the second cluster before the shuffle.
   * @param ci The first cluster index involved in the shuffle.
   * @param cj The second cluster index involved in the shuffle.
   * @return The log prior ratio for the shuffle operation.
   */

  const int n_ci = data.get_cluster_size(ci);
  const int n_cj = data.get_cluster_size(cj);

  double log_acceptance_ratio = 0.0;
  log_acceptance_ratio += n_ci - params.sigma > 0 ? lgamma(n_ci - params.sigma) : 0;
  log_acceptance_ratio += n_cj - params.sigma > 0 ? lgamma(n_cj - params.sigma) : 0;
  log_acceptance_ratio -= size_old_ci - params.sigma > 0 ? lgamma(size_old_ci - params.sigma) : 0;
  log_acceptance_ratio -= size_old_cj - params.sigma > 0 ? lgamma(size_old_cj - params.sigma) : 0;

  // Spatial part
  log_acceptance_ratio += params.coefficient * get_neighbors_cls(ci);
  log_acceptance_ratio += params.coefficient * get_neighbors_cls(cj);
  log_acceptance_ratio -= params.coefficient * get_neighbors_cls(ci, true);
  log_acceptance_ratio -= params.coefficient * get_neighbors_cls(cj, true);

  return log_acceptance_ratio;
}