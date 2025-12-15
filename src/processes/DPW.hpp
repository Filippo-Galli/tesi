/**
 * @file DPW.hpp
 * @brief Dirichlet Process with spatial weights (DPW) implementation for
 * spatially-aware Bayesian nonparametric clustering.
 */

#pragma once

#include "../utils/Process.hpp"
#include "module/spatial_module.hpp"
#include "Eigen/src/Core/Matrix.h"

/**
 * @class DPW
 * @brief Dirichlet Process with spatial Weights class for spatially-aware
 * Bayesian nonparametric clustering.
 *
 * This class extends the Dirichlet Process to incorporate spatial information
 * through an adjacency matrix W. It provides methods for Gibbs sampling and
 * split-merge algorithms that account for spatial dependencies between
 * observations in the clustering process.
 */
class DPW : public Process, protected SpatialModule {

public:
  /**
   * @brief Constructor for the Dirichlet Process with spatial Weights.
   * @param d Reference to the data object containing observations and cluster
   * assignments.
   * @param p Reference to the parameters object containing the adjacency matrix
   * W and spatial coefficient.
   */
  DPW(Data &d, Params &p, Covariates &c) 
    : Process(d, p), 
      SpatialModule(c, d, [this]() -> const Eigen::VectorXi& { return this->old_allocations; }) {};

  /**
   * @name Gibbs Sampling Methods
   * @{
   */

  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   *
   * This method incorporates spatial information by considering the number of
   * neighbors in the target cluster when computing the prior probability.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation to assign.
   * @return The log prior probability of assigning the data point to the
   * existing cluster.
   */
  [[nodiscard]] double gibbs_prior_existing_cluster(int cls_idx,
                                                    int obs_idx) const override;

  /**
    * @brief Computes the log prior probabilities of assigning a data point to every existing cluster.
    * This method is useful for Gibbs sampling over existing clusters.
    * It returns a vector of log prior probabilities for all existing clusters.
    * @param obs_idx The index of the observation to assign.
    * @return A vector of log prior probabilities for assigning the data point to each existing cluster.
  */                                     
  [[nodiscard]] Eigen::VectorXd gibbs_prior_existing_clusters(int obs_idx) const override;

  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   * @return The log prior probability of assigning the data point to a new
   * cluster.
   */
  [[nodiscard]] double gibbs_prior_new_cluster() const override;

  /** @} */

  /**
   * @name Split-Merge Algorithm Methods
   * @{
   */

  /**
   * @brief Computes the prior ratio for a split operation in a spatially-aware
   * split-merge MCMC algorithm.
   *
   * This method accounts for both the Dirichlet Process prior and spatial
   * dependencies when computing the acceptance ratio for splitting clusters.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */
  [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;

  /**
   * @brief Computes the prior ratio for a merge operation in a spatially-aware
   * split-merge MCMC algorithm.
   *
   * This method accounts for both the Dirichlet Process prior and spatial
   * dependencies when computing the acceptance ratio for merging clusters.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */
  [[nodiscard]] double prior_ratio_merge(int size_old_ci,
                                         int size_old_cj) const override;

  /**
   * @brief Computes the prior ratio for a shuffle operation in a
   * spatially-aware split-merge MCMC algorithm.
   *
   * This method accounts for both the Dirichlet Process prior and spatial
   * dependencies when computing the acceptance ratio for shuffling observations
   * between clusters.
   * @param size_old_ci The size of the first cluster before the shuffle.
   * @param size_old_cj The size of the second cluster before the shuffle.
   * @param ci The first cluster index involved in the shuffle.
   * @param cj The second cluster index involved in the shuffle.
   * @return The log prior ratio for the shuffle operation.
   */
  [[nodiscard]] double prior_ratio_shuffle(int size_old_ci, int size_old_cj,
                                           int ci, int cj) const override;

  /** @} */

  /**
   * @note Spatial methods (get_neighbors_obs, get_neighbors_cls) are inherited
   * from SpatialModule with optimized caching implementation.
   */

  /**
   * @brief Updates the parameters of the Dirichlet Process with spatial
   * Weights.
   *
   * This is a null implementation as the DPW has no parameters to update beyond
   * those handled by the base Process class.
   */
  void update_params() override { return; };
};