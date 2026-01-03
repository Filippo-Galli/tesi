/**
 * @file NGGPW.hpp
 * @brief Normalized Generalized Gamma Process with spatial weights (NGGPW)
 * implementation for spatially-aware Bayesian nonparametric clustering.
 */

#pragma once

#include "../samplers/U_sampler/U_sampler.hpp"
#include "../processes/NGGP.hpp"
#include "module/spatial_module.hpp"

/**
 * @class NGGPW
 * @brief Normalized Generalized Gamma Process with spatial Weights class for
 * spatially-aware Bayesian nonparametric clustering.
 *
 * This class extends the Normalized Generalized Gamma Process (NGGP) to
 * incorporate spatial information through an adjacency matrix W. It combines
 * the flexibility of the NGGP in modeling cluster sizes with spatial
 * dependencies between observations.
 */
class NGGPW : public NGGP, protected SpatialModule {

protected:
  /**
   * @name Private Member Variables
   * @{
   */

  /**
   * @brief Reference to the U_sampler instance for updating the latent variable
   * U.
   *
   * This can be any derived class of U_sampler (e.g., RWMH or MALA) that
   * implements the MCMC algorithm for sampling U from its conditional
   * distribution.
   */
  U_sampler &U_sampler_method;
  /** @} */

  /**
   * @name Random Number Generation
   * @{
   */

  /** @brief Random device for seeding. */
  std::random_device rd;

  /** @brief Mersenne Twister random number generator. */
  mutable std::mt19937 gen;

  /** @} */

public:
  /**
   * @brief Constructor for the Normalized Generalized Gamma Process with
   * spatial Weights.
   * @param d Reference to the data object containing observations and cluster
   * assignments.
   * @param p Reference to the parameters object containing NGGP parameters (a,
   * sigma, tau) and spatial information (W, coefficient).
   * @param mh Reference to a U_sampler instance (e.g., RWMH or MALA) for
   * updating the latent variable U via MCMC.
   */
  NGGPW(Data &d, Params &p, Covariates &c, U_sampler &mh)
      : NGGP(d, p, mh), SpatialModule(c, d,
                                     &this->old_allocations_view()),
        U_sampler_method(mh), gen(rd()) {};

  /**
   * @name Gibbs Sampling Methods
   * @{
   */

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
  [[nodiscard]] double
  gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) const override;

  /**
   * @brief Computes the log prior probabilities of assigning a data point to
   * every existing cluster. This method is useful for Gibbs sampling over
   * existing clusters. It returns a vector of log prior probabilities for all
   * existing clusters.
   * @param obs_idx The index of the observation to assign.
   * @return A vector of log prior probabilities for assigning the data point to
   * each existing cluster.
   */
  [[nodiscard]] Eigen::VectorXd
  gibbs_prior_existing_clusters(int obs_idx) const override;

  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   *
   * For NGGPW, this follows the NGGP formulation and is proportional to
   * alpha * sigma * (tau + U)^sigma.
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
   * NGGP-based split-merge MCMC algorithm.
   *
   * This method accounts for both the generalized gamma process prior and
   * spatial dependencies when computing the acceptance ratio for splitting
   * clusters.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */
  [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;

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
  [[nodiscard]] double prior_ratio_merge(int size_old_ci,
                                         int size_old_cj) const override;

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
  [[nodiscard]] double prior_ratio_shuffle(int size_old_ci, int size_old_cj,
                                           int ci, int cj) const override;

  /** @} */

  /**
   * @name Parameter Update Methods
   * @{
   */

  /**
   * @brief Updates the NGGPW parameters by updating the latent variable U.
   *
   * This method delegates the update to the U_sampler instance, which uses
   * an MCMC algorithm (RWMH or MALA) to sample U from its conditional
   * distribution given the current partition.
   *
   * @see U_sampler::update_U(), RWMH::update_U(), MALA::update_U()
   */
  void update_params() override { U_sampler_method.update_U(); };

  /** @} */
};