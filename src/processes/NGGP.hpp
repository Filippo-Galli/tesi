/**
 * @file NGGP.hpp
 * @brief Normalized Generalized Gamma Process (NGGP) implementation for
 * Bayesian nonparametric clustering.
 */

#pragma once

#include "../utils/Process.hpp"
#include <deque>

/**
 * @class NGGP
 * @brief Normalized Generalized Gamma Process class for Bayesian nonparametric
 * clustering.
 *
 * This class implements a Normalized Generalized Gamma Process (NGGP) that
 * extends the Dirichlet Process by incorporating a latent variable U. The NGGP
 * provides more flexibility in modeling cluster sizes and incorporates adaptive
 * behavior through the U parameter.
 */
class NGGP : public Process {

private:
  /**
   * @name Private Member Variables
   * @{
   */

  /** @brief Latent variable U in the NGGP model. */
  double U = 1;

  /** @brief Parameter tau from the NGGP specification. */
  double tau = params.tau;

  /** @} */

  /**
   * @name Private Methods
   * @{
   */

  /**
   * @brief Updates the latent variable U using Metropolis-Hastings.
   */
  void update_U();

  /**
   * @brief Computes the log conditional density of U given the partition.
   * @param u The value of U.
   * @return The log of the unnormalized conditional density.
   */
  double log_conditional_density_U(double u) const;

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

  /** @brief Standard deviation for the Gaussian proposal distribution. */
  const double proposal_std = 0.5;

  /**
   * @name Cached Constants
   * @brief Pre-computed values for computational efficiency.
   * @{
   */

  /** @brief Ratio a/sigma for efficient computation. */
  const double a_over_sigma = params.a / params.sigma;

  /** @brief Pre-computed tau^sigma for efficient computation. */
  const double tau_power_sigma = std::pow(params.tau, params.sigma);

  /** @} */

  /**
   * @name Statistics
   * @{
   */

  /** @brief Counter for accepted U updates (for monitoring acceptance rate). */
  mutable int accepted_U = 0;

  /** @brief total number of iteration */
  int total_iterations = 0;

  /** @} */

public:
  /**
   * @brief Constructor for the Normalized Generalized Gamma Process.
   * @param d Reference to the data object containing observations and cluster
   * assignments.
   * @param p Reference to the parameters object containing NGGP parameters (a,
   * sigma, tau).
   */
  NGGP(Data &d, Params &p) : Process(d, p), gen(rd()) {};

  /**
   * @name Gibbs Sampling Methods
   * @{
   */

  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   *
   * For NGGP, this incorporates the discount parameter sigma, giving
   * probability proportional to (n_k - sigma) where n_k is the cluster size.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation (default: 0, unused in this
   * implementation).
   * @return The log prior probability of assigning the data point to the
   * existing cluster.
   */
  [[nodiscard]] double gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) const override;

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
   *
   * For NGGP, this depends on the latent variable U and is proportional to
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
   * @brief Computes the prior ratio for a split operation in an NGGP-based
   * split-merge MCMC algorithm.
   *
   * This method accounts for the generalized gamma process prior when computing
   * the acceptance ratio for splitting clusters.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */
  [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;

  /**
   * @brief Computes the prior ratio for a merge operation in an NGGP-based
   * split-merge MCMC algorithm.
   *
   * This method accounts for the generalized gamma process prior when computing
   * the acceptance ratio for merging clusters.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */
  [[nodiscard]] double prior_ratio_merge(int size_old_ci,
                                         int size_old_cj) const override;

  /**
   * @brief Computes the prior ratio for a shuffle operation in an NGGP-based
   * split-merge MCMC algorithm.
   *
   * This method accounts for the generalized gamma process prior when computing
   * the acceptance ratio for shuffling observations between clusters.
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
   * @brief Updates the NGGP parameters by updating the latent variable U.
   */
  void update_params() override { update_U(); };

  /** @} */

  /**
   * @name Accessor Methods
   * @{
   */

  /**
   * @brief Gets the current value of the latent variable U.
   * @return The current value of U.
   */
  double get_U() const { return U; }

  /**
   * @brief Gets the number of accepted U updates for monitoring convergence.
   * @return The number of accepted U updates.
   */
  int get_accepted_U() const { return accepted_U; }

  /**
   * @brief Gets the current acceptance rate for U updates.
   * @return The acceptance rate as a fraction in [0, 1].
   */
  double get_acceptance_rate() const {
    return total_iterations > 0 ?  static_cast<double>(accepted_U) / total_iterations : 0.0;
  }

  /** @} */
};