/**
 * @file neal.hpp
 * @brief Neal's Algorithm 3 implementation for collapsed Gibbs sampling
 *
 * This file contains the implementation of Neal's Algorithm 3, a standard
 * collapsed Gibbs sampler for Bayesian nonparametric mixture models. The
 * algorithm provides efficient sequential updates for cluster assignments with
 * automatic cluster creation and deletion.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"

/**
 * @brief Implementation of Neal's Algorithm 3 for collapsed Gibbs sampling
 *
 * This class implements Neal's Algorithm 3, a standard collapsed Gibbs sampler
 * for Bayesian nonparametric mixture models. The algorithm sequentially updates
 * cluster assignments by sampling from the full conditional distribution of
 * each observation, integrating out the cluster parameters.
 *
 * @details Algorithm 3 is characterized by:
 * - **Sequential updates**: One observation at a time
 * - **Collapsed sampling**: Cluster parameters are integrated out analytically
 * - **Full conditionals**: Each assignment is sampled from its exact posterior
 * - **Automatic cluster creation/deletion**: Clusters are created when needed
 * and deleted when empty
 *
 * The algorithm is particularly effective for models where the likelihood can
 * be computed in closed form after integrating out cluster-specific parameters.
 *
 * @reference Neal, R. M. (2000). "Markov Chain Sampling Methods for Dirichlet
 * Process Mixture Models"
 *
 * @see Sampler, Process, Likelihood
 */
class Neal3ZDNAM : public Sampler {
private:
  // ========== Random Number Generation ==========

  /** @brief Mersenne Twister random number generator for sampling operations */
  mutable std::mt19937 gen;

  // ========== Core Algorithm Methods ==========

  /**
   * @brief Sample cluster assignment for a single observation
   *
   * @param index Index of the observation to update
   *
   * @details This method implements the core of Algorithm 3:
   * 1. Remove the observation from its current cluster
   * 2. Compute probabilities for all existing clusters plus a new cluster
   * 3. Sample new assignment from the full conditional distribution
   * 4. Update cluster assignments and clean up empty clusters
   *
   * The probabilities combine prior information from the Process with
   * likelihood information computed by integrating out cluster parameters.
   */
  void step_1_observation(int index);

  int sample_from_log_probs(const std::vector<double>& log_probs);

  std::vector<double> compute_zdnam_probabilities(const std::vector<double>& log_probs, int current_k);

  int sample_from_log_probs_zdnam( const std::vector<double>& log_probs, int current_k);

public:
  // ========== Constructor ==========

  /**
   * @brief Constructor for Neal's Algorithm 3 sampler
   *
   * @param d Reference to Data object containing observations
   * @param p Reference to Params object with hyperparameters
   * @param l Reference to Likelihood object for probability computations
   * @param pr Reference to Process object (DP, NGGP, etc.) defining the prior
   *
   * @details Initializes the Gibbs sampler with all required components.
   * The random number generator is seeded from the inherited random device.
   */
  Neal3ZDNAM(Data &d, Params &p, Likelihood &l, Process &pr)
      : Sampler(d, p, l, pr), gen(rd()) {};

  // ========== MCMC Interface ==========

  /**
   * @brief Perform one complete iteration of Neal's Algorithm 3
   *
   * @details Executes one full sweep of the collapsed Gibbs sampler by
   * sequentially updating the cluster assignment of each observation.
   * The order of updates is typically randomized to avoid systematic bias.
   *
   * After each full sweep:
   * - All observations have been considered for reassignment
   * - Cluster structure may have changed (clusters created/destroyed)
   * - The Markov chain has advanced by one step
   *
   * @see step_1_observation()
   */
  void step() override;
};