#pragma once

#include "NGGPW.hpp"
#include "module/covariate_module.hpp"

class NGGPWx : public NGGPW, protected CovariatesModule {
public:
  NGGPWx(Data &d, Params &p, Covariates &cov, U_sampler &U_sam)
      : NGGPW(d, p, cov, U_sam),
        CovariatesModule(cov, d, [this]() -> const Eigen::VectorXi & {
          return this->old_allocations_view();
        }) {}

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

  /**
   * @brief Observation-specific new-cluster prior including covariate singleton
   * similarity.
   */
  [[nodiscard]] double gibbs_prior_new_cluster_obs(int obs_idx) const override;

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
  void update_params() override { NGGPW::U_sampler_method.update_U(); };

  /** @} */
};