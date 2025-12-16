/**
 * @file covariate_module.hpp
 * @brief Covariate-related computations for clustering processes.
 */

#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Covariates.hpp"
#include "Eigen/Dense"
#include <functional>
#include <utility>
#include <cmath>

/**
 * @class CovariatesModule
 * @brief Module for covariate-related computations within clustering processes.
 *
 * This class implements the product partition model with regression on covariates
 * as described in Müller et al. (2011). It computes similarity measures based on
 * how well observations within a cluster can be explained by a common covariate
 * distribution (Normal conjugate prior).
 *
 * Reference: Müller, P., Quintana, F. (2011)
 * "A Product Partition Model With Regression on Covariates"
 */
class CovariatesModule {
protected:
    /**
     * @name Module References
     * @{
     */

    /** @brief Reference to covariates data containing ages and prior parameters */
    const Covariates &covariates_data;

    /** @brief Reference to data object with cluster assignments */
    const Data &data;

    /**
     * @brief Provider function for accessing old allocation state.
     *
     * Used when computing similarity based on previous cluster assignments
     * (e.g., in split-merge algorithms).
     */
    std::function<const Eigen::VectorXi &()> old_allocations_provider;

    /** @} */

    /**
     * @name Precomputed values
     * @{
     */    

    const double Bv = covariates_data.B * covariates_data.v; ///< Product of prior variance and observation variance

    const double log_B = std::log(covariates_data.B); ///< Log of prior variance
    const double log_v = std::log(covariates_data.v); ///< Log of observation variance

    const double const_term = -0.5 * std::log(2.0 * M_PI); ///< Constant term in log likelihood

    /** @} */

    /**
     * @name Helper Methods
     * @{
     */

    /**
     * @brief Sufficient statistics for covariate likelihood computations.
     */
    struct ClusterStats {
        int n = 0;
        double sum = 0.0;
        double sumsq = 0.0;
    };

    /**
     * @brief Compute cluster statistics for covariate similarity
     *
     * @param cls_idx Index of the cluster
     * @param allocations Current allocation vector
     * @return Sufficient statistics (n, sum, sum of squares)
     */
    ClusterStats compute_cluster_statistics(int cls_idx, const Eigen::VectorXi &allocations) const;

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal conjugate prior model:
     * - Prior on mean μ: N(m, B)
     *
     * @param stats Sufficient statistics (n, sum, sum of squares)
     * @return Log marginal likelihood contribution
     */
    double compute_log_marginal_likelihood_NN(const ClusterStats &stats) const;

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal-InverseGamma conjugate prior model:
     * - Prior on mean μ: N(m, B)
     * - Prior on variance σ²: IG(nu, S0)
     *
     * @param stats Sufficient statistics (n, sum, sum of squares)
     * @return Log marginal likelihood contribution
     */
    double compute_log_marginal_likelihood_NNIG(const ClusterStats &stats) const;

    /** @} */

public:
    /**
     * @brief Constructor for CovariatesModule
     *
     * @param covariates_ Reference to Covariates object with age data and priors
     * @param data_ Reference to Data object with cluster assignments
     * @param old_alloc_provider Optional function to access old allocations
     */
    CovariatesModule(const Covariates &covariates_, const Data &data_,
                     std::function<const Eigen::VectorXi &()> old_alloc_provider = {})
        : covariates_data(covariates_), data(data_), old_allocations_provider(std::move(old_alloc_provider)) {}

    /**
     * @name Similarity Computation Methods
     * @{
     */

    /**
     * @brief Compute covariate similarity contribution for a cluster
     *
     * Computes the log marginal likelihood of the covariates within a cluster
     * under the Normal conjugate model. Higher values indicate
     * that observations in the cluster have similar covariate values.
     *
     * @param cls_idx Index of the cluster (0 to K-1)
     * @param old_allo If true, uses old allocations from old_allocations_provider;
     *                 if false, uses current allocations from data (default: false)
     * @return Log marginal likelihood contribution (similarity score)
     *
     * @details The computation follows Müller et al. (2011):
     * 1. Compute sufficient statistics (n, sum, sum of squares)
     * 2. Update hyperparameters using conjugate update rules
     * 3. Compute log marginal likelihood using updated parameters
     *
     * This value is added to the clustering prior in split-merge moves
     * to encourage clusters with homogeneous covariate values.
     */
    double compute_similarity_cls(int cls_idx, bool old_allo = false) const;

    /**
     * @brief Compute covariate similarity for a single observation in a cluster
     *
     * Computes the predictive contribution when adding observation obs_idx to
     * cluster cls_idx, considering the covariate values.
     *
     * @param obs_idx Index of the observation
     * @param cls_idx Index of the cluster
     * @param old_allo If true, uses old allocations (default: false)
     * @return Log predictive density contribution
     *
     * @details Used in Gibbs sampling to compute the probability of assigning
     * an observation to a cluster based on covariate similarity.
     */
    double compute_similarity_obs(int obs_idx, int cls_idx, bool old_allo = false) const;

    /** @} */
};