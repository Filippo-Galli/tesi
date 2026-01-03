/**
 * @file covariate_module.hpp
 * @brief Covariate-related computations for clustering processes.
 */

#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Covariates.hpp"
#include "../../utils/Module.hpp"
#include "Eigen/Dense"
#include <cmath>
#include <vector>

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
class CovariatesModule : public Module {
protected:
    /**
     * @name Module References
     * @{
     */

    /** @brief Reference to covariates data containing ages and prior parameters */
    const Covariates &covariates_data;

    /** @brief Reference to data object with cluster assignments */
    const Data &data;

    /** @brief Provider function for accessing old cluster members map */
    const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider;

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
    ClusterStats compute_cluster_statistics(const Eigen::Ref<const Eigen::VectorXi> obs) const;

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal-Normal conjugate prior model:
     * - x_i ~ N(μ_j, v)  for i ∈ S_j
     * - Prior on mean μ_j: N(m, B)
     * - Observation variance v is known and fixed
     *
     * @details The marginal likelihood integrates out the cluster-specific mean μ_j.
     * With sufficient statistics:
     * - n_j = |S_j| (cluster size)
     * - x̄_j = (1/n_j) Σ_{i ∈ S_j} x_i (sample mean)
     * - SS = Σ_{i ∈ S_j} (x_i - x̄_j)² (centered sum of squares)
     *
     * The posterior distribution of μ_j is N(m̂_j, τ_j) where:
     * - τ_j = Bv / (v + n_j B)  (posterior variance)
     * - m̂_j = τ_j (n_j x̄_j / v + m / B)  (posterior mean)
     *
     * The log marginal likelihood is:
     *
     * log q(x_j) = -n_j/2 log(2π) - n_j/2 log(v) - 1/2 log(B) + 1/2 log(τ_j)
     *              - SS/(2v) - n_j(x̄_j - m)² / (2(v + n_j B))
     *
     * where log(τ_j) = log(B) + log(v) - log(v + n_j B)
     *
     *
     * @param stats Sufficient statistics for the cluster
     * @return Log marginal likelihood value
     *
     * @note This is marked as __attribute__((hot)) for performance optimization
     *       as it is called frequently in the MCMC sampling loop.
     */
    double compute_log_marginal_likelihood_NN(const ClusterStats &stats) const __attribute__((hot));

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal-InverseGamma conjugate prior model:
     * x ~ N(μ, v_j)
     * - Prior on mean μ: N(m, B*v_j)
     * - Prior on variance v_j: IG(nu, S0)
     *
     * @param stats Sufficient statistics (n, sum, sum of squares)
     * @return Log marginal likelihood contribution
     *
     * @details The marginal likelihood for the NNIG model is:
     * log g(S) = log Γ(ν + n/2) - log Γ(ν) - n/2 log(2π)
     *            - 1/2 log(1 + nB) + ν log(S₀)
     *            - (ν + n/2) log(S₀ + SS/2 + n/(2(1+nB)) (x̄-m)²)
     */
    double compute_log_marginal_likelihood_NNIG(const ClusterStats &stats) const __attribute__((hot));

    /**
     * @brief Compute log marginal likelihood based on model type
     *
     * Chooses between NN and NNIG models based on covariates_data.fixed_v.
     *
     * @param stats Sufficient statistics for the cluster
     * @return Log marginal likelihood value
     */
    inline double compute_log_marginal_likelihood(const ClusterStats &stats) const __attribute__((hot)) {
        if (covariates_data.fixed_v) {
            return compute_log_marginal_likelihood_NN(stats);
        } else {
            return compute_log_marginal_likelihood_NNIG(stats);
        }
    }

    /** @} */

    /**
     * @name Precomputed values
     * @{
     */

    const double Bv; ///< Product of prior variance and observation variance

    const double log_B; ///< Log of prior variance
    const double log_v; ///< Log of observation variance

    const double const_term; ///< Constant term in log likelihood

    const double lgamma_nu; ///< log Gamma(ν) for NNIG model (v ~ IG(ν, S₀))
    const double nu_logS0;  ///< ν log(S₀) for NNIG model (v ~ IG(ν, S₀))

    std::vector<double> log_v_plus_nB; ///< Cache for log(v_plus_nB) for NN
    std::vector<double> lgamma_nu_n;   ///< Cache for lgamma(nu_n) for NNIG

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
                     const Eigen::VectorXi *old_alloc_provider = nullptr,
                     const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : covariates_data(covariates_), data(data_), Module(old_alloc_provider),
          old_cluster_members_provider(old_cluster_members_provider_),
          // Initialize constants here in the list
          Bv(covariates_.B * covariates_.v), log_B(std::log(covariates_.B)), log_v(std::log(covariates_.v)),
          const_term(-0.5 * std::log(2.0 * M_PI)), lgamma_nu(std::lgamma(covariates_.nu)),
          nu_logS0(covariates_.nu * std::log(covariates_.S0)) {

        // Precompute caches for efficiency if needed
        if (covariates_data.fixed_v) {
            log_v_plus_nB.reserve(data_.get_n() + 1);
            for (int n = 0; n <= data_.get_n(); ++n) {
                log_v_plus_nB.push_back(std::log(covariates_data.v + n * covariates_data.B));
            }
        } else {
            lgamma_nu_n.reserve(data_.get_n() + 1);
            for (int n = 0; n <= data_.get_n(); ++n) {
                lgamma_nu_n.push_back(std::lgamma(covariates_data.nu + 0.5 * static_cast<double>(n)));
            }
        }
    }

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
    double compute_similarity_cls(int cls_idx, bool old_allo = false) const override __attribute__((hot));

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
    double compute_similarity_obs(int obs_idx, int cls_idx) const override __attribute__((hot));

    /**
     * @brief Compute covariate similarity contributions for all existing clusters
     *
     * Computes the predictive contributions for adding observation obs_idx
     * to each existing cluster, considering covariate values.
     *
     * @param obs_idx Index of the observation
     * @return Vector of log predictive density contributions for each cluster
     */
    Eigen::VectorXd compute_similarity_obs(int obs_idx) const override __attribute__((hot));

    /** @} */
};