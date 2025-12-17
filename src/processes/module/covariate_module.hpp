/**
 * @file covariate_module.hpp
 * @brief Covariate-related computations for clustering processes.
 */

#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Covariates.hpp"
#include "Eigen/Dense"
#include <functional>
#include <unordered_map>
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

    const double lgamma_nu = std::lgamma(covariates_data.nu); ///< log Gamma(ν) for NNIG model (v ~ IG(ν, S₀))
    const double nu_logS0 =
        covariates_data.nu * std::log(covariates_data.S0); ///< ν log(S₀) for NNIG model (v ~ IG(ν, S₀))

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
        mutable double cached_log_ml = 0.0; ///< Cached log marginal likelihood
        mutable bool log_ml_valid = false;  ///< Whether cached log ML is valid

        // Pre-compute expensive operations for incremental update
        inline void invalidate() { log_ml_valid = false; }
    };

    mutable ClusterStats empty_stats; ///< Empty cluster statistics

    mutable std::unordered_map<int, ClusterStats> cluster_stats_cache; ///< Cache for cluster statistics

    /**
     * @brief Compute cluster statistics for covariate similarity
     *
     * @param cls_idx Index of the cluster
     * @param allocations Current allocation vector
     * @return Sufficient statistics (n, sum, sum of squares)
     */
    ClusterStats compute_cluster_statistics(const std::vector<int> &obs) const;

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal conjugate prior model:
     * - x ~ N(μ, v)
     * - Prior on mean μ: N(m, B)
     * - Known variance v
     *
     * @param stats Sufficient statistics (n, sum, sum of squares)
     * @return Log marginal likelihood contribution
     *
     * @details Following Müller et al. (2011), the marginal likelihood is:
     * log g(S) = -n/2 log(2π) - n/2 log(v) + 1/2 log(B) - 1/2 log(B + nv)
     *            - SS/(2v) - n(x̄ - m)²/(2(B + nv))
     * where SS = Σ(xᵢ - x̄)²
     */
    inline double compute_log_marginal_likelihood_NN(const ClusterStats &stats) const __attribute__((hot)) {
        if (stats.n == 0) {
            return 0.0;
        }

        const double n_dbl = static_cast<double>(stats.n);
        const double inv_n = 1.0 / n_dbl;
        const double xbar = stats.sum * inv_n;

        // Centered sum of squares: SS = Σ(xᵢ - x̄)²
        const double ss = stats.sumsq - n_dbl * xbar * xbar;

        // Posterior variance parameter: B + nv (NOT v + nB!)
        const double B_nv = covariates_data.B + n_dbl * covariates_data.v;

        // Deviation from prior mean
        const double dev = xbar - covariates_data.m;

        // log g(S) = -n/2 log(2π) - n/2 log(v) + 1/2 log(B) - 1/2 log(B + nv)
        //            - SS/(2v) - n(x̄ - m)²/(2(B + nv))
        return n_dbl * (const_term - 0.5 * log_v) + 0.5 * (log_B - std::log(B_nv)) -
               0.5 * (ss / covariates_data.v + n_dbl * dev * dev / B_nv);
    }

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
     *            - 1/2 log(1 + nB) + ν log(S₀)              ← Fix: remove /2
     *            - (ν + n/2) log(S₀ + SS/2 + n/(2(1+nB)) (x̄-m)²)
     */
    inline double compute_log_marginal_likelihood_NNIG(const ClusterStats &stats) const __attribute__((hot)) {
        if (stats.n == 0) {
            return 0.0;
        }

        const double n_dbl = static_cast<double>(stats.n);
        const double inv_n = 1.0 / n_dbl;
        const double xbar = stats.sum * inv_n;

        // Centered sum of squares: SS = Σ(xᵢ - x̄)²
        const double ss = stats.sumsq - n_dbl * xbar * xbar;

        // Posterior parameters
        // Prior: v ~ IG(ν, S0) and μ|v ~ N(m, B v) with B interpreted as a variance multiplier.
        const double nu_n = covariates_data.nu + 0.5 * n_dbl;

        // Deviation from prior mean
        const double dev = xbar - covariates_data.m;

        // Posterior scale parameter:
        // S_n = S₀ + SS/2 + n/(2(1+nB)) (x̄-m)²
        const double one_plus_nB = 1.0 + n_dbl * covariates_data.B;
        const double S_n = covariates_data.S0 + 0.5 * ss + 0.5 * (n_dbl / one_plus_nB) * dev * dev;

        // log g(S) = log Γ(ν + n/2) - log Γ(ν) - n/2 log(2π)
        //            - 1/2 log(1+nB) + ν log(S₀) - (ν + n/2) log(S_n)
        return std::lgamma(nu_n) - lgamma_nu + n_dbl * const_term - 0.5 * std::log(one_plus_nB) + nu_logS0 -
               nu_n * std::log(S_n);
    }

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
    double compute_similarity_obs(int obs_idx, int cls_idx) const __attribute__((hot));

    /** @} */
};