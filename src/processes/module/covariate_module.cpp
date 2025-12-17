/**
 * @file covariate_module.cpp
 * @brief Implementation of `CovariatesModule`.
 */

#include "covariate_module.hpp"

#include <cmath>

CovariatesModule::ClusterStats CovariatesModule::compute_cluster_statistics(const std::vector<int> &obs) const {
    ClusterStats stats;

    for (const auto &i : obs) {
        const double age = covariates_data.ages(i);
        stats.n += 1;
        stats.sum += age;
        stats.sumsq += age * age;
    }

    return stats;
}

double CovariatesModule::compute_similarity_cls(int cls_idx, bool old_allo) const {

    // Try to use cache first (only for current allocations, not old)
    if (!old_allo) {
        auto it = cluster_stats_cache.find(cls_idx);
        if (it != cluster_stats_cache.end()) {
            // Cache hit - check if log ML is already computed
            ClusterStats &stats = it->second;
            if (!stats.log_ml_valid) {
                stats.cached_log_ml = covariates_data.fixed_v ? compute_log_marginal_likelihood_NN(stats)
                                                              : compute_log_marginal_likelihood_NNIG(stats);
                stats.log_ml_valid = true;
            }
            return stats.cached_log_ml;
        }
    }

    // Cache miss or old allocations requested - compute from scratch
    ClusterStats stats;

    if (old_allo && old_allocations_provider) {
        // Use old allocations from provider
        const Eigen::VectorXi &old_allocs = old_allocations_provider();
        for (int i = 0; i < old_allocs.size(); ++i) {
            if (old_allocs(i) == cls_idx) {
                const double age = covariates_data.ages(i);
                stats.n += 1;
                stats.sum += age;
                stats.sumsq += age * age;
            }
        }
    } else {
        // Use current allocations from data
        const auto &assignments = data.get_cluster_assignments_ref(cls_idx);
        for (int i = 0; i < assignments.size(); ++i) {
            const double age = covariates_data.ages(assignments(i));
            stats.n += 1;
            stats.sum += age;
            stats.sumsq += age * age;
        }
    }

    return covariates_data.fixed_v ? compute_log_marginal_likelihood_NN(stats)
                                   : compute_log_marginal_likelihood_NNIG(stats);
}

double CovariatesModule::compute_similarity_obs(int obs_idx, int cls_idx) const {

    const double obs_age = covariates_data.ages(obs_idx);
    const double obs_age_sq = obs_age * obs_age;

    // Single cache lookup
    auto it = cluster_stats_cache.find(cls_idx);
    if (it != cluster_stats_cache.end()) {
        // Cache hit - ensure log ML is computed
        ClusterStats &stats = it->second;
        if (!stats.log_ml_valid) {
            stats.cached_log_ml = covariates_data.fixed_v ? compute_log_marginal_likelihood_NN(stats)
                                                          : compute_log_marginal_likelihood_NNIG(stats);
            stats.log_ml_valid = true;
        }

        // Compute stats with observation added (inline for speed)
        ClusterStats stats_with;
        stats_with.n = stats.n + 1;
        stats_with.sum = stats.sum + obs_age;
        stats_with.sumsq = stats.sumsq + obs_age_sq;

        const double log_ml_with = covariates_data.fixed_v ? compute_log_marginal_likelihood_NN(stats_with)
                                                           : compute_log_marginal_likelihood_NNIG(stats_with);
        return log_ml_with - stats.cached_log_ml;
    }

    // Cache miss - singleton cluster (log_ml_without = 0)
    ClusterStats stats_with;
    stats_with.n = 1;
    stats_with.sum = obs_age;
    stats_with.sumsq = obs_age_sq;

    return covariates_data.fixed_v ? compute_log_marginal_likelihood_NN(stats_with)
                                   : compute_log_marginal_likelihood_NNIG(stats_with);
}