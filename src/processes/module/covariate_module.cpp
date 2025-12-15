/**
 * @file covariate_module.cpp
 * @brief Implementation of `CovariatesModule`.
 */

#include "covariate_module.hpp"

#include <cmath>

CovariatesModule::ClusterStats CovariatesModule::compute_cluster_statistics(int cls_idx,
                                                                            const Eigen::VectorXi &allocations) const {
    ClusterStats stats;

    if (cls_idx == -1) {
        return stats;
    }

    for (int i = 0; i < data.get_n(); ++i) {
        if (allocations(i) == cls_idx) {
            const double age = covariates_data.ages(i);
            stats.n += 1;
            stats.sum += age;
            stats.sumsq += age * age;
        }
    }

    return stats;
}

double CovariatesModule::compute_log_marginal_likelihood(const ClusterStats &stats) const {
    if (stats.n == 0) {
        return 0.0;
    }

    const int n = stats.n;

    const double xbar = stats.sum / static_cast<double>(n);
    const double ss = stats.sumsq - static_cast<double>(n) * xbar * xbar;

    const double v_nB = v + static_cast<double>(n) * B;
    const double tau_j = Bv / v_nB;

    double log_ml = 0.0;
    log_ml += static_cast<double>(n) * const_term;
    log_ml += -0.5 * static_cast<double>(n) * log_v;
    log_ml += 0.5 * (std::log(tau_j) - log_B);
    log_ml -= 0.5 * ss / v;

    const double prior_deviation = static_cast<double>(n) * (xbar - m) * (xbar - m) / v_nB;
    log_ml -= 0.5 * prior_deviation;

    return log_ml;
}

double CovariatesModule::compute_similarity_cls(int cls_idx, bool old_allo) const {
    const Eigen::VectorXi &allocations =
        (old_allo && old_allocations_provider) ? old_allocations_provider() : data.get_allocations();

    const ClusterStats stats = compute_cluster_statistics(cls_idx, allocations);
    return compute_log_marginal_likelihood(stats);
}

double CovariatesModule::compute_similarity_obs(int obs_idx, int cls_idx, bool old_allo) const {
    const Eigen::VectorXi &allocations =
        (old_allo && old_allocations_provider) ? old_allocations_provider() : data.get_allocations();

    const ClusterStats base_stats = compute_cluster_statistics(cls_idx, allocations);
    const double obs_age = covariates_data.ages(obs_idx);

    const double log_ml_without = compute_log_marginal_likelihood(base_stats);

    ClusterStats with_stats = base_stats;
    with_stats.n += 1;
    with_stats.sum += obs_age;
    with_stats.sumsq += obs_age * obs_age;

    const double log_ml_with = compute_log_marginal_likelihood(with_stats);
    return log_ml_with - log_ml_without;
}