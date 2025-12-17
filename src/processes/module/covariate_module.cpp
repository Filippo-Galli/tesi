/**
 * @file covariate_module.cpp
 * @brief Implementation of `CovariatesModule`.
 */

#include "covariate_module.hpp"

CovariatesModule::ClusterStats
CovariatesModule::compute_cluster_statistics(const Eigen::Ref<const Eigen::VectorXi> obs) const {
    ClusterStats stats;

    for (int idx = 0; idx < obs.size(); ++idx) {
        const int obs_idx = obs(idx);
        const double age = covariates_data.ages(obs_idx);
        stats.n += 1;
        stats.sum += age;
        stats.sumsq += age * age;
    }

    return stats;
}

double CovariatesModule::compute_similarity_cls(int cls_idx, bool old_allo) const {
    ClusterStats stats;

    if (old_allo && old_allocations_provider) {
        const Eigen::VectorXi &allocations = old_allocations_provider();

        for (int i = 0; i < allocations.size(); ++i) {
            if (allocations(i) == cls_idx) {
                stats.n += 1;
                const double age = covariates_data.ages(i);
                stats.sum += age;
                stats.sumsq += age * age;
            }
        }
    } else {
        const Eigen::VectorXi &obs = data.get_cluster_assignments_ref(cls_idx);
        stats = compute_cluster_statistics(obs);
    }

    return log_marginal_likelihood_function(stats);
}

double CovariatesModule::compute_similarity_obs(int obs_idx, int cls_idx) const {

    ClusterStats base_stats = compute_cluster_statistics(data.get_cluster_assignments_ref(cls_idx));
    const double obs_age = covariates_data.ages(obs_idx);

    const double log_ml_without = log_marginal_likelihood_function(base_stats);

    base_stats.n += 1;
    base_stats.sum += obs_age;
    base_stats.sumsq += obs_age * obs_age;

    const double log_ml_with = log_marginal_likelihood_function(base_stats);
    return log_ml_with - log_ml_without;
}

Eigen::VectorXd CovariatesModule::compute_similarity_obs(int obs_idx) const {
    const Eigen::VectorXi &allocations = data.get_allocations();
    const int num_clusters = data.get_K();

    // Compute all cluster statistics once
    std::vector<ClusterStats> all_stats(num_clusters);
    for (int i = 0; i < allocations.size(); ++i) {
        int k = allocations(i);
        double val = covariates_data.ages(i);
        all_stats[k].n++;
        all_stats[k].sum += val;
        all_stats[k].sumsq += val * val;
    }

    // Compute log similarities for each cluster
    Eigen::VectorXd log_similarities(num_clusters);
    const double obs_age = covariates_data.ages(obs_idx);
    for (int k = 0; k < num_clusters; ++k) {
        ClusterStats &stats = all_stats[k];
        const double log_ml_without = log_marginal_likelihood_function(stats);

        stats.n += 1;
        stats.sum += obs_age;
        stats.sumsq += obs_age * obs_age;

        const double log_ml_with = log_marginal_likelihood_function(stats);

        log_similarities(k) = log_ml_with - log_ml_without;
    }

    return log_similarities;
}