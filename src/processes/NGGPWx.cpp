/**
 * @file NGGPWx.cpp
 * @brief Implementation of `NGGPWx` process.
 */

#include "./NGGPWx.hpp"

void NGGPWx::adjusted_cluster_statistics_cache() const {

    // Early exit if no changes
    if (NGGPW::data.changed_clusters.empty()) {
        return;
    }

    // Process any pending cluster changes to keep cache in sync
    for (const auto &change : NGGPW::data.changed_clusters) {
        const int index = change[0];
        const int old_cluster = change[1];
        const int new_cluster = change[2];

        // Special marker for full reallocation (index == -2)
        if (index == -2) {
            // Complete cache invalidation - clear everything
            CovariatesModule::cluster_stats_cache.clear();
            continue;
        }

        // Special events from compaction (index == -1):
        if (index == -1) {
            if (new_cluster == -1) {
                // Pure deletion: old_cluster was removed
                CovariatesModule::cluster_stats_cache.erase(old_cluster);
            } else {
                // Relabel: last cluster (new_cluster) was moved to old_cluster position
                auto it = CovariatesModule::cluster_stats_cache.find(new_cluster);
                if (it != CovariatesModule::cluster_stats_cache.end()) {
                    CovariatesModule::cluster_stats_cache[old_cluster] = std::move(it->second);
                    CovariatesModule::cluster_stats_cache.erase(new_cluster);
                } else {
                    // Last cluster had no cache entry; ensure old_cluster is empty too
                    CovariatesModule::cluster_stats_cache.erase(old_cluster);
                }
            }
            continue;
        }

        // Normal point movement: update old and new cluster stats
        // Pre-fetch age value once
        const double age = covariates_data.ages(index);

        // Remove from old cluster stats (if old_cluster is valid and cached)
        if (old_cluster >= 0) {
            auto old_it = CovariatesModule::cluster_stats_cache.find(old_cluster);
            if (old_it != CovariatesModule::cluster_stats_cache.end()) {
                auto &old_stats = old_it->second;
                old_stats.n -= 1;
                old_stats.sum -= age;
                old_stats.sumsq -= age * age;
                old_stats.invalidate();

                // Safety: if cluster is now empty, remove from cache
                if (old_stats.n <= 0) {
                    CovariatesModule::cluster_stats_cache.erase(old_cluster);
                }
            }
        }

        // Add to new cluster stats (if new_cluster is valid)
        if (new_cluster >= 0) {
            auto new_it = CovariatesModule::cluster_stats_cache.find(new_cluster);
            if (new_it != CovariatesModule::cluster_stats_cache.end()) {
                // Existing cluster in cache: update incrementally
                auto &new_stats = new_it->second;
                new_stats.n += 1;
                new_stats.sum += age;
                new_stats.sumsq += age * age;
                new_stats.invalidate();
            } else {
                // New cluster not in cache: compute full statistics
                CovariatesModule::ClusterStats new_stats;
                const auto &members = NGGPW::data.get_cluster_assignments_ref(new_cluster);
                const int size = members.size();
                for (int i = 0; i < size; ++i) {
                    const double member_age = covariates_data.ages(members(i));
                    new_stats.n += 1;
                    new_stats.sum += member_age;
                    new_stats.sumsq += member_age * member_age;
                }
                CovariatesModule::cluster_stats_cache.emplace(new_cluster, std::move(new_stats));
            }
        }
    }

    // Clear changes after adjustment
    NGGPW::data.changed_clusters.clear();
}

double NGGPWx::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {

    adjusted_cluster_statistics_cache();

    // NGGPW gibbs prior for existing cluster
    double log_prior = NGGPW::gibbs_prior_existing_cluster(cls_idx, obs_idx);

    // Add covariate module contribution
    log_prior += CovariatesModule::compute_similarity_obs(obs_idx, cls_idx);

    return log_prior;
}

Eigen::VectorXd NGGPWx::gibbs_prior_existing_clusters(int obs_idx) const {

    // Update cache only once at the start of Gibbs step
    adjusted_cluster_statistics_cache();

    // Get NGGPW gibbs prior for existing clusters
    Eigen::VectorXd log_prior = NGGPW::gibbs_prior_existing_clusters(obs_idx);

    // Add covariate module contributions
    for (int k = 0; k < log_prior.size(); ++k) {
        log_prior(k) += CovariatesModule::compute_similarity_obs(obs_idx, k);
    }

    return log_prior;
}

double NGGPWx::gibbs_prior_new_cluster() const { return NGGPW::gibbs_prior_new_cluster(); }

double NGGPWx::gibbs_prior_new_cluster_obs(int obs_idx) const {
    double log_prior = NGGPW::gibbs_prior_new_cluster();
    // Singleton covariate similarity: log g({x_i}) - log g(∅) = log g({x_i})
    log_prior += CovariatesModule::compute_similarity_obs(obs_idx, -1);
    return log_prior;
}

double NGGPWx::prior_ratio_split(int ci, int cj) const {
    // NGGPW prior ratio for split
    double log_prior_ratio = NGGPW::prior_ratio_split(ci, cj);

    // Covariate similarity ratio: sim(new ci) + sim(new cj) - sim(old merged ci)
    log_prior_ratio += CovariatesModule::compute_similarity_cls(ci, false);
    log_prior_ratio += CovariatesModule::compute_similarity_cls(cj, false);
    log_prior_ratio -= CovariatesModule::compute_similarity_cls(ci, true);

    return log_prior_ratio;
}

double NGGPWx::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
    // NGGPW prior ratio for merge
    double log_prior_ratio = NGGPW::prior_ratio_merge(size_old_ci, size_old_cj);

    // Add covariate module contributions
    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    log_prior_ratio -= CovariatesModule::compute_similarity_cls(old_ci, true);
    log_prior_ratio -= CovariatesModule::compute_similarity_cls(old_cj, true);

    const int new_ci = NGGPW::data.get_allocations()[idx_i];
    log_prior_ratio += CovariatesModule::compute_similarity_cls(new_ci, false);

    return log_prior_ratio;
}

double NGGPWx::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) const {
    // NGGPW prior ratio for shuffle
    double log_prior_ratio = NGGPW::prior_ratio_shuffle(size_old_ci, size_old_cj, ci, cj);

    // Add covariate module contributions
    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    log_prior_ratio -= CovariatesModule::compute_similarity_cls(old_ci, true);
    log_prior_ratio -= CovariatesModule::compute_similarity_cls(old_cj, true);

    const int new_ci = NGGPW::data.get_allocations()[idx_i];
    const int new_cj = NGGPW::data.get_allocations()[idx_j];
    log_prior_ratio += CovariatesModule::compute_similarity_cls(new_ci, false);
    log_prior_ratio += CovariatesModule::compute_similarity_cls(new_cj, false);

    return log_prior_ratio;
}
