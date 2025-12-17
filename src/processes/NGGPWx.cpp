/**
 * @file NGGPWx.cpp
 * @brief Implementation of `NGGPWx` process.
 */

#include "./NGGPWx.hpp"

double NGGPWx::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {

    // NGGPW gibbs prior for existing cluster
    double log_prior = NGGPW::gibbs_prior_existing_cluster(cls_idx, obs_idx);

    // Add covariate module contribution
    log_prior += CovariatesModule::compute_similarity_obs(obs_idx, cls_idx);

    return log_prior;
}

Eigen::VectorXd NGGPWx::gibbs_prior_existing_clusters(int obs_idx) const {

    // Get NGGPW gibbs prior for existing clusters
    Eigen::VectorXd log_prior = NGGPW::gibbs_prior_existing_clusters(obs_idx);

    // Add covariate module contributions
    log_prior += CovariatesModule::compute_similarity_obs(obs_idx);

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
