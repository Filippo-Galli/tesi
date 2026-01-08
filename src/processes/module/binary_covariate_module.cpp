/**
 * @file binary_covariate_module.cpp
 * @brief Implementation of `BinaryCovariatesModule`.
 */

#include "binary_covariate_module.hpp"

double BinaryCovariatesModule::compute_similarity_cls(int cls_idx, bool old_allo) const {
    // Retrieve cluster members based on allocation type
    Eigen::VectorXi cluster_members;
    if (old_allo && old_cluster_members_provider) {
        auto it = old_cluster_members_provider->find(cls_idx);
        if (it != old_cluster_members_provider->end()) {
            cluster_members = Eigen::Map<const Eigen::VectorXi>(it->second.data(), it->second.size());
        }
    } else {
        cluster_members = data.get_cluster_assignments(cls_idx);
    }

    // If the cluster is empty, return 0 similarity
    if (cluster_members.size() == 0) {
        return 0.0;
    }

    // Initialize counts for binary covariates
    int num_covariates = cluster_members.size();
    int counts = 0;

    // Count occurrences of '1' in the binary covariates for the cluster members
    for (int i = 0; i < num_covariates; ++i) {
        int obs_idx = cluster_members(i);
        counts += binary_covariate_data(obs_idx);
    }

    // Compute log marginal likelihood using Beta-Binomial model
    double alpha = beta_prior_alpha;
    double beta = beta_prior_beta;

    double log_marginal_likelihood = log_beta_prior + lgamma(counts + alpha) +
                                     lgamma(num_covariates - counts + beta) -
                                     lgamma(num_covariates + alpha + beta);

    return log_marginal_likelihood;
}

double BinaryCovariatesModule::compute_similarity_obs(int obs_idx, int cls_idx) const {
    
    const auto & cluster_members = data.get_cluster_assignments(cls_idx);

    // Initialize counts for binary covariates
    int num_covariates = cluster_members.size();
    int counts = 0;

    // Count occurrences of '1' in the binary covariates for the cluster members
    for (int i = 0; i < num_covariates; ++i) {
        int obs_idx = cluster_members(i);
        counts += binary_covariate_data(obs_idx);
    }

    const double denominator = num_covariates + beta_prior_alpha + beta_prior_beta;
    const double numerator = (beta_prior_alpha + counts)*binary_covariate_data(obs_idx) + 
                             (beta_prior_beta + num_covariates - counts)*(1 - binary_covariate_data(obs_idx));
    return std::log(numerator / denominator);
}

Eigen::VectorXd BinaryCovariatesModule::compute_similarity_obs(int obs_idx) const{
    Eigen::VectorXd similarities(data.get_K());
    for (int cls_idx = 0; cls_idx < data.get_K(); ++cls_idx) {
        similarities(cls_idx) = compute_similarity_obs(obs_idx, cls_idx);
    }
    return similarities;
}
