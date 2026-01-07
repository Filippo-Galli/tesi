/**
 * @file covariate_module_cache.cpp
 * @brief Implementation of `CovariatesModuleCache`.
 */

#include "covariate_module_cache.hpp"

ClusterStats CovariatesModuleCache::compute_cluster_statistics(const Eigen::Ref<const Eigen::VectorXi> obs) const {
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

double CovariatesModuleCache::compute_similarity_cls(int cls_idx, bool old_allo) const {

    if (old_allo) {
        ClusterStats stats;
        const auto &old_cls_allo = old_cluster_members_provider->at(cls_idx);
        stats = compute_cluster_statistics(Eigen::Map<const Eigen::VectorXi>(old_cls_allo.data(), old_cls_allo.size()));
        return compute_log_marginal_likelihood(stats);

    } else {
        const Eigen::VectorXi &obs = data.get_cluster_assignments_ref(cls_idx);
        auto stats = covariate_cache.get_cluster_stats_ref(cls_idx);
        return compute_log_marginal_likelihood(stats);
    }
}

double CovariatesModuleCache::compute_similarity_obs(int obs_idx, int cls_idx) const {

    ClusterStats base_stats;

    // Handle new cluster case
    if (cls_idx > -1 && cls_idx < data.get_K())
        base_stats = covariate_cache.get_cluster_stats(cls_idx);

    const double obs_age = covariates_data.ages(obs_idx);

    const double log_ml_without = compute_log_marginal_likelihood(base_stats);

    base_stats.n += 1;
    base_stats.sum += obs_age;
    base_stats.sumsq += obs_age * obs_age;

    const double log_ml_with = compute_log_marginal_likelihood(base_stats);
    return log_ml_with - log_ml_without;
}

Eigen::VectorXd CovariatesModuleCache::compute_similarity_obs(int obs_idx) const {

    const int num_clusters = data.get_K();
    // Compute log similarities for each cluster
    Eigen::VectorXd log_similarities(num_clusters);
    const double obs_age = covariates_data.ages(obs_idx);
    for (int k = 0; k < num_clusters; ++k) {
        ClusterStats stats = covariate_cache.get_cluster_stats(k);
        const double log_ml_without = compute_log_marginal_likelihood(stats);

        stats.n += 1;
        stats.sum += obs_age;
        stats.sumsq += obs_age * obs_age;

        const double log_ml_with = compute_log_marginal_likelihood(stats);

        log_similarities(k) = log_ml_with - log_ml_without;
    }

    return log_similarities;
}

double CovariatesModuleCache::compute_log_marginal_likelihood_NNIG(const ClusterStats &stats) const {
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

    double lgamma_nu_n_temp;
    if (stats.n <= lgamma_nu_n.size() - 1) {
        lgamma_nu_n_temp = lgamma_nu_n[stats.n];
    } else {
        lgamma_nu_n_temp = std::lgamma(covariates_data.nu + 0.5 * n_dbl);
    }

    // log g(S) = log Γ(ν + n/2) - log Γ(ν) - n/2 log(2π)
    //            - 1/2 log(1+nB) + ν log(S₀) - (ν + n/2) log(S_n)
    return lgamma_nu_n_temp - lgamma_nu + n_dbl * const_term - 0.5 * std::log(one_plus_nB) + nu_logS0 -
           nu_n * std::log(S_n);
}

double CovariatesModuleCache::compute_log_marginal_likelihood_NN(const ClusterStats &stats) const {

    if (stats.n == 0) {
        return 0.0;
    }

    const double n_dbl = static_cast<double>(stats.n);
    const double inv_n = 1.0 / n_dbl;
    const double xbar = stats.sum * inv_n;

    // Centered sum of squares: SS = Σ(x_i - x̄)²
    const double ss = stats.sumsq - n_dbl * xbar * xbar;

    // v + n_j B term from integrating μ_j ~ N(m, B)
    const double v_plus_nB = covariates_data.v + n_dbl * covariates_data.B;

    // Mean deviation from prior mean
    const double dev = xbar - covariates_data.m;

    // Log of posterior variance: log(τ_j) = log(B) + log(v) - log(v + n_j B)
    double log_v_plus_nB_temp;
    if (stats.n <= log_v_plus_nB.size() - 1) {
        log_v_plus_nB_temp = log_v_plus_nB[stats.n];
    } else {
        log_v_plus_nB_temp = std::log(covariates_data.v + n_dbl * covariates_data.B);
    }

    const double log_tau_j = log_B + log_v - log_v_plus_nB_temp;

    // Compute log marginal likelihood using posterior variance form:
    // log q(x_j) = -n_j/2 log(2π) - n_j/2 log(v) - 1/2 log(B) + 1/2 log(τ_j)
    //              - SS/(2v) - n_j(x̄_j - m)² / (2(v + n_j B))
    return n_dbl * const_term - 0.5 * n_dbl * log_v - 0.5 * log_B + 0.5 * log_tau_j -
           0.5 * (ss / covariates_data.v + n_dbl * dev * dev / v_plus_nB);
}