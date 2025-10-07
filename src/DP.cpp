#include "DP.hpp"

double DP::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) {
    /**
     * @brief Computes the log prior probability of assigning a data point to an existing cluster.
     * @param cls_idx The index of the cluster
     * @return The log prior probability of assigning the data point to its current cluster.
    */

    int cluster_size = data.get_cluster_size(cls_idx);   
    return log(cluster_size);
}

double DP::gibbs_prior_new_cluster() {
    /**
     * @brief Computes the log prior probability of assigning a data point to a new cluster.
     * @return The log prior probability of assigning the data point to a new cluster.
    */
    return log(params.a);
}

double DP::prior_ratio_split(int ci, int cj) {
    /**
     * @brief Computes the prior ratio for a split operation in a split-merge MCMC algorithm.
     * @param ci the first cluster index involved in the split.
     * @param cj the second cluster index involved in the split.
     * @return The log prior ratio for the split operation.
    */

    int n_ci = data.get_cluster_size(ci);
    int n_cj = data.get_cluster_size(cj);

    return log(params.a) - lgamma(n_ci + n_cj) + lgamma(n_ci) + lgamma(n_cj);
}

double DP::prior_ratio_merge(int size_old_ci, int size_old_cj) {
    /**
     * @brief Computes the prior ratio for a merge operation in a split-merge MCMC algorithm.
     * @param size_old_ci the size of the first cluster before the merge.
     * @param size_old_cj the size of the second cluster before the merge.
     * @return The log prior ratio for the merge operation.
    */

    int size_merge = size_old_ci + size_old_cj;
    
    double log_acceptance_ratio = -log(params.a);
    log_acceptance_ratio += lgamma(size_merge);
    log_acceptance_ratio -= lgamma(size_old_ci);
    log_acceptance_ratio -= lgamma(size_old_cj);

    return log_acceptance_ratio;
}

double DP::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) {
    /**
     * @brief Computes the prior ratio for a shuffle operation in a split-merge MCMC algorithm.
     * @param size_old_ci the size of the first cluster before the shuffle.
     * @param size_old_cj the size of the second cluster before the shuffle.
     * @param ci the first cluster index involved in the shuffle.
     * @param cj the second cluster index involved in the shuffle.
     * @return The log prior ratio for the shuffle operation.
    */

    int n_ci = data.get_cluster_size(ci);
    int n_cj = data.get_cluster_size(cj);

    double log_acceptance_ratio = 0.0;
    log_acceptance_ratio += lgamma(n_ci);
    log_acceptance_ratio += lgamma(n_cj);
    log_acceptance_ratio -= lgamma(size_old_ci);
    log_acceptance_ratio -= lgamma(size_old_cj);

    return log_acceptance_ratio;
}

