#include "DPW.hpp"
#include <math.h>

int DPW::get_neighbors(int obs_idx, int cls_idx, Eigen::VectorXi allocations) const {
    /**
     * @brief Returns the number of neighbors for a given observation based on the adjacency matrix W.
     * @param obs_idx The index of the observation.
     * @param cls_idx The index of the cluster to consider for neighbor counting. If -1, counts all neighbors.
     * @return The number of neighbors for the observation.
    */

    if (allocations.size() == 0) 
        allocations = data.get_allocations();
    

    int neighbors = (params.W.row(obs_idx).array() * (allocations.array() == cls_idx).cast<int>().array()).sum();

    return neighbors;
}

double DPW::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) {
    /**
     * @brief Computes the log prior probability of assigning a data point to an existing cluster.
     * @param cls_idx The index of the cluster
     * @param obs_idx The index of the observation. 
     * @return The log prior probability of assigning the data point to its current cluster.
    */

    int cluster_size = data.get_cluster_size(cls_idx);   
    return log(cluster_size) + params.coefficient * get_neighbors(obs_idx, cls_idx);
}

double DPW::gibbs_prior_new_cluster() {
    /**
     * @brief Computes the log prior probability of assigning a data point to a new cluster.
     * @return The log prior probability of assigning the data point to a new cluster.
    */
    return log(params.a);
}

double DPW::prior_ratio_split(int ci, int cj) {
    /**
     * @brief Computes the prior ratio for a split operation in a split-merge MCMC algorithm.
     * @param ci the first cluster index involved in the split.
     * @param cj the second cluster index involved in the split.
     * @return The log prior ratio for the split operation.
    */

    double log_acceptance_ratio = log(params.a);

    int n_ci = data.get_cluster_size(ci);
    int n_cj = data.get_cluster_size(cj);

    log_acceptance_ratio += (n_ci > 0) ? lgamma(n_ci) : 0;
    log_acceptance_ratio += (n_cj > 0) ? lgamma(n_cj) : 0;
    log_acceptance_ratio -= (n_ci + n_cj > 0) ? lgamma(n_ci + n_cj) : 0;

    log_acceptance_ratio += params.coefficient

    return log_acceptance_ratio;
}

double DPW::prior_ratio_merge(int size_old_ci, int size_old_cj) {
    /**
     * @brief Computes the prior ratio for a merge operation in a split-merge MCMC algorithm.
     * @param size_old_ci the size of the first cluster before the merge.
     * @param size_old_cj the size of the second cluster before the merge.
     * @return The log prior ratio for the merge operation.
    */

    double log_acceptance_ratio = -log(params.a);
    int size_merge = size_old_ci + size_old_cj;
    log_acceptance_ratio += (size_merge > 0) ? lgamma(size_merge) : 0;
    log_acceptance_ratio -= (size_old_ci > 0) ? lgamma(size_old_ci) : 0;
    log_acceptance_ratio -= (size_old_cj > 0) ? lgamma(size_old_cj) : 0;

    return log_acceptance_ratio;
}

double DPW::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) {
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
    log_acceptance_ratio += (n_ci > 0) ? lgamma(n_ci) : 0;
    log_acceptance_ratio += (n_cj > 0) ? lgamma(n_cj) : 0;
    log_acceptance_ratio -= (size_old_ci > 0) ? lgamma(size_old_ci) : 0;
    log_acceptance_ratio -= (size_old_cj > 0) ? lgamma(size_old_cj) : 0;

    return log_acceptance_ratio;
}

