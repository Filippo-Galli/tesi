#pragma once

#include "../../utils/ClusterInfo.hpp"

// Important to declare the struct useful.
// HINT: If the name changed modified the file which uses it accordingly.
struct ClusterStats {
    int n = 0;
    double sum = 0.0;
    double sumsq = 0.0;
};

class Covariate_cache : public ClusterInfo {

private:
    std::vector<ClusterStats> cluster_stats;

public:
    const Eigen::VectorXd continuos_covariates;
    
    Covariate_cache(Eigen::VectorXi &initial_allocations,
                    const Eigen::VectorXd &continuos_covariates)
        : continuos_covariates(continuos_covariates) {
        cluster_stats.clear();

        const int K = initial_allocations.maxCoeff() + 1;
        cluster_stats.resize(K);

        // Initialize cluster statistics based on initial allocations
        for (int i = 0; i < initial_allocations.size(); ++i) {
            int cluster = initial_allocations(i);
            if (cluster < 0)
                continue;

            double value = continuos_covariates(i);
            ClusterStats &stats = cluster_stats[cluster];
            stats.n++;
            stats.sum += value;
            stats.sumsq += value * value;
        }
    }

    /**
     * @brief Assigns a point to a cluster
     * @param index Index of the point to reassign
     * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
     * @param old_cluster Previous cluster index of the point
     * @throws std::out_of_range if index or cluster is invalid
     */
    void set_allocation(int index, int cluster, int old_cluster) override;

    /**
     * @brief Get cluster statistics for a specific cluster
     * @param cluster Index of the cluster
     * @return ClusterStats struct
     */
    inline ClusterStats get_cluster_stats(int cluster) const { return cluster_stats[cluster]; }

    /**
     * @brief Get cluster statistics reference for a specific cluster
     * @param cluster Index of the cluster
     * @return Const reference to ClusterStats struct
     */
    inline const ClusterStats &get_cluster_stats_ref(int cluster) const { return cluster_stats[cluster]; }

    /**
     * @brief Recomputes all cluster information from current allocations
     * @param K Current number of clusters
     */
    void recompute(const int K, const Eigen::VectorXi &allocations) override;

    /**
     * @brief Moves cluster information from one cluster to another
     * @param from_cluster Index of the source cluster
     * @param to_cluster Index of the target cluster
     */
    inline void move_cluster_info(int from_cluster, int to_cluster) override {
        cluster_stats[to_cluster] = std::move(cluster_stats[from_cluster]);
    };

    /**
     * @brief Removes information related to a specific cluster
     * @param cluster Index of the cluster to remove
     */
    void remove_info(int cluster) override;
};