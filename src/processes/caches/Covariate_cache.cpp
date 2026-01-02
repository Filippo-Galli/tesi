#include "Covariate_cache.hpp"
#include <sstream>
#include <stdexcept>

void Covariate_cache::set_allocation(int index, int cluster, int old_cluster) {
    // Remove point from old cluster stats
    if (old_cluster != -1) {
        if (old_cluster >= static_cast<int>(cluster_stats.size())) {
            std::stringstream ss;
            ss << "Covariate_cache::set_allocation: old_cluster (" << old_cluster
               << ") out of bounds (size=" << cluster_stats.size() << ")";
            throw std::out_of_range(ss.str());
        }
        ClusterStats &old_stats = cluster_stats[old_cluster];
        double value = covariates.ages(index);
        old_stats.n--;
        old_stats.sum -= value;
        old_stats.sumsq -= value * value;
    }

    // Add point to new cluster stats
    if (cluster != -1) {
        if (cluster < 0) {
            throw std::out_of_range("Covariate_cache::set_allocation: cluster index negative");
        }
        // Ensure cluster_stats vector is large enough
        if (cluster >= static_cast<int>(cluster_stats.size())) {
            cluster_stats.resize(cluster + 1);
        }

        ClusterStats &new_stats = cluster_stats[cluster];
        double value = covariates.ages(index);
        new_stats.n++;
        new_stats.sum += value;
        new_stats.sumsq += value * value;
    }
}

void Covariate_cache::recompute(const int K, const Eigen::VectorXi &allocations) {
    // Clear existing stats
    cluster_stats.clear();

    // Get the number of clusters from covariates allocations
    cluster_stats.resize(K);

    // Recompute stats from scratch
    for (int i = 0; i < allocations.size(); ++i) {
        int cluster = allocations(i);
        if (cluster < 0)
            continue; // Skip unallocated points

        double value = covariates.ages(i);

        ClusterStats &stats = cluster_stats[cluster];
        stats.n++;
        stats.sum += value;
        stats.sumsq += value * value;
    }
}

void Covariate_cache::move_cluster_info(int from_cluster, int to_cluster) {
    if (from_cluster >= static_cast<int>(cluster_stats.size()) ||
        to_cluster >= static_cast<int>(cluster_stats.size())) {
        std::stringstream ss;
        ss << "Covariate_cache::move_cluster_info: cluster index out of bounds. "
           << "from: " << from_cluster << ", to: " << to_cluster << ", size: " << cluster_stats.size();
        throw std::out_of_range(ss.str());
    }
    // move stats from from_cluster to to_cluster
    cluster_stats[to_cluster] = std::move(cluster_stats[from_cluster]);
}

void Covariate_cache::remove_info(int cluster) {
    // Remove stats for the specified cluster
    if (cluster >= 0 && cluster < static_cast<int>(cluster_stats.size())) {
        if (cluster == static_cast<int>(cluster_stats.size()) - 1) {
            cluster_stats.pop_back();
        } else {
            cluster_stats.erase(cluster_stats.begin() + cluster);
        }
    } else {
        throw std::out_of_range("Covariate_cache::remove_info: cluster index out of bounds");
    }
}