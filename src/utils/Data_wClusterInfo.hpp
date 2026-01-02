#pragma once

#include "Data.hpp"
#include "ClusterInfo.hpp"

class Data_wClusterInfo : public Data {

protected:
    ClusterInfo &cluster_info;

    void compact_cluster(int old_cluster);

public:
    Data_wClusterInfo(const Params &p, ClusterInfo &ci, const Eigen::VectorXi &initial_allocations = Eigen::VectorXi())
        : Data(p, initial_allocations), cluster_info(ci) {}

    /**
     * @brief Assigns a point to a cluster
     * @param index Index of the point to reassign
     * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
     * @throws std::out_of_range if index or cluster is invalid
     */
    void set_allocation(int index, int cluster);

    /**
     * @brief Sets all cluster allocations at once
     * @param new_allocations Vector of cluster assignments for all points
     * @throws std::invalid_argument if vector size doesn't match number of points
     */
    void set_allocations(const Eigen::VectorXi &new_allocations);

    /**
     * @brief Restores allocations, cluster memberships, and cluster count from a saved state
     *
     * @param old_allocations Vector holding the saved allocations (swapped into place)
     * @param old_cluster_members Map of saved cluster memberships (swapped into place)
     * @param old_K Number of clusters in the saved state
     */
    void restore_state(Eigen::VectorXi &old_allocations, std::unordered_map<int, std::vector<int>> &old_cluster_members,
                       int old_K) override;

    virtual ~Data_wClusterInfo() = default;
};