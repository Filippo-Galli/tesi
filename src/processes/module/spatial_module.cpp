#include "spatial_module.hpp"

int SpatialModule::get_neighbors_obs(int obs_idx, int cls_idx) const {
  /**
   * @brief Returns the number of neighbors for a given observation in a
   * specific cluster.
   *
   * This method counts the neighbors of an observation based on the adjacency
   * matrix W, considering only neighbors that belong to the specified cluster.
   * @param obs_idx The index of the observation.
   * @param cls_idx The index of the cluster to consider for neighbor counting.
   * @return The number of neighbors for the observation in the specified
   * cluster.
   */

  int neighbors = 0;
  Eigen::RowVectorXi row = params_module.W.row(obs_idx);
  for (int i = 0; i < row.size(); ++i) {
    int cluster_i = data_module.get_cluster_assignment(i);
    if (row(i) == 1 && cluster_i != -1 && cluster_i == cls_idx) {
      neighbors += 1;
    }
  }

  return neighbors;
}

int SpatialModule::get_neighbors_cls(int cls_idx, bool old_allo) const {
  /**
   * @brief Returns the total number of neighbors for all observations in a
   * given cluster.
   *
   * This method computes the sum of all neighbor connections within a cluster,
   * which is used in the spatial component of the prior calculations.
   * @param cls_idx The index of the cluster.
   * @param old_allo If true, uses the old allocations for neighbor counting;
   * otherwise, uses current allocations.
   * @return The total number of neighbors for the cluster.
   */

  const Eigen::VectorXi &allocations_reference =
      (old_allo && old_allocations_provider) ? old_allocations_provider()
                                             : data_module.get_allocations();
  Eigen::VectorXi obs_in_cluster = (allocations_reference.array() == cls_idx).cast<int>();

  const int total_neighbors = obs_in_cluster.dot(params_module.W * obs_in_cluster) / 2; // Divide by 2 to avoid double counting
  return total_neighbors;
}

Eigen::VectorXi SpatialModule::get_neighbors_obs(int obs_idx) const {
  /**
   * @brief Returns the number of neighbors for a given observation
   * regardless of cluster membership.
   * This method counts the total number of neighbors of an observation
   * based on the adjacency matrix W.
   * @param obs_idx The index of the observation.
   * @return The total number of neighbors for the observation for all clusters.
   */

  Eigen::VectorXi cluster_adjacency =
      Eigen::VectorXi::Zero(data_module.get_K());

  // Extract the adjacency edges of index
  Eigen::RowVectorXi row = params_module.W.row(obs_idx);
  for (int i = 0; i < row.size(); ++i) {
    int cluster_i = data_module.get_cluster_assignment(i);
    if (row(i) == 1 && cluster_i != -1) {
      cluster_adjacency(cluster_i) += 1;
    }
  }

  return cluster_adjacency;
}