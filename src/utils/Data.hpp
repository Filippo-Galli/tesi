/**
 * @file Data.hpp
 * @brief Data structure for managing point distances and cluster allocations
 */

#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_map>

// compilation time verbosity level
#ifndef VERBOSITY_LEVEL
#define VERBOSITY_LEVEL 0
#endif

/**
 * @class Data
 * @brief Manages distance matrices and cluster allocations for points
 *
 * This class stores a distance matrix between points and tracks their
 * cluster assignments. It provides methods to query distances, manage
 * cluster allocations, and maintain cluster membership information.
 */
class Data {
private:
  Eigen::MatrixXd D; ///< Distance matrix between all points
  int n;             ///< Number of points

  Eigen::VectorXi allocations; ///< Cluster allocation for each point
  int K;                       ///< Current number of clusters

  /// Maps cluster indices to vectors of point indices in that cluster
  std::unordered_map<int, std::vector<int>> cluster_members;

  /**
   * @brief Removes an empty cluster and compacts cluster indices
   * @param old_cluster Index of the cluster to remove
   */
  void compact_cluster(int old_cluster);

public:
  /**
   * @brief Constructs a Data object with a distance matrix
   * @param distances Square matrix of pairwise distances between points
   * @param initial_allocations Optional initial cluster assignments (default:
   * all in one cluster)
   * @throws std::invalid_argument if distance matrix is not square
   */
  Data(const Eigen::MatrixXd &distances,
       const Eigen::VectorXi &initial_allocations = Eigen::VectorXi());

  // Getters

  /**
   * @brief Gets the distance between two points
   * @param i Index of first point
   * @param j Index of second point
   * @return Distance between points i and j
   */
  double get_distance(int i, int j) const;

  /**
   * @brief Gets the total number of points
   * @return Number of points
   */
  int get_n() const { return n; }

  /**
   * @brief Gets the current number of clusters
   * @return Number of clusters
   */
  int get_K() const { return K; }

  /**
   * @brief Gets the cluster allocations vector
   * @return Reference to the allocations vector
   */
  const Eigen::VectorXi &get_allocations() const { return allocations; }

  /**
   * @brief Gets the size of a specific cluster
   * @param cluster_index Index of the cluster
   * @return Number of points in the cluster (0 if cluster doesn't exist)
   */
  int get_cluster_size(unsigned cluster_index) const {
    auto it = cluster_members.find(cluster_index);
    return (cluster_index < K && it != cluster_members.end())
               ? it->second.size()
               : 0;
  }

  /**
   * @brief Gets the cluster assignment of a specific point
   * @param index Index of the point
   * @return Cluster index the point is assigned to
   * @throws std::out_of_range if index is out of bounds
   */
  int get_cluster_assignment(int index) const {
    if (index < 0 || index >= n) {
      throw std::out_of_range("Index out of bounds in get_cluster_assignment");
    }
    return allocations(index);
  }

  /**
   * @brief Gets all point indices assigned to a specific cluster
   * @param cluster Index of the cluster
   * @return Vector of point indices in the cluster
   * @throws std::out_of_range if cluster index is invalid
   */
  Eigen::VectorXi get_cluster_assignments(int cluster) const;

  /**
   * @brief Gets all point indices assigned to a specific cluster (map form)
   * @param cluster Index of the cluster
   * @return Map to vector of point indices in the cluster
   * @throws std::out_of_range if cluster index is invalid
   */
  Eigen::Map<const Eigen::VectorXi> get_cluster_assignments_ref(int cluster) const;

  // Setters

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
};