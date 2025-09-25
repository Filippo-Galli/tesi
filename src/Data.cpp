#include "Data.hpp"
#include "Eigen/src/Core/Matrix.h"
#include <iostream>

Data::Data(const Eigen::MatrixXd &distances,
           const Eigen::VectorXi &initial_allocations)
    : D(distances), allocations(initial_allocations) {
  /**
   * @brief Constructor initializes the Data object with a distance matrix.
   * @details It also initializes the allocations vector based on the
   * first_allocation parameter.
   * @details If first_allocation is "all-in-one", all points are allocated to a
   * single cluster.
   * @details If first_allocation is "sequential", each point is allocated to
   * its own cluster.
   * @details Otherwise, it throws an invalid_argument exception.
   * @param distances The distance matrix to initialize the Data object.
   * @param param The parameters for the model.
   * @param first_allocation The initial allocation strategy for the clusters.
   * @throws std::invalid_argument if the distance matrix is not square or if
   * the first_allocation is invalid.
   */


  if (distances.rows() != distances.cols()) {
    throw std::invalid_argument("Distance matrix must be square");
  }

  n = distances.rows();

  // Initialize K (number of clusters) based on initial allocations
  if (allocations.size() == 0) {
    std::cout << "[INFO] No initial allocations provided, starting with all "
                 "points in one cluster."
              << std::endl;
    // If no initial allocations provided, start with all points in one cluster
    allocations = Eigen::VectorXi::Zero(n);
  }

  // Find the maximum cluster index to determine K
  K = allocations.maxCoeff() + 1;

  // Allocate cluster mapping
  for (int i = 0; i < n; ++i) {
    int cluster_id = allocations(i);
    if (cluster_id >= 0) {
      cluster_members[cluster_id].push_back(i);
    }
  }
}

double Data::get_distance(int i, int j) const {
  /**
   * @brief Gets the distance between two points.
   * @details It checks if the indices are within bounds and returns the
   * distance.
   * @param i The index of the first point.
   * @param j The index of the second point.
   * @return The distance between the two points.
   * @throws std::out_of_range if the indices are out of bounds.
   */

  #if VERBOSITY_LEVEL >= 1
  if (i < 0 || i >= D.rows() || j < 0 || j >= D.cols()) {
    throw std::out_of_range("Index out of bounds in get_distance");
  }
  #endif

  return D(i, j);
}

Eigen::VectorXi Data::get_cluster_assignments(int cluster) const {
  /**
   * @brief Gets the assignments of points to a specific cluster.
   * @details It returns a vector containing the indices of points assigned to
   * the specified cluster.
   * @param cluster The index of the cluster (0 to K-1).
   * @return A vector of indices of points assigned to the specified cluster.
   * @throws std::out_of_range if the cluster index is out of bounds.
   */

  #if VERBOSITY_LEVEL >= 1
  if (cluster > K || cluster < 0) {
    throw std::out_of_range("Index out of bounds in get_cluster_assignments");
  }
  #endif

  // Check if the cluster exists in the map before accessing it
  if (cluster_members.find(cluster) == cluster_members.end()) {
    return Eigen::VectorXi::Zero( 0); // Return empty vector for non-existent clusters
  }

  return Eigen::Map<const Eigen::VectorXi>(cluster_members.at(cluster).data(),
                                           cluster_members.at(cluster).size());
}

void Data::compact_cluster(int old_cluster) {
  /**
   * @brief Compacts the clusters by removing an empty cluster and shifting
   * allocations.
   * @details It removes the specified old_cluster and shifts allocations of
   * points in the last cluster to the old_cluster. It also updates the number
   * of clusters K.
   * @param old_cluster The index of the cluster to be removed (0 to K-1).
   * @throws std::out_of_range if the old_cluster index is out of bounds.
   */

  // Shift allocations of the last cluster to the old cluster
  for (int i = 0; i < n; ++i) {
    if (allocations(i) == K - 1) {
      allocations(i) = old_cluster;
    }
  }

  // Move members of the last cluster to the old cluster
  cluster_members.at(old_cluster) = std::move(cluster_members[K - 1]);

  cluster_members.erase(K - 1);
  K--; // Decrease the number of clusters
}

void Data::set_allocation(int index, int cluster) {
  /**
   * @brief Sets the allocation of a point to a specific cluster.
   * @details It updates the allocations vector and calls update_cluster_sizes
   * to refresh cluster sizes.
   * @param index The index of the point to be allocated.
   * @param cluster The cluster index to which the point should be allocated (0
   * to K-1 for existing, K for new).
   */

  // Bounds checking for index
  #if VERBOSITY_LEVEL >= 1
  if (index < 0 || index >= n) {
    throw std::out_of_range("Index out of bounds in set_allocation");
  }
  #endif

  int old_cluster = allocations(index);
  int final_cluster;

  // Invalidating index: cluster = -1
  if (cluster == -1) {
    final_cluster = -1; // Mark as unallocated
    allocations(index) = final_cluster;

    // Remove index from old cluster members (check if old_cluster exists)
    if (old_cluster != -1 &&
        cluster_members.find(old_cluster) != cluster_members.end()) {
      // Remove index from old cluster members
      cluster_members[old_cluster].erase(
          std::remove(cluster_members[old_cluster].begin(),
                      cluster_members[old_cluster].end(), index),
          cluster_members[old_cluster].end());

      // If the old cluster becomes empty, we need to remove it
      if (cluster_members.find(old_cluster) != cluster_members.end() && cluster_members.at(old_cluster).size() == 0)
        compact_cluster(old_cluster);
      else if (cluster_members.find(old_cluster) == cluster_members.end())
        throw std::runtime_error("Inconsistent state: old_cluster not found in "
                                 "cluster_members during deallocation");
    }

    return;
  }

  // label an unallocated index
  if (old_cluster == -1) {

    // Setting to an existing cluster
    if (cluster >= 0 && cluster < K) {
      final_cluster = cluster;
      allocations(index) = final_cluster;
      // Add index to new cluster members
      cluster_members[final_cluster].push_back(index);
    }
    // Setting to a new cluster
    else if (cluster == K) {
      final_cluster = K;
      allocations(index) = final_cluster; // Assign to new cluster
      K++;                                // Increase the number of clusters
      // Add index to new cluster members
      cluster_members[final_cluster].push_back(index);
    } else {
      throw std::out_of_range( "Invalid cluster index in set_allocation for unallocated point");
    }
    return;
  }

  // Remove index from old cluster members (check if old_cluster exists)
  if (old_cluster != -1 &&
      cluster_members.find(old_cluster) != cluster_members.end()) {
    // Remove index from old cluster members
    cluster_members[old_cluster].erase(
        std::remove(cluster_members[old_cluster].begin(),
                    cluster_members[old_cluster].end(), index),
        cluster_members[old_cluster].end());
  }

  // If the point is already allocated to a cluster
  if (cluster < K) {
    final_cluster = cluster;
    allocations(index) = final_cluster; // Update allocation
    cluster_members[final_cluster].push_back(index);

    // If the old cluster becomes empty, we need to remove it
    if (old_cluster != -1 &&
        cluster_members.find(old_cluster) != cluster_members.end() &&
        cluster_members.at(old_cluster).size() == 0)
      compact_cluster(old_cluster);

  } else if (cluster == K) {
    final_cluster = K;
    allocations(index) = final_cluster; // Assign to new cluster
    K++;                                // Increase the number of clusters
    cluster_members[final_cluster].push_back(index);
  } else {
    throw std::out_of_range(
        "Invalid cluster index in set_allocation for allocated point");
  }
}

void Data::set_allocations(const Eigen::VectorXi &new_allocations) {
  /**
   * @brief Sets the allocations vector to a new vector and updates cluster
   * sizes.
   * @details It checks if the new allocations vector is of the correct size
   * and updates the allocations and cluster sizes accordingly.
   * @param new_allocations The new allocations vector to set.
   * @throws std::invalid_argument if the new allocations vector is not of the
   * correct size.
   */

  #if VERBOSITY_LEVEL >= 1
  if (new_allocations.size() != n) {
    throw std::invalid_argument(
        "New allocations vector must be of size n, set_allocations failed");
  }
  #endif

  allocations = new_allocations;

  // Update cluster members mapping
  cluster_members.clear();
  for (int i = 0; i < n; ++i)
    cluster_members[allocations(i)].push_back(i);

  // Update K based on the new allocations
  K = allocations.maxCoeff() + 1;
}