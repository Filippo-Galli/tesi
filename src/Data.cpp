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

  #if VERBOSITY_LEVEL >= 1
  if (old_cluster < 0 || old_cluster >= K) {
    throw std::out_of_range("old_cluster index out of bounds in compact_cluster");
  }
  #endif

  // Early exit if there's only one cluster or if compacting the last cluster
  if (K <= 1 || old_cluster == K - 1) {
    // Just remove the cluster and decrease K
    cluster_members.erase(old_cluster);
    K--;
    return;
  }

  auto last_cluster_it = cluster_members.find(K - 1);
  bool last_cluster_exists = (last_cluster_it != cluster_members.end());

  // Shift allocations of the last cluster to the old cluster
  if (last_cluster_exists && !last_cluster_it->second.empty()) {
    for (int point_index : last_cluster_it->second) {
      allocations(point_index) = old_cluster;
    }
    
    // Move members of the last cluster to the old cluster using move semantics
    cluster_members[old_cluster] = std::move(last_cluster_it->second);
  } else {
    // Last cluster is empty or doesn't exist, just clear the old cluster
    cluster_members[old_cluster].clear();
  }

  // Remove the last cluster
  cluster_members.erase(K - 1);
  K--; // Decrease the number of clusters
}

void Data::set_allocation(int index, int cluster) {
  /**
   * @brief Sets the allocation of a point to a specific cluster.
   * @details It updates the allocations vector and calls compact_cluster
   * to remove empty clusters when necessary.
   * @param index The index of the point to be allocated.
   * @param cluster The cluster index to which the point should be allocated (0
   * to K-1 for existing, K for new, -1 for unallocated).
   */

  // Bounds checking for index
  #if VERBOSITY_LEVEL >= 1
  if (index < 0 || index >= n) {
    throw std::out_of_range("Index out of bounds in set_allocation");
  }
  #endif

  // Validate cluster parameter early
  #if VERBOSITY_LEVEL >= 1
  if (cluster < -1 || cluster > K) {
    throw std::out_of_range("Invalid cluster index in set_allocation");
  }
  #endif

  int old_cluster = allocations(index);
  
  // Early exit if no change needed
  if (old_cluster == cluster) {
    return;
  }

  auto old_cluster_it = cluster_members.end();
  bool old_cluster_exists = false;

  // Check if old_cluster exists in the map
  if (old_cluster != -1) {
    old_cluster_it = cluster_members.find(old_cluster);
    old_cluster_exists = (old_cluster_it != cluster_members.end());
    
    #if VERBOSITY_LEVEL >= 1
    if (!old_cluster_exists) {
      throw std::runtime_error("Inconsistent state: old_cluster not found in cluster_members");
    }
    #endif
  }

  // Remove from old cluster first (if applicable)
  if (old_cluster_exists) {
    auto& old_members = old_cluster_it->second;
    // Use find + erase instead of remove + erase for single element (more efficient)
    auto it = std::find(old_members.begin(), old_members.end(), index);
    if (it != old_members.end()) {
      old_members.erase(it);
    }
  }

  // Update allocation
  allocations(index) = cluster;

  // Handle new cluster assignment
  if (cluster == -1) {
    // Point becomes unallocated - nothing more to do after removal above
  } 
  else if (cluster == K) {
    // New cluster creation
    cluster_members[K].push_back(index);
    K++;
  } 
  else {
    // Existing cluster assignment - cache the target cluster iterator
    auto target_cluster_it = cluster_members.find(cluster);
    if (target_cluster_it != cluster_members.end()) {
      target_cluster_it->second.push_back(index);
    } 
    else {
      // Create the cluster if it doesn't exist
      cluster_members[cluster].push_back(index);
    }
  }

  // Check if old cluster became empty and needs compaction
  if (old_cluster_exists && old_cluster_it->second.empty()) {
    compact_cluster(old_cluster);
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