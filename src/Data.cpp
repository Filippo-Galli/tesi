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
    K = 1;
  } else {
    // Find the maximum cluster index to determine K
    K = allocations.maxCoeff() + 1;
  }

  update_cluster_sizes(); // initialize cluster sizes
}

void Data::update_cluster_sizes() {
  /**
   * @brief Updates the sizes of each cluster based on the current allocations.
   * @details It counts the number of points in each cluster and updates the
   * cluster_sizes vector.
   */
  cluster_sizes = Eigen::VectorXi::Zero(K);
  for (int i = 0; i < n; ++i) {
    cluster_sizes(allocations(i))++;
  }
}

void Data::update_cluster_sizes(unsigned cluster) {
  /**
   * @brief Updates the size of a specific cluster.
   * @details It counts the number of points in the specified cluster and
   * updates the cluster_sizes
   * @param cluster The index of the cluster to update (0 to K-1).
   * @throws std::out_of_range if the cluster index is out of bounds.
   */

  double count = 0;
  for (int i = 0; i < n; ++i) {
    if (allocations(i) == cluster)
      count++;
  }
  cluster_sizes(cluster) = count;
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
  if (i < 0 || i >= D.rows() || j < 0 || j >= D.cols()) {
    throw std::out_of_range("Index out of bounds");
  }
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

  Eigen::VectorXi assignments = Eigen::VectorXi::Zero(cluster_sizes(cluster));
  int count = 0;

  for (int i = 0; i < n; ++i) {
    if (allocations(i) == cluster) {
      assignments(count++) = i;
    }
  }

  return assignments;
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

  int old_cluster = allocations(index);

  int final_cluster;

  // Invalidating the index : -1
  if (cluster == -1) {
    final_cluster = -1; // Mark as unallocated
    cluster_sizes(old_cluster)--;

    // If the old cluster becomes empty, we need to remove it
    if (cluster_sizes(old_cluster) == 0) {
      // Shift allocations of the last cluster to the old cluster
      for (int i = 0; i < n; ++i) {
        if (allocations(i) == K - 1) {
          allocations(i) = old_cluster;
        }
      }

      cluster_sizes(old_cluster) =
          cluster_sizes(K - 1);            // Update the size of the old cluster
      K--;                                 // Decrease the number of clusters
      cluster_sizes.conservativeResize(K); // Resize the cluster sizes vector
    }
    allocations(index) = final_cluster;
    return;
  }

  // label an unallocated index
  if (old_cluster == -1) {

    // Setting to an existing cluster
    if (cluster >= 0 && cluster < K) {
      final_cluster = cluster;
      allocations(index) = final_cluster;
      cluster_sizes(final_cluster)++;
    }
    // Setting to a new cluster
    else if (cluster == K) {
      final_cluster = K;
      allocations(index) = final_cluster;  // Assign to new cluster
      K++;                                 // Increase the number of clusters
      cluster_sizes.conservativeResize(K); // Resize the cluster sizes vector
      cluster_sizes(final_cluster) = 1;    // Initialize size of new cluster
    } else {
      throw std::out_of_range("Invalid cluster index");
    }
    return;
  }

  // If the point is already allocated to a cluster
  if (cluster < K) {
    final_cluster = cluster;
    allocations(index) = final_cluster; // Update allocation
    cluster_sizes(old_cluster)--;       // Decrease size of old cluster
    cluster_sizes(final_cluster)++;     // Increase size of new cluster

    // If the old cluster becomes empty, we need to remove it
    if (cluster_sizes(old_cluster) == 0) {
      // Shift allocations of the last cluster to the old cluster
      for (int i = 0; i < n; ++i) {
        if (allocations(i) == K - 1) {
          allocations(i) = old_cluster;
        }
      }
      cluster_sizes(old_cluster) =
          cluster_sizes(K - 1);            // Update the size of the old cluster
      K--;                                 // Decrease the number of clusters
      cluster_sizes.conservativeResize(K); // Resize the cluster sizes vector
    }
  } else if (cluster == K) {
    final_cluster = K;
    allocations(index) = final_cluster;  // Assign to new cluster
    K++;                                 // Increase the number of clusters
    cluster_sizes.conservativeResize(K); // Resize the cluster sizes vector
    cluster_sizes(final_cluster) = 1;    // Initialize size of new cluster
  } else {
    throw std::out_of_range("Invalid cluster index");
  }
}

void Data::print() const {
  /**
   * @brief Prints the contents of the Data object.
   * @details It prints the distances matrix, number of points, allocations
   * vector, and number of clusters.
   */

  // std::cout << "Distances Matrix:\n" << D << std::endl;
  // std::cout << "Number of points: " << n << std::endl;

  std::cout << std::endl;
  std::cout << "Allocations:\n" << allocations.transpose() << std::endl;
  std::cout << "Number of clusters: " << K << std::endl;
  std::cout << "Cluster Sizes:\n" << cluster_sizes.transpose() << std::endl;
}