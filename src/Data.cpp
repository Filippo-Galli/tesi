#include "Data.hpp"
#include "Eigen/src/Core/Matrix.h"
#include <iostream>

Data::Data(const Eigen::MatrixXd& distances, const std::string & first_allocation) : D(distances) {
    /**
    * @brief Constructor initializes the Data object with a distance matrix.
    * @details It also initializes the allocations vector based on the first_allocation parameter.
    * @details If first_allocation is "all-in-one", all points are allocated to a single cluster.
    * @details If first_allocation is "sequential", each point is allocated to its own cluster.
    * @details Otherwise, it throws an invalid_argument exception.
    * @param distances The distance matrix to initialize the Data object.
    * @param param The parameters for the model.
    * @param first_allocation The initial allocation strategy for the clusters.
    * @throws std::invalid_argument if the distance matrix is not square or if the first_allocation is invalid.
    */

    if (distances.rows() != distances.cols()) {
        throw std::invalid_argument("Distance matrix must be square");
    }

    n = distances.rows();

    if(first_allocation == "all-in-one") {
        allocations = Eigen::VectorXi::Constant(n, 0); // all points in one cluster
        K = 1; // only one cluster
    } else if (first_allocation == "sequential") {
        allocations = Eigen::VectorXi::LinSpaced(n, 0, n - 1); // sequential allocations
        K = n; // each point in its own cluster
    } else {
        throw std::invalid_argument("Invalid first allocation strategy");
    }

    update_cluster_sizes(); // initialize cluster sizes
}

void Data::update_cluster_sizes() {
    /**
    * @brief Updates the sizes of each cluster based on the current allocations.
    * @details It counts the number of points in each cluster and updates the cluster_sizes vector.
    */
    cluster_sizes = Eigen::VectorXi::Zero(K);
    for (int i = 0; i < n; ++i) {
        cluster_sizes(allocations(i))++;
    }
}

void Data::update_cluster_sizes(unsigned cluster){
    double count = 0;
    for (int i = 0; i < n; ++i) {
        if(allocations(i) == cluster)
            count++;
    }
    cluster_sizes(cluster) = count;
}

double Data::get_distance(int i, int j) const {
    /**
    * @brief Gets the distance between two points.
    * @details It checks if the indices are within bounds and returns the distance.
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
    * @details It returns a vector containing the indices of points assigned to the specified cluster.
    * @param cluster The index of the cluster (0 to K-1).
    * @return A vector of indices of points assigned to the specified cluster.
    * @throws std::out_of_range if the cluster index is out of bounds.
    */
    
    if (cluster < 0 || cluster >= K) {
        throw std::out_of_range("Cluster index out of bounds");
    }

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
    * @details It updates the allocations vector and calls update_cluster_sizes to refresh cluster sizes.
    * @param index The index of the point to be allocated.
    * @param cluster The cluster index to which the point should be allocated (0 to K-1 for existing, K for new).
    */
    
    int old_cluster = allocations(index);
    
    int final_cluster;
    
    if (cluster < K && get_cluster_size(old_cluster) > 1) {
        // Moving to existing cluster, old cluster remains non-empty
        cluster_sizes(cluster)++;
        cluster_sizes(old_cluster)--;
        final_cluster = cluster;
    } 
    else if (cluster < K && get_cluster_size(old_cluster) == 1) {

        // Relabel the last cluster as the old cluster
        int count = 0;
        for(int i = 0; i < n; ++i) {
            if (allocations(i) == K - 1 ) {
                allocations(i) = old_cluster; // Relabel the point
                count++;
            }
        }

        if(cluster == K - 1) {
            final_cluster = old_cluster; // No change, just relabel
            cluster_sizes(old_cluster) = count +1; // cluster moved + the final assignment
        } else {
            final_cluster = cluster; // Move to the new cluster
            cluster_sizes(old_cluster) = count; // Old cluster loses one member
            cluster_sizes(cluster) += 1; // New cluster gains one member
        }
        
        // Reduce number of clusters
        K--;
        cluster_sizes.conservativeResize(K);
        
    } 
    else if (cluster == K && get_cluster_size(old_cluster) > 1) {
        // Creating new cluster, old cluster remains non-empty
        cluster_sizes.conservativeResize(K + 1);
        cluster_sizes(K) = 1;  // New cluster gets 1 member
        cluster_sizes(old_cluster)--;  // Old cluster loses 1 member
        final_cluster = K;  // The new cluster index
        K++;  // Increment K after setting final_cluster
    } 
    else if (cluster == K && get_cluster_size(old_cluster) == 1) {
        // Extending and then shrinking the cluster is uselss
        final_cluster = old_cluster;  // No change, just relabel
    } 
    else {
        throw std::invalid_argument("Invalid cluster transition: cluster=" + std::to_string(cluster) + 
                                ", K=" + std::to_string(K) + ", old_cluster=" + std::to_string(old_cluster));
    }
    
    allocations(index) = final_cluster;
}

void Data::print() const {
    /** 
    * @brief Prints the contents of the Data object.
    * @details It prints the distances matrix, number of points, allocations vector, and number of clusters.
    */

    //std::cout << "Distances Matrix:\n" << D << std::endl;
    //std::cout << "Number of points: " << n << std::endl;

    std::cout << std::endl;
    std::cout << "Allocations:\n" << allocations.transpose() << std::endl;
    std::cout << "Number of clusters: " << K << std::endl;
    std::cout << "Cluster Sizes:\n" << cluster_sizes.transpose() << std::endl;
}