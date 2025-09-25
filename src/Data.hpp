#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <unordered_map>

// compilation time verbosity level
#ifndef VERBOSITY_LEVEL
#define VERBOSITY_LEVEL 0
#endif

class Data {
private:
  Eigen::MatrixXd D; // distances matrix
  int n;             // number of points

  Eigen::VectorXi allocations; // cluster allocations of each point
  int K;                       // number of clusters
  std::unordered_map<int, std::vector<int>>
      cluster_members; // members of each cluster

  void compact_cluster(int old_cluster);

public:
  Data(const Eigen::MatrixXd &distances,
       const Eigen::VectorXi &initial_allocations = Eigen::VectorXi());

  // Getters
  double get_distance(int i, int j) const;
  int get_n() const { return n; }
  int get_K() const { return K; }
  const Eigen::VectorXi &get_allocations() const { return allocations; }
  int get_cluster_size(unsigned cluster_index) const {
    return (cluster_index < K &&
            cluster_members.find(cluster_index) != cluster_members.end())
               ? cluster_members.at(cluster_index).size()
               : 0;
  }
  int get_cluster_assignment(int index) const {
    if (index < 0 || index >= n) {
      throw std::out_of_range("Index out of bounds in get_cluster_assignment");
    }
    return allocations(index);
  }

  Eigen::VectorXi get_cluster_assignments(int cluster) const;

  // Setters
  void set_allocation(int index, int cluster);
  void set_allocations(const Eigen::VectorXi &new_allocations);
};