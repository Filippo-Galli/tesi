#include "spatial_module.hpp"

int SpatialModule::get_neighbors_obs(int obs_idx, int cls_idx) const {
  int neighbors = 0;
  const std::vector<int> &row = neighbor_cache[obs_idx];

  // Iterate through cached neighbors and count those in the specified cluster
  for (size_t i = 0; i < row.size(); ++i) {
    int cluster_i = data_module.get_cluster_assignment(row[i]);
    if (cluster_i != -1 && cluster_i == cls_idx) {
      neighbors += 1;
    }
  }

  return neighbors;
}

void SpatialModule::neighbor_cache_compute() {
  const int N = data_module.get_n();
  neighbor_cache.resize(N);

  // For each observation, store indices of its neighbors (where W(i,j) == 1)
  for (int obs_idx = 0; obs_idx < N; ++obs_idx) {
    Eigen::RowVectorXi row = params_module.W.row(obs_idx);
    neighbor_cache[obs_idx].reserve(row.sum());

    for (int j = 0; j < row.size(); ++j) {
      if (row(j) == 1) {
        neighbor_cache[obs_idx].push_back(j);
      }
    }
  }
}

int SpatialModule::get_neighbors_cls(int cls_idx, bool old_allo) const {
  const Eigen::VectorXi &allocations = (old_allo && old_allocations_provider) ? old_allocations_provider() : data_module.get_allocations();
  
  int total_neighbors = 0;
  const int N = data_module.get_n();
  
  // Only iterate over points in this cluster
  for (int i = 0; i < N; ++i) {
    if (allocations(i) != cls_idx) continue;
    
    // Count neighbors also in this cluster
    const std::vector<int>& neighbors = neighbor_cache[i];
    for (int j : neighbors) {
      if (allocations(j) == cls_idx && j > i) {  // j > i avoids double counting
        ++total_neighbors;
      }
    }
  }
  
  return total_neighbors;
}

Eigen::VectorXi SpatialModule::get_neighbors_obs(int obs_idx) const {
  Eigen::VectorXi cluster_adjacency = Eigen::VectorXi::Zero(data_module.get_K());

  // Use cached neighbor indices instead of iterating over full adjacency matrix
  const std::vector<int> &row = neighbor_cache[obs_idx];

  // Count neighbors in each cluster
  for (size_t i = 0; i < row.size(); ++i) {
    int neighbor_idx = row[i];
    int cluster_i = data_module.get_cluster_assignment(neighbor_idx);

    if (cluster_i != -1) {
      cluster_adjacency(cluster_i) += 1;
    }
  }

  return cluster_adjacency;
}