#include "DP_splitmerge_W_MartinezMena.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Rcpp/iostream/Rstreambuf.h"

void DPSplitMergeWMartinezMena::choose_indeces() {
  /**
   * @brief Randomly choose two distinct indices i and j from the data points.
   *        Identify their clusters ci and cj.
   *        Prepare the launch_state and S vectors for the split-merge
   * operation.
   */

  std::uniform_int_distribution<> dis(0, data.get_n() - 1);

  idx_i = dis(gen);
  do {
    idx_j = dis(gen);
  } while (idx_j == idx_i); // Ensure i and j are distinct

  // Get the clusters of the chosen indices
  ci = data.get_cluster_assignment(idx_i);
  cj = data.get_cluster_assignment(idx_j);

  // Get number of points in clusters ci and cj
  int size_ci = data.get_cluster_size(ci);
  int size_cj = data.get_cluster_size(cj);

  // Initialize launch_state with current allocations of points in clusters ci
  // and cj while S with their indices
  size_t launch_state_size = ci == cj ? size_ci - 2 : size_ci + size_cj - 2; // Exclude points i and j from the launch state

  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);
  original_allocations = data.get_allocations(); // Store original allocations in case of rejection

  // Properly collect all points from clusters ci and cj
  int s_idx = 0;
  for (int idx = 0; idx < data.get_n(); ++idx) {
    if (idx == idx_i || idx == idx_j)
      continue; // Skip points i and j

    int temp_cluster = data.get_cluster_assignment(idx);
    if (temp_cluster == ci || temp_cluster == cj) {
      S(s_idx) = idx;
      launch_state(s_idx) = temp_cluster;
      s_idx++;
    }
  }

// Ensure we collected the expected number of points
#if VERBOSITY_LEVEL >= 1
  if (s_idx != launch_state_size) { // since s_idx is zero-based
    // Rcpp::Rcout << "[ERROR] s_idx = " << s_idx << ", launch_state_size = " <<
    // launch_state_size << std::endl;
    throw std::runtime_error(
        "Mismatch in expected cluster sizes during split-merge initialization");
  }
#endif
}

int DPSplitMergeWMartinezMena::cluster_neighbors(int cluster) {
  /**
   * @brief Count the number of neighbors for a given cluster based on the
   * adjacency matrix W.
   * @param cluster The cluster index for which to count neighbors.
   * @return The number of neighbors for the specified cluster.
   */

  int count = 0;

  auto cluster_allocations = data.get_cluster_assignments(cluster);

  for(auto && idx : cluster_allocations) {
    count += (params.W.row(idx).array() * (data.get_allocations().array() != cluster).cast<int>().array()).sum();
  }
  return count;
}

void DPSplitMergeWMartinezMena::restricted_gibbs(int iterations, bool only_probabilities) {
  /**
   * @brief Perform restricted Gibbs sampling on the points in S to propose new
   * allocations for iter iterations.
   * @param iterations Number of Gibbs sampling iterations to perform.
   * @param only_probabilities If true, only compute the probabilities without
   * changing allocations (used in merge move).
   */

  for (int i = 0; i < iterations; ++i) {
    for (int idx = 0; idx < S.size(); ++idx) {
      int point_idx = S(idx);
      int current_cluster = launch_state(idx);

      // Remove point from its current cluster
      data.set_allocation(point_idx, -1); // Temporarily unassign the point

      // Compute probabilities for each cluster (ci and cj)
      Eigen::Vector2d log_probs;
      int c_i_size_minus_idx = data.get_cluster_size(ci);
      int c_j_size_minus_idx = data.get_cluster_size(cj);

      // "DP prior" term: cluster size (excluding the current point)
      log_probs(0) = likelihood.point_loglikelihood_cond(point_idx, ci);
      log_probs(1) = likelihood.point_loglikelihood_cond(point_idx, cj);
      log_probs(0) += (c_i_size_minus_idx > 0) ? log(c_i_size_minus_idx) : -std::numeric_limits<double>::infinity();
      log_probs(1) += (c_j_size_minus_idx > 0) ? log(c_j_size_minus_idx) : -std::numeric_limits<double>::infinity();

      // "W" term: adjacency information
      int neighbors_in_ci = (params.W.row(point_idx).array() * (data.get_allocations().array() == ci).cast<int>().array()).sum();
      int neighbors_in_cj = (params.W.row(point_idx).array() * (data.get_allocations().array() == cj).cast<int>().array()).sum();

      log_probs(0) += params.coefficient * neighbors_in_ci; // Add log of number of neighbors in ci if > 0
      log_probs(1) += params.coefficient * neighbors_in_cj; // Add log of number of neighbors in cj if > 0

      // Normalize to get probabilities
      double max_log_prob = log_probs.maxCoeff();
      Eigen::Vector2d probs = (log_probs.array() - max_log_prob).exp();
      probs /= probs.sum();

      if (!only_probabilities) {
        // Sample new cluster based on computed probabilities
        std::discrete_distribution<int> dist(probs.data(),
                                             probs.data() + probs.size());
        int new_cluster_idx = dist(gen);
        int new_cluster = (new_cluster_idx == 0) ? ci : cj;

        // Assign point to the new cluster
        data.set_allocation(point_idx, new_cluster);

        // if last iteration, accumulate the log probability of the move
        if (i == iterations - 1)
          log_split_gibbs_prob += log(probs(new_cluster_idx));
      } 
      else {
        // Just restore the previous allocation
        data.set_allocation(point_idx, current_cluster);

        // accumulate the log probability of the move usually for the merge move
        log_merge_gibbs_prob += log(probs((current_cluster == ci) ? 0 : 1));
      }
    }
  }
}

double DPSplitMergeWMartinezMena::compute_acceptance_ratio_merge(double likelihood_old_ci,
                                                     double likelihood_old_cj, 
                                                     int old_ci_neighbors, 
                                                     int old_cj_neighbors) {
  /**
   * @brief Compute the log acceptance ratio for a merge move.
   * @param likelihood_old_ci The log likelihood of cluster ci before the merge.
   * @param likelihood_old_cj The log likelihood of cluster cj before the merge.
   * @param old_ci_neighbors The number of neighbors for cluster ci before the merge.
   * @param old_cj_neighbors The number of neighbors for cluster cj before the merge.
   * @return The log acceptance ratio for the merge move.
   */

  // Prior ratio
  int size_old_ci = (launch_state.array() == ci).count();
  int size_old_cj = (launch_state.array() == cj).count();
  double log_acceptance_ratio = -log(params.a);
  log_acceptance_ratio += (S.size() > 0) ? lgamma(S.size()) : 0;
  log_acceptance_ratio -= (size_old_ci > 0) ? lgamma(size_old_ci) : 0;
  log_acceptance_ratio -= (size_old_cj > 0) ? lgamma(size_old_cj) : 0;

  // Add W term to the prior ratio
  int new_cluster_neighbors = cluster_neighbors(ci);
  log_acceptance_ratio += params.coefficient*new_cluster_neighbors;
  log_acceptance_ratio -= params.coefficient*old_ci_neighbors;
  log_acceptance_ratio -= params.coefficient*old_cj_neighbors;

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
  log_acceptance_ratio -= likelihood_old_ci;
  log_acceptance_ratio -= likelihood_old_cj;

  // Proposal ratio
  restricted_gibbs(1, true); // only compute probabilities
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void DPSplitMergeWMartinezMena::merge_move() {
  /**
   * @brief Propose a merge move by combining clusters ci and cj.
   *        Compute the acceptance ratio and decide whether to accept or reject
   * the move.
   */

  // Reset log probabilities
  log_merge_gibbs_prob = 0;

  double likelihood_old_ci = likelihood.cluster_loglikelihood(ci);
  double likelihood_old_cj = likelihood.cluster_loglikelihood(cj);
  int old_ci_neighbors = cluster_neighbors(ci);
  int old_cj_neighbors = cluster_neighbors(cj);

  data.set_allocation(idx_j, ci); // Temporarily assign j to ci for likelihood computation

  // Propose new allocations by merging clusters ci and cj
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state(idx) == cj) {
      data.set_allocation(S(idx), ci); // Merge cj into ci
    }
  }

  // Compute acceptance ratio
  double acceptance_ratio = compute_acceptance_ratio_merge(likelihood_old_ci, likelihood_old_cj, old_ci_neighbors, old_cj_neighbors);

  // Accept or reject the move
  std::uniform_real_distribution<> dis(0.0, 1.0);
  if (log(dis(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
}

void DPSplitMergeWMartinezMena::split_move() {
  /**
   * @brief Propose a split move by dividing cluster ci into two clusters.
   *        Compute the acceptance ratio and decide whether to accept or reject
   * the move.
   */

  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci); // likelihood of the merged cluster old configuration
  int old_cluster_neighbors = cluster_neighbors(ci);

  log_split_gibbs_prob = 0; // reset the log probability of the split move

  data.set_allocation(idx_j, data.get_K()); // Create a new cluster for point j if
                                        // they are in the same cluster
  cj = data.get_cluster_assignment(idx_j);  // Update cj to the new cluster index

  // Allocate randomly points in S to either ci or cj
  std::uniform_int_distribution<> dis(0, 1);
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    int new_cluster = (dis(gen) == 0) ? ci : cj;
    data.set_allocation(S(idx), new_cluster);
  }

  // Perform restricted Gibbs sampling to refine the allocations
  restricted_gibbs(30);

  // Compute acceptance ratio
  double acceptance_ratio = compute_acceptance_ratio_split(likelihood_old_cluster, old_cluster_neighbors);

  // Accept or reject the move
  std::uniform_real_distribution<> dis2(0.0, 1.0);
  if (log(dis2(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
}

double DPSplitMergeWMartinezMena::compute_acceptance_ratio_split(double likelihood_old_cluster, 
                                              int old_cluster_neighbors) {
  /**
   * @brief Compute the log acceptance ratio for a split move.
   * @param likelihood_old_cluster The log likelihood of the original cluster before the split.
   * @param likelihood_old_cluster The log likelihood of the original cluster
   * @return The log acceptance ratio for the split move.
   */

  // Prior ratio
  double log_acceptance_ratio = log(params.a);
  log_acceptance_ratio += (data.get_cluster_size(cj) > 0) ? lgamma(data.get_cluster_size(cj)) : 0;
  log_acceptance_ratio += (data.get_cluster_size(ci) > 0) ? lgamma(data.get_cluster_size(ci)) : 0;
  log_acceptance_ratio -= (S.size() > 0) ? lgamma(S.size()) : 0;
  
  // Add W term to the prior ratio
  int new_ci_neighbors = cluster_neighbors(ci);
  int new_cj_neighbors = cluster_neighbors(cj);
  log_acceptance_ratio += params.coefficient*new_ci_neighbors; 
  log_acceptance_ratio += params.coefficient*new_cj_neighbors;
  log_acceptance_ratio -= params.coefficient*old_cluster_neighbors; 

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
  log_acceptance_ratio += likelihood.cluster_loglikelihood(cj);
  log_acceptance_ratio -= likelihood_old_cluster;

  // Proposal ratio
  log_acceptance_ratio -= log_split_gibbs_prob;

  return log_acceptance_ratio;
}

void DPSplitMergeWMartinezMena::shuffle(){
  /**
   * @brief Perform a shuffle move to refine the allocations of points in S.
   *        This move helps to improve mixing by allowing points to switch
   * clusters.
   */

  log_split_gibbs_prob = 0; // reset the log probability of the split move - use it to store the proposal prob

  if(data.get_K() < 2)
    return; // No point in shuffling if there's only one cluster

  // Choose random two distinct clusters 
  std::uniform_int_distribution<> dis(0, data.get_K() - 1);
  ci = dis(gen);
  do{
    cj = dis(gen);
  } while(cj == ci);
  idx_i = data.get_cluster_assignments(ci)[0]; // pick a random point from cluster ci
  idx_j = data.get_cluster_assignments(cj)[0]; // pick a random point 

  int old_ci_neighbors = cluster_neighbors(ci);
  int old_cj_neighbors = cluster_neighbors(cj);
  double likelihood_old_ci = likelihood.cluster_loglikelihood(ci);
  double likelihood_old_cj = likelihood.cluster_loglikelihood(cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  // Create S and launch_state for points in clusters ci and cj
  int size_ci = data.get_cluster_size(ci);
  int size_cj = data.get_cluster_size(cj);
  int launch_state_size = size_ci + size_cj - 2; // Exclude points i and j from the launch state
  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);
  original_allocations = data.get_allocations(); // Store original allocations in case of rejection

  // Create S and launch_state
  int s_idx = 0;
  for (int idx = 0; idx < data.get_n(); ++idx) {
    if (idx == idx_i || idx == idx_j)
      continue; // Skip points i and j

    int temp_cluster = data.get_cluster_assignment(idx);
    if (temp_cluster == ci || temp_cluster == cj) {
      S(s_idx) = idx;
      launch_state(s_idx) = temp_cluster;
      s_idx++;
    }
  }

  // Use restricted gibbs to refine the allocations
  restricted_gibbs(10);

  // Compute acceptance ratio
  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(likelihood_old_ci, likelihood_old_cj, old_ci_neighbors, old_cj_neighbors, old_ci_size, old_cj_size);

  // Accept or reject the move
  std::uniform_real_distribution<> acceptance_ratio_dis(0.0, 1.0);
  if (log(acceptance_ratio_dis(gen)) > log_acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
}

double DPSplitMergeWMartinezMena::compute_acceptance_ratio_shuffle(double likelihood_old_ci, double likelihood_old_cj, 
                                                                  int old_ci_neighbors, int old_cj_neighbors,
                                                                  int old_ci_size, int old_cj_size) {
  /**
   * @brief Compute the log acceptance ratio for a shuffle move.
   * @return The log acceptance ratio for the shuffle move.
   */

  // Prior ratio
  double log_acceptance_ratio = 0;
  log_acceptance_ratio += (data.get_cluster_size(cj) > 0) ? lgamma(data.get_cluster_size(cj)) : 0;
  log_acceptance_ratio += (data.get_cluster_size(ci) > 0) ? lgamma(data.get_cluster_size(ci)) : 0;
  log_acceptance_ratio -= (old_ci_size > 0) ? lgamma(old_ci_size) : 0;
  log_acceptance_ratio -= (old_cj_size > 0) ? lgamma(old_cj_size) : 0;

  // Add W term to the prior ratio
  int new_ci_neighbors = cluster_neighbors(ci);
  int new_cj_neighbors = cluster_neighbors(cj);
  log_acceptance_ratio += params.coefficient*new_ci_neighbors; 
  log_acceptance_ratio += params.coefficient*new_cj_neighbors;
  log_acceptance_ratio -= params.coefficient*old_ci_neighbors; 
  log_acceptance_ratio -= params.coefficient*old_cj_neighbors; 

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
  log_acceptance_ratio += likelihood.cluster_loglikelihood(cj);
  log_acceptance_ratio -= likelihood_old_ci;
  log_acceptance_ratio -= likelihood_old_cj; 

  // Proposal ratio
  log_acceptance_ratio += log_split_gibbs_prob;
  restricted_gibbs(1, true); // only compute probabilities
  log_acceptance_ratio -= log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void DPSplitMergeWMartinezMena::step() {
  /**
   * @brief Perform a single split-merge MCMC step.
   *        Randomly choose two indices and decide whether to propose a split or
   * merge move.
   */

  choose_indeces();

  if (ci == cj) {
    split_move();
  } else {
    merge_move();
  }

  shuffle();
}
