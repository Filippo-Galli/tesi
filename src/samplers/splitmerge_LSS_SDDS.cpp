/**
 * @file splitmerge_LSS_SDDS.cpp
 * @brief Implementation of Locality Sensitive Sampling (LSS) with SDDS
 * Split-Merge sampler
 *
 * This file contains the complete implementation of the SplitMerge_LSS_SDDS
 * class, which provides an optimized variant of split-merge sampling using:
 * - Locality sensitive sampling for anchor point selection
 * - SDDS (Smart-split, Dumb-merge, Dumb-split, Smart-merge) for adaptive move
 * selection
 * - Sequential allocation for "smart" proposal generation
 * - Simple random allocation for "dumb" proposals
 * - Adaptive pairing strategy for computational efficiency
 *
 * This approach offers computational advantages for large datasets while
 * maintaining theoretical guarantees of split-merge samplers.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_LSS_SDDS.hpp"
#include <random>

void SplitMerge_LSS_SDDS::choose_indeces(bool similarity) {
  /**
   * @brief Select two distinct data points using locality sensitive sampling
   *
   * @param similarity If true, weights proportional to distance (prefers
   * similar points); if false, weights proportional to 1/distance (prefers
   * dissimilar points)
   *
   * @details First selects idx_i uniformly at random. Then selects idx_j based
   * on distance from idx_i:
   * - If similarity=true: weights proportional to distance (prefers similar
   * points)
   * - If similarity=false: weights proportional to 1/distance (prefers
   * dissimilar points)
   *
   * After selection, identifies clusters ci and cj, and prepares launch_state
   * and S vectors containing all other points from these clusters for the
   * split-merge operation. The vectors are shuffled to ensure random processing
   * order.
   */

  std::uniform_int_distribution<> dis(0, data.get_n() - 1);

  idx_i = dis(gen);

  auto distances = params.D.row(idx_i);
  double distance_sum = 0.0;
  std::vector<double> probs(data.get_n());

  if (similarity) {
    for (auto idx = 0; idx < data.get_n(); ++idx) {
      if (idx == idx_i)
        probs[idx] = 0.0;
      else
        probs[idx] = distances(idx);
      distance_sum += probs[idx];
    }

    for (auto idx = 0; idx < data.get_n(); ++idx) {
      probs[idx] /= distance_sum;
    }
  } else {
    for (auto idx = 0; idx < data.get_n(); ++idx) {
      if (idx == idx_i)
        probs[idx] = 0.0;
      else
        probs[idx] = 1 / distances(idx);
      distance_sum += probs[idx];
    }

    for (auto idx = 0; idx < data.get_n(); ++idx) {
      probs[idx] /= distance_sum;
    }
  }

  std::discrete_distribution<> dis1(probs.begin(), probs.end());
  do {
    idx_j = dis1(gen);
  } while (idx_j == idx_i); // Ensure i and j are distinct

  // Get the clusters of the chosen indices
  ci = data.get_cluster_assignment(idx_i);
  cj = data.get_cluster_assignment(idx_j);

  // Get number of points in clusters ci and cj
  int size_ci = data.get_cluster_size(ci);
  int size_cj = data.get_cluster_size(cj);

  // Initialize launch_state with current allocations of points in clusters ci
  // and cj while S with their indices
  size_t launch_state_size =
      ci == cj ? size_ci - 2
               : size_ci + size_cj -
                     2; // Exclude points i and j from the launch state

  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);
  original_allocations =
      data.get_allocations(); // Store original allocations in case of rejection

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

  // Shuffle launch_state and S in unison using Fisher-Yates algorithm
  // This maintains alignment between the two vectors while shuffling in-place
  for (int i = launch_state_size - 1; i > 0; --i) {
    std::uniform_int_distribution<> dis_shuffle(0, i);
    int j = dis_shuffle(gen);
    if (i != j) {
      std::swap(launch_state(i), launch_state(j));
      std::swap(S(i), S(j));
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

void SplitMerge_LSS_SDDS::sequential_allocation(int iterations,
                                                bool only_probabilities,
                                                bool sequential) {
  /**
   * @brief Perform sequential allocation or restricted Gibbs sampling on points
   * in S
   *
   * @param iterations Number of allocation passes to perform
   * @param only_probabilities If true, compute probabilities without changing
   * allocations
   * @param sequential If true, use sequential allocation; if false, use
   * restricted Gibbs
   *
   * @details This function allocates points between clusters ci and cj using
   * one of two modes:
   * - Sequential allocation (sequential=true): Unallocates all points in S at
   * once, then processes them one by one, computing conditional probabilities
   * based on current state
   * - Restricted Gibbs (sequential=false): Unallocates and reallocates points
   * one at a time
   *
   * For each point, computes conditional log-probabilities for assignment to ci
   * or cj based on likelihood and prior. If only_probabilities=false, samples
   * new assignments and accumulates log_split_gibbs_prob. If
   * only_probabilities=true, retains original assignments and accumulates
   * log_merge_gibbs_prob (for reverse move probability).
   */

  const int S_size = S.size();

  for (int i = 0; i < iterations; ++i) {

    // Unallocate all points in S
    for (int idx = 0; idx < S_size && sequential; ++idx) {
      data.set_allocation(S(idx), -1); // Unallocate point
    }

    for (int idx = 0; idx < S_size; ++idx) {
      int point_idx = S(idx);
      int current_cluster = launch_state(idx);

      if (!sequential) {
        current_cluster = data.get_cluster_assignment(point_idx);
        // Unallocate point only if using restricted Gibbs sampling
        data.set_allocation(point_idx, -1);
      }

      // Compute probabilities for each cluster (ci and cj)
      Eigen::Vector2d log_probs;

      log_probs(0) = likelihood.point_loglikelihood_cond(point_idx, ci);
      log_probs(1) = likelihood.point_loglikelihood_cond(point_idx, cj);

      log_probs(0) += process.gibbs_prior_existing_cluster(ci, point_idx);
      log_probs(1) += process.gibbs_prior_existing_cluster(cj, point_idx);

      // Normalize log probabilities using log-sum-exp trick
      double max_log_prob = log_probs.maxCoeff();
      double log_sum =
          max_log_prob + log((log_probs.array() - max_log_prob).exp().sum());

      if (!only_probabilities) {
        // Compute normalized probabilities only for sampling
        Eigen::Vector2d probs = (log_probs.array() - max_log_prob).exp();
        // probs /= probs.sum();

        // Sample new cluster based on computed probabilities
        std::discrete_distribution<int> dist(probs.data(),
                                             probs.data() + probs.size());
        int new_cluster_idx = dist(gen);
        int new_cluster = (new_cluster_idx == 0) ? ci : cj;

        // Assign point to the new cluster
        data.set_allocation(point_idx, new_cluster);

        // if last iteration, accumulate the log probability of the move
        if (i == iterations - 1)
          log_split_gibbs_prob += log_probs(new_cluster_idx) - log_sum;
      } else {
        // Just restore the previous allocation
        data.set_allocation(point_idx, current_cluster);

        // accumulate the log probability of the move
        int cluster_idx = (current_cluster == ci) ? 0 : 1;
        log_merge_gibbs_prob += log_probs(cluster_idx) - log_sum;
      }
    }
  }
}

double
SplitMerge_LSS_SDDS::compute_acceptance_ratio_merge(double likelihood_old_ci,
                                                    double likelihood_old_cj) {
  /**
   * @brief Compute the log acceptance ratio for a merge move
   *
   * @param likelihood_old_ci Log-likelihood of cluster ci before merge
   * @param likelihood_old_cj Log-likelihood of cluster cj before merge
   * @return Log acceptance ratio for the merge move
   *
   * @details Computes log(α) = log(prior_ratio) + log(likelihood_ratio) +
   * log(proposal_ratio)
   * - Prior ratio: accounts for change from two clusters to one
   * - Likelihood ratio: L(merged_ci) - L(old_ci) - L(old_cj)
   * - Proposal ratio: log_merge_gibbs_prob (probability of reverse split move)
   *   Note: For dumb merge, log_merge_gibbs_prob = 0
   */

  // Prior ratio
  int size_old_ci = (original_allocations.array() == ci).count();
  int size_old_cj = (original_allocations.array() == cj).count();
  double log_acceptance_ratio =
      process.prior_ratio_merge(size_old_ci, size_old_cj);

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
  log_acceptance_ratio -= likelihood_old_ci;
  log_acceptance_ratio -= likelihood_old_cj;

  // Proposal ratio (only included for smart merge)
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS_SDDS::smart_merge_move() {
  /**
   * @brief Propose a smart merge move using sequential allocation
   *
   * @details Merges clusters ci and cj into ci with intelligent proposal:
   * 1. Records old cluster likelihoods
   * 2. Assigns both anchor points (idx_i, idx_j) to ci
   * 3. Initially assigns all other points to ci
   * 4. Uses sequential_allocation(only_probabilities=true) to compute the
   * probability of arriving at the merged configuration (for reverse split
   * probability)
   * 5. Computes acceptance ratio including proposal probability
   * 6. Accepts/rejects via Metropolis-Hastings
   *
   * Higher computational cost than dumb_merge but better proposals through
   * sequential allocation refinement.
   */

  log_merge_gibbs_prob =
      launch_state.size() * rand_split_prob; // reverse dumb split prob
  log_split_gibbs_prob = 0;

  // Get old cluster likelihoods
  double likelihood_old_ci = likelihood.cluster_loglikelihood(ci);
  double likelihood_old_cj = likelihood.cluster_loglikelihood(cj);

  // Assign both anchor points to ci
  data.set_allocation(idx_j, ci);

  // Initially assign all points from both clusters to ci
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    data.set_allocation(S(idx), ci);
  }

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_merge(likelihood_old_ci, likelihood_old_cj);

  // Accept or reject the move
  std::uniform_real_distribution<> dis(0.0, 1.0);
  if (log(dis(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_merge++;
}

void SplitMerge_LSS_SDDS::dumb_merge_move() {
  /**
   * @brief Propose a dumb merge move with direct merging
   *
   * @details Simple merge without sequential allocation:
   * 1. Records old cluster likelihoods
   * 2. Sets log_merge_gibbs_prob = 0 (no proposal probability)
   * 3. Directly assigns idx_j to ci
   * 4. Reassigns all points from cj to ci
   * 5. Computes acceptance ratio without proposal term
   * 6. Accepts/rejects via Metropolis-Hastings
   *
   * Computationally faster than smart_merge but may have lower acceptance rates
   * due to simpler proposal mechanism.
   */

  sequential_allocation(1, true); // reverse smart split proposal probabilities

  double likelihood_old_ci = likelihood.cluster_loglikelihood(ci);
  double likelihood_old_cj = likelihood.cluster_loglikelihood(cj);

  // Direct merge: assign all points to ci
  data.set_allocation(idx_j, ci);

  for (int idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state(idx) == cj) {
      data.set_allocation(S(idx), ci); // Merge cj into ci
    }
  }

  // Compute acceptance ratio (without sequential allocation proposal)
  double acceptance_ratio =
      compute_acceptance_ratio_merge(likelihood_old_ci, likelihood_old_cj);

  // Accept or reject the move
  std::uniform_real_distribution<> dis(0.0, 1.0);
  if (log(dis(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_merge++;
}

void SplitMerge_LSS_SDDS::smart_split_move() {
  /**
   * @brief Propose a smart split move using sequential allocation
   *
   * @details Intelligently splits cluster ci into ci and new cluster cj:
   * 1. Records old single cluster likelihood
   * 2. Creates new cluster cj and assigns idx_j to it (idx_i remains in ci)
   * 3. Randomly initializes allocation of other points to ci or cj
   * 4. Refines via sequential_allocation(1) which:
   *    - Processes points in random order
   *    - Computes conditional probabilities for each assignment
   *    - Accumulates log_split_gibbs_prob for proposal ratio
   * 5. Computes acceptance ratio including proposal probability
   * 6. Accepts/rejects via Metropolis-Hastings
   *
   * Higher computational cost than dumb_split but better mixing through
   * intelligent sequential allocation.
   */

  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);

  log_split_gibbs_prob = 0; // reset the log probability of the split move

  data.set_allocation(idx_j, data.get_K()); // Create a new cluster for point j
  cj = data.get_cluster_assignment(idx_j); // Update cj to the new cluster index

  // Allocate randomly points in S to either ci or cj
  std::uniform_int_distribution<> dis(0, 1);
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    int new_cluster = (dis(gen) == 0) ? ci : cj;
    data.set_allocation(S(idx), new_cluster);
  }

  // Perform sequential allocation to refine the allocations
  sequential_allocation(1);

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  // Accept or reject the move
  std::uniform_real_distribution<> dis2(0.0, 1.0);
  if (log(dis2(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_split++;
}

void SplitMerge_LSS_SDDS::dumb_split_move() {
  /**
   * @brief Propose a dumb split move with random allocation
   *
   * @details Simple split without sequential allocation refinement:
   * 1. Records old single cluster likelihood
   * 2. Sets log_split_gibbs_prob = 0 (no proposal probability)
   * 3. Creates new cluster cj and assigns idx_j to it (idx_i remains in ci)
   * 4. Randomly allocates other points to ci or cj with equal probability
   * (50/50)
   * 5. Computes acceptance ratio without proposal term
   * 6. Accepts/rejects via Metropolis-Hastings
   *
   * Computationally faster than smart_split but may have lower acceptance rates
   * due to purely random initial allocation.
   */

  // old likelihood cluster
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);

  log_split_gibbs_prob =
      launch_state.size() * rand_split_prob; // dumb split prob

  data.set_allocation(idx_j, data.get_K()); // Create a new cluster for point j
  cj = data.get_cluster_assignment(idx_j); // Update cj to the new cluster index

  // Randomly allocate points in S to either ci or cj
  std::uniform_int_distribution<> dis(0, 1);
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    int new_cluster = (dis(gen) == 0) ? ci : cj;
    data.set_allocation(S(idx), new_cluster);
  }

  // Compute acceptance ratio (without sequential allocation proposal)
  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  // Accept or reject the move
  std::uniform_real_distribution<> dis2(0.0, 1.0);
  if (log(dis2(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_split++;
}

double SplitMerge_LSS_SDDS::compute_acceptance_ratio_split(
    double likelihood_old_cluster) {
  /**
   * @brief Compute the log acceptance ratio for a split move
   *
   * @param likelihood_old_cluster Log-likelihood of the original cluster before
   * split
   * @return Log acceptance ratio for the split move
   *
   * @details Computes log(α) = log(prior_ratio) + log(likelihood_ratio) -
   * log(proposal_ratio)
   * - Prior ratio: accounts for change from one cluster to two
   * - Likelihood ratio: L(new_ci) + L(new_cj) - L(old_ci)
   * - Proposal ratio: log_split_gibbs_prob (probability of forward split path)
   *   Note: For dumb split, log_split_gibbs_prob = 0
   */

  // Prior ratio
  double log_acceptance_ratio = process.prior_ratio_split(ci, cj);

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
  log_acceptance_ratio += likelihood.cluster_loglikelihood(cj);
  log_acceptance_ratio -= likelihood_old_cluster;

  // Proposal ratio (only included for smart split)
  // For dumb split, log_split_gibbs_prob = 0, so this term vanishes
  log_acceptance_ratio -= log_split_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS_SDDS::shuffle() {
  /**
   * @brief Perform a shuffle move to redistribute points between two existing
   * clusters
   *
   * @details Refines allocations between clusters ci and cj while maintaining
   * both:
   * 1. Checks if at least 2 clusters exist (returns early if not)
   * 2. Records old cluster likelihoods and sizes
   * 3. Computes reverse proposal probability via
   * sequential_allocation(only_probabilities=true) stored in
   * log_merge_gibbs_prob
   * 4. Performs sequential_allocation(only_probabilities=false) to:
   *    - Reallocate points between ci and cj
   *    - Compute forward proposal probability in log_split_gibbs_prob
   * 5. Computes acceptance ratio with bidirectional proposal probabilities
   * 6. Accepts/rejects via Metropolis-Hastings
   *
   * This move improves mixing by allowing refined redistribution between
   * existing clusters without changing the total number of clusters.
   */

  log_split_gibbs_prob = 0; // reset the log probability of the split move - use
                            // it to store the proposal prob
  log_merge_gibbs_prob = 0;

  if (data.get_K() < 2)
    return; // No point in shuffling if there's only one cluster

  // Get number of points in clusters ci and cj and likelihoods
  double likelihood_old_ci = likelihood.cluster_loglikelihood(ci);
  double likelihood_old_cj = likelihood.cluster_loglikelihood(cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  // Compute probabilities
  sequential_allocation(1, true, true); // only compute probabilities

  // Use restricted gibbs to refine the allocations
  sequential_allocation(1, false, true);

  // Compute acceptance ratio
  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(
      likelihood_old_ci, likelihood_old_cj, old_ci_size, old_cj_size);

  // Accept or reject the move
  std::uniform_real_distribution<> acceptance_ratio_dis(0.0, 1.0);
  if (log(acceptance_ratio_dis(gen)) >
      log_acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_shuffle++;
}

double SplitMerge_LSS_SDDS::compute_acceptance_ratio_shuffle(
    double likelihood_old_ci, double likelihood_old_cj, int old_ci_size,
    int old_cj_size) {
  /**
   * @brief Compute the log acceptance ratio for a shuffle move
   *
   * @param likelihood_old_ci Log-likelihood of cluster ci before shuffle
   * @param likelihood_old_cj Log-likelihood of cluster cj before shuffle
   * @param old_ci_size Size of cluster ci before shuffle
   * @param old_cj_size Size of cluster cj before shuffle
   * @return Log acceptance ratio for the shuffle move
   *
   * @details Computes log(α) = log(prior_ratio) + log(likelihood_ratio) +
   * log(proposal_ratio)
   * - Prior ratio: accounts for cluster size changes (shuffle maintains 2
   * clusters)
   * - Likelihood ratio: L(new_ci) + L(new_cj) - L(old_ci) - L(old_cj)
   * - Proposal ratio: log_merge_gibbs_prob - log_split_gibbs_prob
   *   (ratio of reverse to forward proposal probabilities)
   */

  // Prior ratio
  double log_acceptance_ratio =
      process.prior_ratio_shuffle(old_ci_size, old_cj_size, ci, cj);

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
  log_acceptance_ratio += likelihood.cluster_loglikelihood(cj);
  log_acceptance_ratio -= likelihood_old_ci;
  log_acceptance_ratio -= likelihood_old_cj;

  // Proposal ratio
  log_acceptance_ratio -= log_split_gibbs_prob;
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS_SDDS::choose_clusters_shuffle() {
  /**
   * @brief Select two distinct clusters and anchor points for shuffle move
   *
   * @details For shuffle moves:
   * 1. Returns early if fewer than 2 clusters exist
   * 2. Uniformly selects two distinct clusters ci and cj
   * 3. Uniformly selects random point idx_i from cluster ci
   * 4. Uniformly selects random point idx_j from cluster cj
   * 5. Stores original allocations for potential rejection
   * 6. Prepares launch_state and S vectors containing all points from ci and cj
   *    (excluding idx_i and idx_j which serve as anchors)
   *
   * Unlike choose_indeces(), this method doesn't use locality sensitive
   * sampling since shuffle moves operate on existing cluster structure.
   */

  if (data.get_K() < 2) {
    // std::cout << "Not enough clusters to perform shuffle." << std::endl;
    return;
  }

  std::uniform_int_distribution<> dis(0, data.get_K() - 1);

  ci = dis(gen);
  do {
    cj = dis(gen);
  } while (cj == ci); // Ensure ci and cj are distinct

  // Choose random points from clusters ci and cj
  std::uniform_int_distribution<> dis_idx_i(0, data.get_cluster_size(ci) - 1);
  idx_i = data.get_cluster_assignments(ci)[dis_idx_i(gen)];

  std::uniform_int_distribution<> dis_idx_j(0, data.get_cluster_size(cj) - 1);
  idx_j = data.get_cluster_assignments(cj)[dis_idx_j(gen)];

  // Store original allocations in case of rejection
  original_allocations = data.get_allocations();

  // Pre-allocate launch_state and S
  const int size_ci = data.get_cluster_size(ci);
  const int size_cj = data.get_cluster_size(cj);
  const int launch_state_size =
      size_ci + size_cj - 2; // Exclude points i and j from the launch state
  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);

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
}

void SplitMerge_LSS_SDDS::step() {
  /**
   * @brief Perform a single LSS-SDDS split-merge MCMC step
   *
   * @details Implements the SDDS strategy (Smart-split, Dumb-merge, Dumb-split,
   * Smart-merge):
   *
   * 1. Randomly chooses between two modes with equal probability (50/50):
   *    - Mode 0 (dissimilarity): Selects dissimilar points (1/distance
   * weighting)
   *    - Mode 1 (similarity): Selects similar points (distance weighting)
   *
   * 2. Selects anchor points using choose_indeces() with appropriate similarity
   * flag
   *
   * 3. Updates process state with old allocations and anchor indices
   *
   * 4. Determines and executes move using SDDS pairing:
   *    - **Dissimilarity mode (move_type=0):**
   *      * If ci==cj: **smart_split_move()** - use sequential allocation for
   * quality split
   *      * If ci!=cj: **dumb_merge_move()** - use simple merge for efficiency
   *    - **Similarity mode (move_type=1):**
   *      * If ci!=cj: **dumb_split_move()** - use simple split for efficiency
   *      * If ci==cj: **smart_merge_move()** - use sequential allocation for
   * quality merge
   *
   * 5. Optionally performs shuffle move if shuffle_bool is enabled
   *
   * The SDDS strategy intelligently allocates computational resources:
   * expensive sequential allocation is used for splits when dissimilar points
   * are selected (where splits are likely) and for merges when similar points
   * are selected (where merges are likely), while using simpler proposals
   * elsewhere.
   */

  std::discrete_distribution<> dist({0.5, 0.5});
  const int move_type =
      dist(gen); // 0 for dissimilarity (split), 1 for similarity (merge)
  choose_indeces(move_type == 0 ? false : true);
  process.set_old_allocations(
      data.get_allocations()); // Update old allocations in the process
  process.set_idx_i(idx_i);
  process.set_idx_j(idx_j);

  // Determine which move to perform based on strategy
  if (move_type == 0) {
    // smart split - dumb merge
    if (ci == cj) {
      smart_split_move();
    } else {
      dumb_merge_move();
    }
  } else {
    // smart merge - dumb split
    if (ci != cj) {
      dumb_split_move();
    } else {
      smart_merge_move();
    }
  }

  if (shuffle_bool) {
    choose_clusters_shuffle();
    process.set_old_allocations(
        data.get_allocations()); // Update old allocations in the process
    process.set_idx_i(idx_i);
    process.set_idx_j(idx_j);
    shuffle();
  }
}