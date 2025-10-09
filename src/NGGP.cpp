/**
 * @file NGGP.cpp
 * @brief Implementation of Normalized Generalized Gamma Process
 *
 * This file contains the implementation of the NGGP class, which provides
 * enhanced flexibility over the standard Dirichlet Process through the
 * generalized gamma construction. NGGP allows for more flexible cluster
 * size distributions and tail behavior control.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "NGGP.hpp"
#include <random>

double NGGP::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {
  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   *
   * For NGGP, this incorporates the discount parameter sigma, giving
   * probability proportional to (n_k - sigma) where n_k is the cluster size.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation (unused in this
   * implementation).
   * @return The log prior probability of assigning the data point to the
   * existing cluster.
   */

  int cluster_size = data.get_cluster_size(cls_idx);
  return (cluster_size - params.sigma > 0)
             ? log(cluster_size - params.sigma)
             : std::numeric_limits<double>::lowest();
}

double NGGP::gibbs_prior_new_cluster() const {
  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   *
   * For NGGP, this depends on the latent variable U and is proportional to
   * alpha * sigma * (tau + U)^sigma.
   * @return The log prior probability of assigning the data point to a new
   * cluster.
   */
  return log_a + params.sigma * log(params.tau + U);
}

double NGGP::prior_ratio_split(int ci, int cj) const {
  /**
   * @brief Computes the prior ratio for a split operation in an NGGP-based
   * split-merge MCMC algorithm.
   *
   * This method accounts for the generalized gamma process prior when computing
   * the acceptance ratio for splitting clusters.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */

  const int n_ci = data.get_cluster_size(ci);
  const int n_cj = data.get_cluster_size(cj);

  double log_acceptance_ratio = 0.0;
  log_acceptance_ratio += log_a;
  log_acceptance_ratio += params.sigma * log(params.tau + U);
  log_acceptance_ratio -= lgamma(n_ci + n_cj - params.sigma);
  log_acceptance_ratio += lgamma(n_ci - params.sigma);
  log_acceptance_ratio += lgamma(n_cj - params.sigma);

  return log_acceptance_ratio;
}

double NGGP::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
  /**
   * @brief Computes the prior ratio for a merge operation in an NGGP-based
   * split-merge MCMC algorithm.
   *
   * This method accounts for the generalized gamma process prior when computing
   * the acceptance ratio for merging clusters.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */

  const int size_merge = size_old_ci + size_old_cj;

  double log_acceptance_ratio = -log_a;
  log_acceptance_ratio -= params.sigma * log(params.tau + U);
  log_acceptance_ratio += lgamma(size_merge - params.sigma);
  log_acceptance_ratio -= lgamma(size_old_ci - params.sigma);
  log_acceptance_ratio -= lgamma(size_old_cj - params.sigma);

  return log_acceptance_ratio;
}

double NGGP::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci,
                                 int cj) const {
  /**
   * @brief Computes the prior ratio for a shuffle operation in an NGGP-based
   * split-merge MCMC algorithm.
   *
   * This method accounts for the generalized gamma process prior when computing
   * the acceptance ratio for shuffling observations between clusters.
   * @param size_old_ci The size of the first cluster before the shuffle.
   * @param size_old_cj The size of the second cluster before the shuffle.
   * @param ci The first cluster index involved in the shuffle.
   * @param cj The second cluster index involved in the shuffle.
   * @return The log prior ratio for the shuffle operation.
   */

  const int n_ci = data.get_cluster_size(ci);
  const int n_cj = data.get_cluster_size(cj);

  double log_acceptance_ratio = 0.0;
  log_acceptance_ratio += lgamma(n_ci - params.sigma);
  log_acceptance_ratio += lgamma(n_cj - params.sigma);
  log_acceptance_ratio -= lgamma(size_old_ci - params.sigma);
  log_acceptance_ratio -= lgamma(size_old_cj - params.sigma);

  return log_acceptance_ratio;
}

// void NGGP::update_U() {
//     /**
//      * @brief Updates U using Metropolis-Hastings with change of variable V =
//      log(U).
//      * @details Uses a Gaussian proposal kernel with mean V and variance 1/4.
//      * The conditional density f_{V|π}(v) is log-concave, making MH
//      efficient.
//     */

//     // Current value V = log(U)
//     double V_current = std::log(U);

//     // Propose new V' from N(V, 1/4)
//     std::normal_distribution<double> proposal(V_current, proposal_std);
//     const double V_proposed = proposal(gen);

//     // Compute log conditional densities (unnormalized)
//     const double log_density_current = log_conditional_density_V(V_current);
//     const double log_density_proposed =
//     log_conditional_density_V(V_proposed);

//     // Compute acceptance ratio (log scale)
//     double log_acceptance_ratio = log_density_proposed - log_density_current;

//     // Accept/reject
//     std::uniform_real_distribution<double> unif(0.0, 1.0);
//     if (std::log(unif(gen)) < log_acceptance_ratio) { // Accept
//         U = std::exp(V_proposed);
//         accepted_U++;
//     }
// }

void NGGP::update_U() {
  /**
   * @brief Updates the latent variable U using slice sampling on V = log(U).
   *
   * This method uses slice sampling to update the latent variable U by working
   * with the log-transformed variable V = log(U) for better numerical
   * stability. The slice sampling algorithm is more efficient than
   * Metropolis-Hastings for this log-concave density.
   */
  // Slice sampling for V = log(U)
  double V_current = std::log(U);

  // Step 1: Sample vertical level
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double log_density_current = log_conditional_density_V(V_current);
  double log_y = log_density_current + std::log(unif(gen));

  // Step 2: Create initial interval around V_current
  double w = 2.0; // Width parameter - tune this
  double L = V_current - w * unif(gen);
  double R = L + w;

  // Step 3: Step out to find slice boundaries
  while (log_conditional_density_V(L) > log_y)
    L -= w;
  while (log_conditional_density_V(R) > log_y)
    R += w;

  // Step 4: Shrinkage sampling
  double V_proposed;
  do {
    V_proposed = L + (R - L) * unif(gen);
    if (log_conditional_density_V(V_proposed) > log_y) {
      U = std::exp(V_proposed);
      accepted_U++;
      return;
    }
    // Shrink interval
    if (V_proposed < V_current)
      L = V_proposed;
    else
      R = V_proposed;
  } while (true);
}

double NGGP::log_conditional_density_V(double v) const {
  /**
   * @brief Computes the log conditional density of V = log(U) given the
   * partition.
   *
   * The conditional density is:
   * f_{V|π}(v) ∝ (e^v)^n / (e^v + τ)^{n-a|π|} * exp(-(a/σ)((e^v+τ)^σ - τ^σ))
   *
   * @param v The value of V = log(U).
   * @return The log of the unnormalized conditional density.
   */

  const double exp_v = std::exp(v);
  const int n = data.get_n();
  const int K = data.get_K();

  // Compute log density components
  // log(e^{vn}) = vn
  const double term1 = v * n;

  // log((e^v + τ)^{n-a|π|}) = (n - a*K) * log(e^v + τ)
  const double term2 = -(n - params.a * K) * std::log(exp_v + tau);

  // -(a/σ)((e^v+τ)^σ - τ^σ)
  const double term3 =
      -a_over_sigma * (std::pow(exp_v + tau, params.sigma) - tau_power_sigma);

  return term1 + term2 + term3;
}
