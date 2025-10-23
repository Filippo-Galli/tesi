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
  return (cluster_size - params.sigma > 0) ? log(cluster_size - params.sigma) : std::numeric_limits<double>::lowest();
}

Eigen::VectorXd NGGP::gibbs_prior_existing_clusters(int obs_idx) const {
  /**
   * @brief Computes the log prior probabilities of assigning a data point to
   * all existing clusters.
   *
   * This method incorporates spatial information by considering the number of
   * neighbors in each target cluster when computing the prior probabilities.
   * @param obs_idx The index of the observation to assign.
   * @return A vector of log prior probabilities for assigning the data point to
   * each existing cluster.
   */

  Eigen::VectorXd priors = Eigen::VectorXd::Zero(data.get_K());

  // Compute prior for each existing cluster
  for (int k = 0; k < data.get_K(); ++k) {
    const int cluster_size = data.get_cluster_size(k);
    double prior = cluster_size  - params.sigma > 0 ? log(cluster_size - params.sigma)
                                    : std::numeric_limits<double>::lowest();
    priors(k) = prior;
  }

  return priors;
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
  log_acceptance_ratio -= n_ci + n_cj - params.sigma > 0 ? lgamma(n_ci + n_cj - params.sigma) : 0;
  log_acceptance_ratio += n_ci - params.sigma > 0 ? lgamma(n_ci - params.sigma) : 0;
  log_acceptance_ratio += n_cj - params.sigma > 0 ? lgamma(n_cj - params.sigma) : 0;

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
  log_acceptance_ratio += size_merge - params.sigma > 0 ? lgamma(size_merge - params.sigma) : 0;
  log_acceptance_ratio -= size_old_ci - params.sigma > 0 ? lgamma(size_old_ci - params.sigma) : 0;
  log_acceptance_ratio -= size_old_cj - params.sigma > 0 ? lgamma(size_old_cj - params.sigma) : 0;

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
  log_acceptance_ratio += n_ci - params.sigma > 0 ? lgamma(n_ci - params.sigma) : 0;
  log_acceptance_ratio += n_cj - params.sigma > 0 ? lgamma(n_cj - params.sigma) : 0;
  log_acceptance_ratio -= size_old_ci - params.sigma > 0 ? lgamma(size_old_ci - params.sigma) : 0;
  log_acceptance_ratio -= size_old_cj - params.sigma > 0 ? lgamma(size_old_cj - params.sigma) : 0;

  return log_acceptance_ratio;
}
void NGGP::update_U() {
    /**
     * @brief Updates U using Adaptive Metropolis-Hastings.
     * @details Proposal standard deviation is adapted during burn-in to target
     * an acceptance rate of 0.234 (Roberts et al., 2001).
     */
     
    total_iterations++;

    // Current value
    double U_current = U;

    // Propose new U from N(U_current, proposal_std^2)
    std::normal_distribution<double> proposal(U_current, proposal_std);
    const double U_proposed = proposal(gen);

    // Only accept positive proposals
    if (U_proposed <= 0) {
      return;
    }

    // Compute log conditional densities (unnormalized)
    const double log_density_current = log_conditional_density_U(U_current);
    const double log_density_proposed = log_conditional_density_U(U_proposed);

    // Compute acceptance ratio (log scale)
    double log_acceptance_ratio = log_density_proposed - log_density_current;

    // Accept/reject
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    if (std::log(unif(gen)) < log_acceptance_ratio) {
      U = U_proposed;
      accepted_U++;
    }
}

double NGGP::log_conditional_density_U(double u) const {
  /**
   * @brief Computes the log conditional density of V = log(U) given the
   * partition.
   *
   * The conditional density is:
   * f_{V|π}(v) ∝ (u)^(n-1) / (u + τ)^{n-sigma|π|} * exp(-(a/σ)((u+τ)^σ - τ^σ))
   *
   * @param u The value of U
   * @return The log of the unnormalized conditional density.
   */
;
  const int n = data.get_n();
  const int K = data.get_K();

  // Compute log density components
  // log(e^{vn}) = vn
  const double term1 = (n - 1) * log(u);

  // log((e^v + τ)^{n-sigma|π|}) = - (n - sigma*K) * log(e^v + τ)
  const double term2 = - (n - params.sigma * K) * std::log(u + tau);

  // -(a/σ)((e^v+τ)^σ - τ^σ)
  const double term3 = - a_over_sigma * (std::pow(u + tau, params.sigma) - tau_power_sigma);

  // Rcpp::Rcout << "U: " << u << " Log cond density components: "
  //           << " term1: " << term1
  //           << " term2: " << term2
  //           << " term3: " << term3 << std::endl;
  
  return term1 + term2 + term3;
}
