/**
 * @file zdnam_sampling.cpp
 * @brief ZDNAM (Zero-self Downward Nested Antithetic Modification) implementation
 * 
 * This implements Algorithm 5 from Neal's paper "Modifying Gibbs Sampling to 
 * Avoid Self Transitions" for improved sampling efficiency.
 */

#include "neal_ZDNAM.hpp"
#include <algorithm>
#include <numeric>

/**
 * @brief Compute ZDNAM transition probabilities from current state
 * @param log_probs Log probabilities for each possible state (unnormalized)
 * @param current_k Current cluster assignment (-1 if unallocated)
 * @return Vector of transition probabilities
 */
std::vector<double> Neal3ZDNAM::compute_zdnam_probabilities(const std::vector<double>& log_probs, int current_k) {
    
    int m = log_probs.size();
    std::vector<double> p(m, 0.0);
    
    // Convert log probabilities to normalized probabilities (π)
    double max_loglik = *std::max_element(log_probs.begin(), log_probs.end());
    std::vector<double> pi(m);
    double sum_pi = 0.0;
    for (int i = 0; i < m; ++i) {
        pi[i] = exp(log_probs[i] - max_loglik);
        sum_pi += pi[i];
    }
    for (int i = 0; i < m; ++i) {
        pi[i] /= sum_pi;
    }
    
    // Handle case where current state has zero probability
    if (current_k >= 0 && pi[current_k] <= 0.0) {
        // Fall back to regular Gibbs sampling
        return pi;
    }
    
    // Find ordering by non-increasing probability (σ for DNAM/ZDNAM)
    std::vector<int> sigma(m);
    std::iota(sigma.begin(), sigma.end(), 0);
    std::sort(sigma.begin(), sigma.end(), [&pi](int a, int b) { return pi[a] > pi[b]; });
    
    // Handle case where most probable value has probability >= 1/2
    if (pi[sigma[0]] >= 0.5 || m <= 2) {
        if (current_k == sigma[0]) {
            // Self-transition probability
            p[current_k] = (2.0 * pi[current_k] - 1.0) / pi[current_k];
            // Transition to others
            for (int i = 0; i < m; ++i) {
                if (i != current_k) {
                    p[i] = std::min(1.0, pi[i] / pi[current_k]);
                }
            }
        } else {
            // Always transition to most probable
            p[sigma[0]] = 1.0;
        }
        return p;
    }
    
    // Main ZDNAM algorithm (Algorithm 5 from paper)
    double s = 1.0;  // Sum of probabilities for non-focal values
    double f = 1.0;  // Sum of transition probs from current to non-focal values
    int i = 0;
    
    // Process focal values until we reach current value or need special handling
    while (f > 0 && i < m && sigma[i] != current_k) {
        // Check if next step needs special construction
        if (i + 1 < m) {
            double q = pi[sigma[i]];
            s -= q;
            
            if (i + 2 < m) {
                double q2 = pi[sigma[i + 1]];
                double s2 = std::max(0.0, s - q2);
                
                // Check if special construction is needed
                if (q2 >= s2 && s2 > 0.0) {
                    // Special ZDNAM construction to avoid self-transition
                    double A = (q + q2 - s2) / 2.0;
                    double B = (q - q2 + s2) / (2.0 * s2);
                    double C = (s2 + q2 - q) / (2.0 * s2);
                    
                    // Set transition probabilities
                    p[sigma[i]] = f * B;
                    p[sigma[i + 1]] = f * C;
                    
                    // Remaining values get proportional probabilities
                    for (int j = i + 2; j < m; ++j) {
                        if (sigma[j] == current_k) {
                            for (int k = i + 2; k < m; ++k) {
                                p[sigma[k]] = f * B * pi[sigma[k]] / q;
                            }
                            return p;
                        }
                    }
                    
                    return p;
                }
            }
            
            // Regular DNAM step
            if (s > 0.0) {
                p[sigma[i]] = std::min(f, (q / s) * f);
                f -= p[sigma[i]];
            }
        }
        ++i;
    }
    
    // Handle when current value becomes focal
    if (i < m && sigma[i] == current_k) {
        double q = pi[current_k];
        s -= q;
        
        if (s > 0.0) {
            p[current_k] = 0.0;  // Zero self-transition
            for (int j = i + 1; j < m; ++j) {
                p[sigma[j]] = std::min(f, (pi[sigma[j]] / s) * f);
            }
        } else {
            p[current_k] = f;
        }
    }
    
    return p;
}

/**
 * @brief Sample from log probabilities using ZDNAM
 * @param log_probs Log probabilities (unnormalized)
 * @param current_k Current state (-1 for no current state)
 * @return Sampled index
 */
int Neal3ZDNAM::sample_from_log_probs_zdnam( const std::vector<double>& log_probs, int current_k) {
    
    std::vector<double> probs;
    
    if (current_k >= 0 && current_k < log_probs.size()) {
        // Use ZDNAM transition probabilities
        probs = compute_zdnam_probabilities(log_probs, current_k);
    } else {
        // No current state, use standard normalized probabilities
        double max_loglik = *std::max_element(log_probs.begin(), log_probs.end());
        probs.resize(log_probs.size());
        double sum_probs = 0.0;
        
        for (size_t i = 0; i < log_probs.size(); ++i) {
            probs[i] = exp(log_probs[i] - max_loglik);
            sum_probs += probs[i];
        }
        for (size_t i = 0; i < log_probs.size(); ++i) {
            probs[i] /= sum_probs;
        }
    }
    
    // Sample using roulette wheel
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    double u = unif(gen);
    double cumsum = 0.0;
    
    for (size_t i = 0; i < probs.size(); ++i) {
        cumsum += probs[i];
        if (u < cumsum) {
            return i;
        }
    }
    
    // Fallback for numerical issues
    return probs.size() - 1;
}

/**
 * @brief Updated step_1_observation using ZDNAM
 */
void Neal3ZDNAM::step_1_observation(int index) {
    // Get current allocation (before removing)
    int current_cluster = data.get_cluster_assignment(index);
    
    // Set unallocated
    data.set_allocation(index, -1);
    

    std::vector<double> log_likelihoods(data.get_K() + 1, 0.0);    
    // Existing clusters
    for (int k = 0; k < data.get_K(); ++k) {
        log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k);
        log_likelihoods[k] += process.gibbs_prior_existing_cluster(k, index);
    }
    
    // New cluster
    log_likelihoods[data.get_K()] = likelihood.point_loglikelihood_cond(index, data.get_K());
    log_likelihoods[data.get_K()] += process.gibbs_prior_new_cluster();
    
    // Sample using ZDNAM (pass current cluster to avoid self-transitions)
    int sampled_cluster = sample_from_log_probs_zdnam(
        log_likelihoods, 
        current_cluster
    );
    
    // Set new allocation
    data.set_allocation(index, sampled_cluster);
}

void Neal3ZDNAM::step() {
  /**
   * @brief Performs a single step of the DP Neal 2 algorithm for all the
   * dataset.
   */

  // Create and shuffle a vector of indices
  std::vector<int> indices(data.get_n());
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), gen);

  for (int j = 0; j < indices.size(); ++j) {
    step_1_observation(indices[j]);
  }

  // Rcpp::Rcout << std::endl << "Element per cluster: " << std::endl;
  // for (int k = 0; k < data.get_K(); ++k) {
  //     Rcpp::Rcout << "\tCluster " << k << ": " << data.get_cluster_size(k) <<
  //     std::endl;
  // }
  // Rcpp::Rcout <<
  // "-----------------------------------------------------------------------------"
  // << std::endl;
}