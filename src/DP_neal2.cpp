#include "DP_neal2.hpp"
#include <algorithm>
#include <random>

void DPNeal2::step(int index) {
    /**
     * @brief Performs a step in the DPNeal2 sampling process.
     * @details This method is responsible for updating the allocations of the data points based on the current state of the model.
     * @param index The index of the data point to update.
    */
    
    // Set unallocated the index
    data.set_allocation(index, -1);

    // for each cluster, compute the log likelihood of the point being in that cluster
    std::vector<double> log_likelihoods(data.get_K() + 1, 0.0);
    for (int k = 0; k < data.get_K(); ++k) 
        log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k);

    // Compute the log likelihood of the point being in a new cluster
    log_likelihoods[data.get_K()] = likelihood.point_loglikelihood_cond(index, data.get_K());
    
    // multiply by the prior probability of the cluster
    for (int k = 0; k < data.get_K(); ++k) {
        log_likelihoods[k] += log(data.get_cluster_size(k));
    }
    log_likelihoods[data.get_K()] += log(params.a);

    // DEBUG: Print the log likelihoods for each cluster
    // std::cout << "[DEBUG] Log likelihoods for each cluster: ";
    // for (const auto& log_likelihood : log_likelihoods) {
    //     std::cout << log_likelihood << " ";
    // }
    // std::cout << std::endl;

    // Normalize the log likelihoods
    double max_loglik = *std::max_element(log_likelihoods.begin(), log_likelihoods.end());
    std::vector<double> probs(log_likelihoods);

    for (double& prob : probs) {
        prob = exp(prob - max_loglik);
    }
    double sum_probs = std::accumulate(probs.begin(), probs.end(), 0.0);
    for (double& prob : probs) {
        prob /= sum_probs;
    }

    // DEBUG: Print the probabilities for each cluster
    // std::cout << "[DEBUG]Probabilities for each cluster: ";
    // for (const auto& prob : probs) {
    //     std::cout << prob << " ";
    // }
    // std::cout << std::endl;

    // Sample a cluster based on the probabilities
    std::discrete_distribution<int> dist(probs.begin(), probs.end());
    int sampled_cluster = dist(gen);

    // Set the allocation for the data point
    data.set_allocation(index, sampled_cluster);
}