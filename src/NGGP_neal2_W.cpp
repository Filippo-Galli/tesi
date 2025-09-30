#include "NGGP_neal2_W.hpp"
#include <algorithm>
#include <random>
#include <Rcpp.h>
using namespace Rcpp;

void NGGPNeal2W::step_1_observation(int index) {
    /**
     * @brief Performs a step in the NGGPNeal2W sampling process.
     * @details This method is responsible for updating the allocations of the data points based on the current state of the model.
     * @param index The index of the data point to update.
    */
    
    // Set unallocated the index
    data.set_allocation(index, -1);

    // Compute adjacentcy edge for each cluster
    Eigen::VectorXi cluster_adjacency = Eigen::VectorXi::Zero(data.get_K());
    
    // Extract the adjacency edges of index
    Eigen::RowVectorXi row = params.W.row(index);
    for(int i = 0; i < row.size(); ++i) {
        int cluster_i = data.get_cluster_assignment(i);
        if (row(i) == 1 && cluster_i != -1) {
            cluster_adjacency(cluster_i) += 1;
        }
    }

    // for each cluster, compute the log likelihood of the point being in that cluster
    std::vector<double> log_likelihoods(data.get_K() + 1, 0.0);
    for (int k = 0; k < data.get_K(); ++k) 
        log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k);

    // Compute the log likelihood of the point being in a new cluster
    log_likelihoods[data.get_K()] = likelihood.point_loglikelihood_cond(index, data.get_K());
    
    // multiply by the prior probability of the cluster
    for (int k = 0; k < data.get_K(); ++k) {
        log_likelihoods[k] += log(data.get_cluster_size(k) - params.sigma) + params.coefficient * cluster_adjacency(k);
    }
    log_likelihoods[data.get_K()] += log(params.a);
    log_likelihoods[data.get_K()] += log(params.sigma) * log(params.tau + U);
    
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
    // Rcout << "[DEBUG] probs for each cluster: ";
    // for (const auto& prob : probs) {
    //     Rcout << prob << " ";
    // }
    // Rcout << std::endl;

    //Rcpp::Rcout << "[DEBUG] Probabilities for each cluster: " << Eigen::Map<Eigen::VectorXd>(probs.data(), probs.size()).transpose() << std::endl;

    // Sample a cluster based on the probabilities
    std::discrete_distribution<int> dist(probs.begin(), probs.end());
    int sampled_cluster = dist(gen);

    // Set the allocation for the data point
    data.set_allocation(index, sampled_cluster);

    //Rcpp::Rcout << "[DEBUG] Data point " << index << " assigned to cluster " << sampled_cluster << std::endl << std::endl;
}

void NGGPNeal2W::update_U() {
    /**
     * @brief Updates U using Metropolis-Hastings with change of variable V = log(U).
     * @details Uses a Gaussian proposal kernel with mean V and variance 1/4.
     * The conditional density f_{V|π}(v) is log-concave, making MH efficient.
    */
    
    // Current value V = log(U)
    double V_current = std::log(U);
    
    // Propose new V' from N(V, 1/4)
    std::normal_distribution<double> proposal(V_current, 0.5); // std = sqrt(1/4) = 0.5
    double V_proposed = proposal(gen);
    
    // Compute log conditional densities (unnormalized)
    double log_density_current = log_conditional_density_V(V_current);
    double log_density_proposed = log_conditional_density_V(V_proposed);
    
    // Compute acceptance ratio (log scale)
    double log_acceptance_ratio = log_density_proposed - log_density_current;
    
    // Accept/reject
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    if (std::log(unif(gen)) < log_acceptance_ratio) { // Accept
        U = std::exp(V_proposed);
    }
}

double NGGPNeal2W::log_conditional_density_V(double v) const {
    /**
     * @brief Computes log of the conditional density f_{V|π}(v).
     * @details f_{V|π}(v) ∝ (e^vn) / (e^v + τ)^{n-a|π|} * e^{-(a/σ)((e^v+τ)^σ - τ^σ)}
     * @param v The value of V = log(U)
     * @return Log of the unnormalized conditional density
    */
    
    double exp_v = std::exp(v);
    int n = data.get_n();
    int K = data.get_K();
    double a = params.a;
    double sigma = params.sigma;
    
    // Compute log density components
    // log(e^{vn}) = vn
    double term1 = v * n;
    
    // log((e^v + τ)^{n-a|π|}) = (n - a*K) * log(e^v + τ)
    double term2 = -(n - a * K) * std::log(exp_v + tau);
    
    // -(a/σ)((e^v+τ)^σ - τ^σ)
    double term3 = -(a / sigma) * (std::pow(exp_v + tau, sigma) - std::pow(tau, sigma));
    
    return term1 + term2 + term3;
}

void NGGPNeal2W::update_tau() {
    /**
     * @brief Updates tau using Metropolis-Hastings with change of variable W = log(tau).
     * @details Uses a Gaussian proposal kernel with mean W and variance 1/4.
     * The conditional density is computed in log scale for numerical stability.
    */
    
    // Current value W = log(tau)
    double W_current = std::log(tau);
    
    // Propose new W' from N(W, 1/4)
    std::normal_distribution<double> proposal(W_current, 0.5); // std = sqrt(1/4) = 0.5
    double W_proposed = proposal(gen);
    
    // Compute log conditional densities (unnormalized)
    double log_density_current = log_conditional_density_W(W_current);
    double log_density_proposed = log_conditional_density_W(W_proposed);
    
    // Compute acceptance ratio (log scale)
    double log_acceptance_ratio = log_density_proposed - log_density_current;
    
    // Accept/reject
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    if (std::log(unif(gen)) < log_acceptance_ratio) {
        tau = std::exp(W_proposed);
    }
    // Otherwise tau remains unchanged
}

double NGGPNeal2W::log_conditional_density_W(double w) const {
    /**
     * @brief Computes log of the conditional density f_{W|a,σ,U,π}(w) where W = log(tau).
     * @details P[dτ|a,σ,U,π] ∝ τ^{α_τ-1} e^{-τβ_τ} * e^{-(a/σ)((U+τ)^σ - τ^σ)} / (τ^σ|π|(U+τ)^{n-σ|π|})
     * @param w The value of W = log(tau)
     * @return Log of the unnormalized conditional density
    */
    
    double tau_val = std::exp(w);
    int n = data.get_n();
    int K = data.get_K();
    double a = params.a;
    double sigma = params.sigma;
    
    // Compute log density components
    // log(τ^{α_τ-1}) = (α_τ - 1) * log(τ) = (α_τ - 1) * w
    double term1 = (alpha_tau - 1) * w;
    
    // log(e^{-τβ_τ}) = -τβ_τ
    double term2 = -tau_val * beta_tau;
    
    // -(a/σ)((U+τ)^σ - τ^σ)
    double term3 = -(a / sigma) * (std::pow(U + tau_val, sigma) - std::pow(tau_val, sigma));
    
    // -log(τ^σ|π|) = -σ|π|log(τ) = -σK * w
    double term4 = -sigma * K * w;
    
    // -log((U+τ)^{n-σ|π|}) = -(n - σK) * log(U + τ)
    double term5 = -(n - sigma * K) * std::log(U + tau_val);
    
    return term1 + term2 + term3 + term4 + term5;
}

void NGGPNeal2W::step() {
    /**
     * @brief Performs a single step of the DP Neal 2 algorithm for all the dataset.
    */

    for (int j = 0; j < data.get_n(); ++j) {
        step_1_observation(j);
    }

    update_U();
    update_tau();
}