#include "NGGPW.hpp"
#include <random>

int NGGPW::get_neighbors_obs(int obs_idx, int cls_idx) const {
    /**
     * @brief Returns the number of neighbors for a given observation based on the adjacency matrix W.
     * @param obs_idx The index of the observation.
     * @param cls_idx The index of the cluster to consider for neighbor counting. If -1, counts all neighbors.
     * @return The number of neighbors for the observation.
    */
    
    //int neighbors = (params.W.row(obs_idx).array() * (data.get_allocations().array() == cls_idx).cast<int>().array()).sum();
    int neighbors = 0;
    Eigen::RowVectorXi row = params.W.row(obs_idx);
    for(int i = 0; i < row.size(); ++i) {
        int cluster_i = data.get_cluster_assignment(i);
        if (row(i) == 1 && cluster_i != -1 && cluster_i == cls_idx) {
            neighbors += 1;
        }
    }

    return neighbors;
}

int NGGPW::get_neighbors_cls(int cls_idx, bool old_allo) const {
    /**
     * @brief Returns the total number of neighbors for all observations in a given cluster.
     * @param cls_idx The index of the cluster.
     * @param old_allo If true, uses the old allocations for neighbor counting; otherwise, uses current allocations.
     * @return The total number of neighbors for the cluster.
    */

    Eigen::VectorXi allocations_to_use = old_allo ? old_allocations : data.get_allocations();
    Eigen::VectorXi obs_in_cluster = (allocations_to_use.array() == cls_idx).cast<int>();
    const int total_neighbors = (params.W * obs_in_cluster).sum();
    return total_neighbors;
}

double NGGPW::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {
    /**
     * @brief Computes the log prior probability of assigning a data point to an existing cluster.
     * @param cls_idx The index of the cluster
     * @return The log prior probability of assigning the data point to its current cluster.
    */

    const int cluster_size = data.get_cluster_size(cls_idx);
    double prior = params.coefficient * get_neighbors_obs(obs_idx, cls_idx);
    prior = cluster_size - params.sigma> 0 ? prior + log(cluster_size - params.sigma) : std::numeric_limits<double>::lowest();
    return prior;
}

double NGGPW::gibbs_prior_new_cluster() const {
    /**
     * @brief Computes the log prior probability of assigning a data point to a new cluster.
     * @return The log prior probability of assigning the data point to a new cluster.
    */
    return log_a + params.sigma * log(params.tau + U);
}

double NGGPW::prior_ratio_split(int ci, int cj) const {
    /**
     * @brief Computes the prior ratio for a split operation in a split-merge MCMC algorithm.
     * @param ci the first cluster index involved in the split.
     * @param cj the second cluster index involved in the split.
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

    log_acceptance_ratio += params.coefficient*get_neighbors_cls(ci);
    log_acceptance_ratio += params.coefficient*get_neighbors_cls(cj);
    log_acceptance_ratio -= params.coefficient*get_neighbors_cls(ci, true);

    return log_acceptance_ratio;
}

double NGGPW::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
    /**
     * @brief Computes the prior ratio for a merge operation in a split-merge MCMC algorithm.
     * @param size_old_ci the size of the first cluster before the merge.
     * @param size_old_cj the size of the second cluster before the merge.
     * @return The log prior ratio for the merge operation.
    */

    const int size_merge = size_old_ci + size_old_cj;
    
    double log_acceptance_ratio = - log_a;
    log_acceptance_ratio -= params.sigma * log(params.tau + U);
    log_acceptance_ratio += lgamma(size_merge - params.sigma);
    log_acceptance_ratio -= lgamma(size_old_ci - params.sigma);
    log_acceptance_ratio -= lgamma(size_old_cj - params.sigma);

    // Spatial part
    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    log_acceptance_ratio -= params.coefficient*get_neighbors_cls(old_ci, true);
    log_acceptance_ratio -= params.coefficient*get_neighbors_cls(old_cj,true);
    
    const int new_ci = data.get_allocations()[idx_i]; 
    log_acceptance_ratio += params.coefficient*get_neighbors_cls(new_ci);


    return log_acceptance_ratio;
}

double NGGPW::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) const {
    /**
     * @brief Computes the prior ratio for a shuffle operation in a split-merge MCMC algorithm.
     * @param size_old_ci the size of the first cluster before the shuffle.
     * @param size_old_cj the size of the second cluster before the shuffle.
     * @param ci the first cluster index involved in the shuffle.
     * @param cj the second cluster index involved in the shuffle.
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

// void NGGPW::update_U() {
//     /**
//      * @brief Updates U using Metropolis-Hastings with change of variable V = log(U).
//      * @details Uses a Gaussian proposal kernel with mean V and variance 1/4.
//      * The conditional density f_{V|π}(v) is log-concave, making MH efficient.
//     */
    
//     // Current value V = log(U)
//     double V_current = std::log(U);
    
//     // Propose new V' from N(V, 1/4)
//     std::normal_distribution<double> proposal(V_current, proposal_std); 
//     const double V_proposed = proposal(gen);
    
//     // Compute log conditional densities (unnormalized)
//     const double log_density_current = log_conditional_density_V(V_current);
//     const double log_density_proposed = log_conditional_density_V(V_proposed);
    
//     // Compute acceptance ratio (log scale)
//     double log_acceptance_ratio = log_density_proposed - log_density_current;
    
//     // Accept/reject
//     std::uniform_real_distribution<double> unif(0.0, 1.0);
//     if (std::log(unif(gen)) < log_acceptance_ratio) { // Accept
//         U = std::exp(V_proposed);
//         accepted_U++;
//     }
// }

void NGGPW::update_U() {
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
    while (log_conditional_density_V(L) > log_y) L -= w;
    while (log_conditional_density_V(R) > log_y) R += w;
    
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
        if (V_proposed < V_current) L = V_proposed;
        else R = V_proposed;
    } while (true);
}

double NGGPW::log_conditional_density_V(double v) const {
    /**
     * @brief Computes log of the conditional density f_{V|π}(v).
     * @details f_{V|π}(v) ∝ (e^vn) / (e^v + τ)^{n-a|π|} * e^{-(a/σ)((e^v+τ)^σ - τ^σ)}
     * @param v The value of V = log(U)
     * @return Log of the unnormalized conditional density
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
    const double term3 = - a_over_sigma * (std::pow(exp_v + tau, params.sigma) - tau_power_sigma);
    
    return term1 + term2 + term3;
}

