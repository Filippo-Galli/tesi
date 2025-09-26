#include "DP_splitmerge.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Rcpp/iostream/Rstreambuf.h"

void DPSplitMerge::choose_indeces() {
    
    std::uniform_int_distribution<> dis(0, data.get_n() - 1);
    
    i = dis(gen);
    do {
        j = dis(gen);
    } while (j == i); // Ensure i and j are distinct

    // Get the clusters of the chosen indices
    ci = data.get_cluster_assignment(i);
    cj = data.get_cluster_assignment(j);

    // Get number of points in clusters ci and cj
    int size_ci = data.get_cluster_size(ci);
    int size_cj = data.get_cluster_size(cj);

    // Initialize launch_state with current allocations of points in clusters ci and cj while S with their indices
    size_t launch_state_size = ci == cj ? size_ci - 2 : size_ci + size_cj - 2; // Exclude points i and j from the launch state   
    
    launch_state.resize(launch_state_size); 
    S.resize(launch_state_size);
    original_allocations = data.get_allocations(); // Store original allocations in case of rejection

    // Properly collect all points from clusters ci and cj
    int s_idx = 0;
    for (int idx = 0; idx < data.get_n(); ++idx) {
        if(idx == i || idx == j) 
            continue; // Skip points i and j

        int temp_cluster = data.get_cluster_assignment(idx);
        if (temp_cluster == ci || temp_cluster == cj) {
            S(s_idx) = idx;
            launch_state(s_idx) = temp_cluster;
            s_idx++;
        }
    }

    Rcpp::Rcout << "[DEBUG] Chosen points for split-merge: " << i << " and " << j << std::endl;
    Rcpp::Rcout << "[DEBUG] Launch state initialized with size " << launch_state_size << ": " << launch_state.transpose() << std::endl;
    Rcpp::Rcout << "[DEBUG] S initialized with indices: " << S.transpose() << std::endl;
    
    // Ensure we collected the expected number of points
    if (s_idx != launch_state_size) { // since s_idx is zero-based
        Rcpp::Rcout << "[ERROR] s_idx = " << s_idx << ", launch_state_size = " << launch_state_size << std::endl;
        throw std::runtime_error("Mismatch in expected cluster sizes during split-merge initialization");
    }
}

void DPSplitMerge::utils_S_filtering(Eigen::VectorXi &S_allocations, int cluster) const {

    // Count points in the specified cluster
    int count = 0;
    for (int i = 0; i < launch_state.size(); ++i) {
        if (launch_state(i) == cluster) {
            count++;
        }
    }
    
    // Resize and fill the output vector
    S_allocations.resize(count);
    int j = 0;
    for (int i = 0; i < launch_state.size(); ++i) {
        if (launch_state(i) == cluster) {
            S_allocations(j++) = S(i);
        }
    }
}

void DPSplitMerge::restricted_gibbs(int iterations, bool only_probabilities){
    /**
    * @brief Perform restricted Gibbs sampling on the points in S to propose new allocations for iter iterations.
    */

    for(int i = 0; i < iterations; ++i){
        for(int idx = 0; idx < S.size(); ++idx){
            int point_idx = S(idx);
            int current_cluster = launch_state(idx);

            // Remove point from its current cluster
            data.set_allocation(point_idx, -1); // Temporarily unassign the point

            // Compute probabilities for each cluster (ci and cj)
            Eigen::Vector2d log_probs;
            int c_i_size_minus_idx = data.get_cluster_size(ci);
            int c_j_size_minus_idx = data.get_cluster_size(cj);
            log_probs(0) = log(c_i_size_minus_idx) + likelihood.point_loglikelihood_cond(point_idx, ci);
            log_probs(1) = log(c_j_size_minus_idx) + likelihood.point_loglikelihood_cond(point_idx, cj);

            // Normalize to get probabilities
            double max_log_prob = log_probs.maxCoeff();
            Eigen::Vector2d probs = (log_probs.array() - max_log_prob).exp();
            probs /= probs.sum();

            if(!only_probabilities){
                // Sample new cluster based on computed probabilities
                std::discrete_distribution<int> dist(probs.data(), probs.data() + probs.size());
                int new_cluster_idx = dist(gen);
                int new_cluster = (new_cluster_idx == 0) ? ci : cj;

                // Assign point to the new cluster
                data.set_allocation(point_idx, new_cluster);

                // if last iteration, accumulate the log probability of the move
                if(i == iterations - 1)
                    log_split_gibbs_prob += log(probs(new_cluster_idx));
            }
            else{
                // Just restore the previous allocation
                data.set_allocation(point_idx, current_cluster);

                // accumulate the log probability of the move usually for the merge move
                log_merge_gibbs_prob += log(probs((current_cluster == ci) ? 0 : 1));
            }
        }
    }
}

double DPSplitMerge::compute_acceptance_ratio_merge(double likelihood_old_ci, double likelihood_old_cj) {
    
    // Prior ratio
    int size_old_ci = (launch_state.array() == ci).count();
    int size_old_cj = (launch_state.array() == cj).count();
    double log_acceptance_ratio = -log(params.alpha);
    log_acceptance_ratio += (S.size() != 0) ? lgamma(S.size()) : 0;
    log_acceptance_ratio -= (size_old_ci != 0) ? lgamma(size_old_ci) : 0;
    log_acceptance_ratio -= (size_old_ci != 0) ? lgamma(size_old_ci) : 0;
    Rcpp::Rcout << "\t[DEBUG] Prior ratio step = " << log_acceptance_ratio << std::endl;

    // Likelihood ratio
    log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
    Rcpp::Rcout << "\t[DEBUG] Likelihood ratio step merged cluster = " << log_acceptance_ratio << std::endl;
    log_acceptance_ratio -= likelihood_old_ci;
    log_acceptance_ratio -= likelihood_old_cj;
    Rcpp::Rcout << "\t[DEBUG] Likelihood step minus old separated clusters = " << log_acceptance_ratio << std::endl;

    // Proposal ratio
    restricted_gibbs(1, true); // only compute probabilities
    log_acceptance_ratio += log_merge_gibbs_prob;

    Rcpp::Rcout << "\t [DEBUG] Merge move: acceptance_ratio = " << log_acceptance_ratio << std::endl;
    return log_acceptance_ratio;
}

void DPSplitMerge::merge_move() {
    // Reset log probabilities
    log_merge_gibbs_prob = 0;

    double likelihood_old_ci = likelihood.cluster_loglikelihood(ci); 
    double likelihood_old_cj = likelihood.cluster_loglikelihood(cj);
    Rcpp::Rcout << "\t [DEBUG] Old likelihoods: ci = " << likelihood_old_ci << ", cj = " << likelihood_old_cj << std::endl;
    
    data.set_allocation(j, ci); // Temporarily assign j to ci for likelihood computation

    // Propose new allocations by merging clusters ci and cj
    for (int idx = 0; idx < launch_state.size(); ++idx) {
        if (launch_state(idx) == cj) {
            data.set_allocation(S(idx), ci); // Merge cj into ci
        }
    }

    // Compute acceptance ratio
    double acceptance_ratio = compute_acceptance_ratio_merge(likelihood_old_ci, likelihood_old_cj);

    // Accept or reject the move
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (log(dis(gen)) > acceptance_ratio) { // move not accepted
        Rcpp::Rcout << "\t [DEBUG] Merge move not accepted." << std::endl;
        Rcpp::Rcout << "\t [DEBUG] Restoring previous allocations, launch state: " << launch_state.transpose() << std::endl;
        data.set_allocations(original_allocations);
    }
    else
        Rcpp::Rcout << "\t [DEBUG] Merge move accepted." << std::endl;
}

void DPSplitMerge::split_move() {

    double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci); // likelihood of the merged cluster old configuration
    Rcpp::Rcout << "\t [DEBUG] Old likelihood of cluster " << ci << " = " << likelihood_old_cluster << std::endl;

    log_split_gibbs_prob = 0; // reset the log probability of the split move

    data.set_allocation(j, data.get_K()); // Create a new cluster for point j if they are in the same cluster
    cj = data.get_cluster_assignment(j); // Update cj to the new cluster index
    Rcpp::Rcout << "\t [DEBUG] number of clusters after adding new cluster: " << data.get_K() << std::endl;

    // Allocate randomly points in S to either ci or cj
    std::uniform_int_distribution<> dis(0, 1);
    for (int idx = 0; idx < launch_state.size(); ++idx) {
        int new_cluster = (dis(gen) == 0) ? ci : cj;
        data.set_allocation(S(idx), new_cluster);
    }
    // Perform restricted Gibbs sampling to refine the allocations
    restricted_gibbs(30); 
    Rcpp::Rcout << "\t [DEBUG] New allocations after split: " << data.get_allocations().transpose() << std::endl;

    // Compute acceptance ratio
    double acceptance_ratio = compute_acceptance_ratio_split(likelihood_old_cluster);

    // Accept or reject the move
    std::uniform_real_distribution<> dis2(0.0, 1.0);
    if (log(dis2(gen)) > acceptance_ratio) { // move not accepted
        Rcpp::Rcout << "\t [DEBUG] Split move not accepted." << std::endl;
        data.set_allocations(original_allocations);
    }
    else
        Rcpp::Rcout << "\t [DEBUG] Split move accepted." << std::endl;
}

double DPSplitMerge::compute_acceptance_ratio_split(double likelihood_old_cluster) {
    
    // Prior ratio
    double log_acceptance_ratio = log(params.alpha);
    log_acceptance_ratio -= (S.size() != 0) ? lgamma(S.size()) : 0;
    log_acceptance_ratio += (data.get_cluster_size(ci) != 0) ? lgamma(data.get_cluster_size(ci)) : 0;
    log_acceptance_ratio += (data.get_cluster_size(cj) != 0) ? lgamma(data.get_cluster_size(cj)) : 0;
    Rcpp::Rcout << "\t[DEBUG] Prior ratio step = " << log_acceptance_ratio << std::endl;
    
    // Likelihood ratio
    log_acceptance_ratio += likelihood.cluster_loglikelihood(ci);
    log_acceptance_ratio += likelihood.cluster_loglikelihood(cj);
    Rcpp::Rcout << "\t[DEBUG] Likelihood ratio step split clusters = " << log_acceptance_ratio << std::endl;
    log_acceptance_ratio -= likelihood_old_cluster;
    Rcpp::Rcout << "\t[DEBUG] Likelihood step minus old merged cluster = " << log_acceptance_ratio << std::endl;

    // Proposal ratio
    log_acceptance_ratio -= log_split_gibbs_prob;
    Rcpp::Rcout << "\t[DEBUG] Split move: acceptance_ratio = " << log_acceptance_ratio << std::endl;

    return log_acceptance_ratio;
}

void DPSplitMerge::step(){

    // Choose two distinct indices i and j
    choose_indeces();

    Rcpp::Rcout << "[DEBUG] Chosen indices: i = " << i << ", j = " << j << " (clusters: " << ci << ", " << cj << ")" << std::endl;
    
    if (ci == cj) {
        Rcpp::Rcout << "[DEBUG] Performing split move." << std::endl;
        // Perform a split move
        split_move();
        Rcpp::Rcout << "[DEBUG] Step completed. Current number of clusters: " << data.get_K()<<std::endl;
    } else {
        Rcpp::Rcout << "[DEBUG] Performing merge move." << std::endl;
        // Perform a merge move
        merge_move();
        Rcpp::Rcout << "[DEBUG] Step completed. Current number of clusters: " << data.get_K()<<std::endl;
    }

    Rcpp::Rcout << "[DEBUG] Allocations: \n " << data.get_allocations().transpose() << std::endl;
    Rcpp::Rcout << "----------------------------------------" << std::endl << std::endl;

}
