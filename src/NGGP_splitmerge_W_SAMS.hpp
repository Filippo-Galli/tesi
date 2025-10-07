#pragma once

#include "Likelihood.hpp"
#include "Sampler.hpp"
#include <random>
#include <Eigen/Dense>

class NGGPSplitMergeWSAMS : public Sampler {
  private:
    mutable std::mt19937 gen;

    // Chosen indexes
    int i, j;
    // Clusters of the chosen indexes
    int ci, cj;

    // Launch state
    Eigen::VectorXi launch_state;
    Eigen::VectorXi S; // Indexes of points in clusters ci and cj

    // Original allocations of points in clusters ci and cj
    Eigen::VectorXi original_allocations;

    // NGGP parameters
    double U = 1;
    double tau = params.tau;

    // Hyperparameters for tau - set to 0.1 both to have a vague prior
    double alpha_tau;
    double beta_tau;
    bool update_tau_hyper = false;

    // restricted gibbs sampling prob
    double log_split_gibbs_prob = 0;
    double log_merge_gibbs_prob = 0;

    void update_U();
    double log_conditional_density_V(double v) const;

    void update_tau();
    double log_conditional_density_W(double w) const;

    void choose_indeces();
    void sequential_allocation(int iterations, bool only_probabilities = false);
    
    void split_move();
    double compute_acceptance_ratio_split(double likelihood_old_cluster, int old_cluster_neighbors);
    
    void merge_move();
    double compute_acceptance_ratio_merge(double likelihood_old_ci, double likelihood_old_cj, int old_ci_neighbors, int old_cj_neighbors);

    int cluster_neighbors(int cluster);

  public:
    NGGPSplitMergeWSAMS(Data &d, Params &p, Likelihood &l, bool update_tau = false, double alpha = 0.1, double beta = 0.1)
        : Sampler(d, p, l), gen(rd()), update_tau_hyper(update_tau),
        alpha_tau(alpha), beta_tau(beta){};

    double get_U() const { return U; }

    void step() override;
};