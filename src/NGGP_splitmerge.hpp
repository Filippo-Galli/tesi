#pragma once

#include "Likelihood.hpp"
#include "Sampler.hpp"
#include <random>
#include <Eigen/Dense>

class NGGPSplitMerge : public Sampler {
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

    // restricted gibbs sampling prob
    double log_split_gibbs_prob = 0;
    double log_merge_gibbs_prob = 0;

    // NGGP parameters
    double U = 1;
    double tau = params.tau;

    // Hyperparameters for tau - set to 0.1 both to have a vague prior
    double alpha_tau = 0.1;
    double beta_tau = 0.1;


    void update_U();
    double log_conditional_density_V(double v) const;

    void update_tau();
    double log_conditional_density_W(double w) const;

    void choose_indeces();
    void restricted_gibbs(int iterations, bool only_probabilities = false);
    
    void split_move();
    double compute_acceptance_ratio_split(double likelihood_old_cluster);
    
    void merge_move();
    double compute_acceptance_ratio_merge(double likelihood_old_ci, double likelihood_old_cj);

  public:
    NGGPSplitMerge(Data &d, Params &p, Likelihood &l)
        : Sampler(d, p, l), gen(rd()){};

    void step() override;
};