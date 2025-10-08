#pragma once

#include "Likelihood.hpp"
#include "Sampler.hpp"
#include <random>
#include <Eigen/Dense>

class SplitMerge_SAMS : public Sampler {
  private:
    mutable std::mt19937 gen;

    // Chosen indexes
    int idx_i, idx_j;
    // Clusters of the chosen indexes
    int ci, cj;

    // use or not shuffle from Mena and Martinez (2014)
    bool shuffle_bool = false;

    // Launch state
    Eigen::VectorXi launch_state;
    Eigen::VectorXi S; // Indexes of points in clusters ci and cj

    // Original allocations of points in clusters ci and cj
    Eigen::VectorXi original_allocations;

    // restricted gibbs sampling prob
    double log_split_gibbs_prob = 0;
    double log_merge_gibbs_prob = 0;

    void choose_indeces();
    void choose_clusters_shuffle();
    void sequential_allocation(int iterations, bool only_probabilities = false);
    
    void split_move();
    double compute_acceptance_ratio_split(double likelihood_old_cluster);
    
    void merge_move();
    double compute_acceptance_ratio_merge(double likelihood_old_ci, double likelihood_old_cj);

    void shuffle();
    double compute_acceptance_ratio_shuffle(double likelihood_old_ci, double likelihood_old_cj, 
                                           int old_ci_size, int old_cj_size);

  public:
    SplitMerge_SAMS(Data &d, Params &p, Likelihood &l, Process &pr, bool shuffle)
      : Sampler(d, p, l, pr), shuffle_bool(shuffle), gen(rd()){};

    void step() override;
};