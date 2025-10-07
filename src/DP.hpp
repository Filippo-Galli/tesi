#pragma once 

#include "Process.hpp"

class DP : public Process {

  public:
    DP(Data &d, Params &p) : Process(d, p){};

    // Gibbs sampling methods
    double gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) override;
    double gibbs_prior_new_cluster() override;

    // Useful for split-merge algorithms
    double prior_ratio_split(int ci, int cj) override;
    double prior_ratio_merge(int size_old_ci, int size_old_cj) override;
    double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) override;

};