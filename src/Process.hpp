#pragma once

#include "Likelihood.hpp"
#include "Params.hpp"
#include "Data.hpp"
#include <Eigen/Dense>

class Process {
  protected:
    Data &data;
    Params &params;
    Eigen::VectorXi old_allocations; // to store old allocations in case of rejection
    int idx_i, idx_j; // indices of the two observations involved in the split-merge move

  public:
    Process(Data &d, Params &p) : data(d), params(p){ old_allocations = data.get_allocations(); };

    // Gibbs sampling methods
    virtual double gibbs_prior_existing_cluster(int cls_idx, int obs_idx) = 0;
    virtual double gibbs_prior_new_cluster() = 0;

    // Useful for split-merge algorithms
    virtual double prior_ratio_split(int ci, int cj) = 0;
    virtual double prior_ratio_merge(int size_old_ci, int size_old_cj) = 0;
    virtual double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) = 0;

    // Set old allocations
    void set_old_allocations(const Eigen::VectorXi &new_allocations) { old_allocations = new_allocations; };
    
    void set_idx_i(int i) { idx_i = i; };
    void set_idx_j(int j) { idx_j = j; };
    
    virtual ~Process() {}

};