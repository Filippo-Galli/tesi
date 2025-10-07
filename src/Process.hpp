#pragma once

#include "Likelihood.hpp"
#include "Params.hpp"
#include "Data.hpp"

class Process {
  protected:
    Data &data;
    Params &params;

  public:
    Process(Data &d, Params &p) : data(d), params(p){};

    // Gibbs sampling methods
    virtual double gibbs_prior_existing_cluster(int cls_idx, int obs_idx) = 0;
    virtual double gibbs_prior_new_cluster() = 0;

    // Useful for split-merge algorithms
    virtual double prior_ratio_split(int ci, int cj) = 0;
    virtual double prior_ratio_merge(int size_old_ci, int size_old_cj) = 0;
    virtual double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) = 0;

    virtual ~Process() {}

};