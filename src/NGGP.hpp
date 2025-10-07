#pragma once 

#include "Process.hpp"
#include <random>

class NGGP : public Process {

  private:
  
  double U = 1;
  double tau = params.tau; 

  void update_U();
  double log_conditional_density_V(double v) const;

  std::random_device rd;
  mutable std::mt19937 gen;

  public:
    NGGP(Data &d, Params &p) : Process(d, p), gen(rd()){};

    // Gibbs sampling methods
    double gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) override;
    double gibbs_prior_new_cluster() override;

    // Useful for split-merge algorithms
    double prior_ratio_split(int ci, int cj) override;
    double prior_ratio_merge(int size_old_ci, int size_old_cj) override;
    double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) override;

    void update_params() override { update_U(); };
    double get_U() const { return U; }
};