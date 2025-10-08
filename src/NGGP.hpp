#pragma once 

#include "Process.hpp"
#include <random>

class NGGP : public Process {

  private:
  
  double U = 100;
  double tau = params.tau; 

  void update_U();
  double log_conditional_density_V(double v) const;

  std::random_device rd;
  mutable std::mt19937 gen;

  // Standard deviation for the Gaussian proposal distribution in the MH update of U
  static constexpr double proposal_std = 0.5;

  // Save constant terms for efficiency
  const double a_over_sigma = params.a / params.sigma; 
  const double tau_power_sigma = std::pow(params.tau, params.sigma);
  int accepted_U = 0;

  public:
    NGGP(Data &d, Params &p) : Process(d, p), gen(rd()){};

    // Gibbs sampling methods
    [[nodiscard]] double gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) const override;
    [[nodiscard]] double gibbs_prior_new_cluster() const override;

    // Useful for split-merge algorithms
    [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;
    [[nodiscard]] double prior_ratio_merge(int size_old_ci, int size_old_cj) const override;
    [[nodiscard]] double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) const override;

    void update_params() override { update_U(); };
    double get_U() const { return U; }
    int get_accepted_U() const { return accepted_U; }
};