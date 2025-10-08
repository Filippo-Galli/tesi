#pragma once 

#include "Process.hpp"
#include <Eigen/Dense>

class DPW : public Process {

  public:
    DPW(Data &d, Params &p) : Process(d, p){};

    // Gibbs sampling methods
    [[nodiscard]] double gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const  override;
    [[nodiscard]] double gibbs_prior_new_cluster() const override;

    // Useful for split-merge algorithms
    [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;
    [[nodiscard]] double prior_ratio_merge(int size_old_ci, int size_old_cj) const override;
    [[nodiscard]] double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) const override;

    // spatial methods
    int get_neighbors_obs(int obs_idx, int cls_idx) const;
    int get_neighbors_cls(int cls_idx, bool old_allo = false) const;

    // Null update_param
    void update_params() override {
      return;
    };

};