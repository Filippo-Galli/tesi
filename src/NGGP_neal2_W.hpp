#pragma once

#include "Likelihood.hpp"
#include "Params.hpp"
#include "Sampler.hpp"
#include <random>

class NGGPNeal2W : public Sampler {
  private:
    mutable std::mt19937 gen;

    double U = 1;
    double tau = params.tau;

    // Hyperparameters for tau - set to 0.1 both to have a vague prior
    double alpha_tau;
    double beta_tau;
    bool update_tau_hyper = false;

    void step_1_observation(int index);

    void update_U();
    double log_conditional_density_V(double v) const;

    void update_tau();
    double log_conditional_density_W(double w) const;


  public:
    NGGPNeal2W(Data &d, Params &p, Likelihood &l, bool update_tau = false, double alpha = 0.1, double beta = 0.1)
        : Sampler(d, p, l), gen(rd()), update_tau_hyper(update_tau),
        alpha_tau(alpha), beta_tau(beta){};

    double get_U() const { return U; }

    void step() override;
};