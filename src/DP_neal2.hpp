#pragma once

#include "Likelihood.hpp"
#include "Marginal_sampler.hpp"
#include <random>

class DPNeal2 : public MarginalSampler {
  private:
    mutable std::mt19937 gen;

    void step_1_observation(int index);

  public:
    DPNeal2(Data &d, Params &p, Likelihood &l)
        : MarginalSampler(d, p, l), gen(rd()){};

    void step() override;
};