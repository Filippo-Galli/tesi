#pragma once

#include "Likelihood.hpp"
#include "Sampler.hpp"
#include <random>

class DPNeal2 : public Sampler {
  private:
    mutable std::mt19937 gen;

    void step_1_observation(int index);

  public:
    DPNeal2(Data &d, Params &p, Likelihood &l)
        : Sampler(d, p, l), gen(rd()){};

    void step() override;
};