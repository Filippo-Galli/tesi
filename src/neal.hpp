#pragma once

#include "Likelihood.hpp"
#include "Sampler.hpp"
#include "Process.hpp"
#include <random>

class Neal3 : public Sampler {
  private:
    mutable std::mt19937 gen;
    void step_1_observation(int index);

  public:
    Neal3(Data &d, Params &p, Likelihood &l, Process &pr)
        : Sampler(d, p, l, pr), gen(rd()){};

    void step() override;
};