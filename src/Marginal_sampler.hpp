#pragma once

#include "Likelihood.hpp"
#include "Sampler.hpp"

class MarginalSampler : public Sampler {
    public:
    MarginalSampler(Data& d, Params& p, Likelihood& l) : Sampler(d, p, l){};

    void virtual step(int index) = 0;
};