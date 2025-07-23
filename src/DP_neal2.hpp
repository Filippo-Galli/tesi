#pragma once

#include "Likelihood.hpp"
#include "Marginal_sampler.hpp"

class DPNeal2 : public MarginalSampler {
    public:
        DPNeal2(Data& d, Params& p, Likelihood& l) : MarginalSampler(d, p, l) {};

        void step(int index) override;

};