#pragma once

#include "Data.hpp"
#include "Params.hpp"
#include "Likelihood.hpp"

class Sampler {
    protected: 
        Data& data;
        Params& params;
        Likelihood& likelihood;
    public:
        Sampler(Data& d, Params& p, Likelihood& l) : data(d), params(p), likelihood(l) {}
};