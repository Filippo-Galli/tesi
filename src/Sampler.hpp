#pragma once

#include "Data.hpp"
#include "Params.hpp"
#include "Likelihood.hpp"
#include <random>

class Sampler {
    protected: 
        Data& data;
        Params& params;
        Likelihood& likelihood;
        std::random_device rd;
    public:
        Sampler(Data& d, Params& p, Likelihood& l) : data(d), params(p), likelihood(l) {};

        void virtual step() = 0;
};