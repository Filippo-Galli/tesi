#pragma once

#include "Data.hpp"
#include "Params.hpp"
#include "Likelihood.hpp"
#include "Process.hpp"
#include <random>

class Sampler {
    protected: 
        Data& data;
        Params& params;
        Likelihood& likelihood;
        Process& process;
        
        std::random_device rd;
    public:
        Sampler(Data& d, Params& p, Likelihood& l, Process &pr) : data(d), params(p), likelihood(l), process(pr) {};

        void virtual step() = 0;
};