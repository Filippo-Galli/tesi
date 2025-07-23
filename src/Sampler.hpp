#pragma once

#include "Data.hpp"
#include "Params.hpp"

class Sampler {
    private: 
        Data& data;
        Params& params;
    public:
        Sampler(Data& d, Params& p) : data(d), params(p) {}

        virtual void step() = 0;
};