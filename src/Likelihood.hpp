#pragma once

#include "Data.hpp"
#include "Params.hpp"

class Likelihood {
    private:
        const Data& data; // Data object containing distances and allocations
        const Params& params; // Parameters for the model

    public:
        Likelihood(const Data& data, const Params& param) : data(data), params(param) {}

        // Calculate the likelihood of the current allocations
        double total_likelihood() const;

        // Calculate the likelihood of a specific cluster
        double cluster_loglikelihood(int cluster_index) const;

        // Calculate the likelihood of a specific point
        double point_loglikelihood(int point_index) const;
};