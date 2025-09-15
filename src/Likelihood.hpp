#pragma once

#include "Data.hpp"
#include "Params.hpp"

class Likelihood {
    private:
        const Data& data; // Data object containing distances and allocations
        const Params& params; // Parameters for the model

        // Useful fixed values
        const double lgamma_delta1 = lgamma(params.delta1);
        const double log_gamma_alpha = log(params.beta) * params.alpha - lgamma(params.alpha);
        const double lgamma_delta2 = lgamma(params.delta2);
        const double log_gamma_zeta = log(params.gamma) * params.zeta - lgamma(params.zeta);

        double compute_cohesion(int point_index, int cluster_index, const Eigen::VectorXi cls_ass_k, int n_k) const;
        double compute_repulsion(int point_index, int cluster_index, const Eigen::VectorXi cls_ass_k, int n_k) const;

    public:
        Likelihood(const Data& data, const Params& param) : data(data), params(param) {}

        // Calculate the likelihood of a specific cluster
        double cluster_loglikelihood(int cluster_index) const;

        // Calculate the conditional likelihood of a specific point
        double point_loglikelihood_cond(int point_index, int cluster_index) const;
};