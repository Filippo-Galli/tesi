#include "Likelihood.hpp"
#include "cmath"
#include <cmath>
#include <iostream>
#include <math.h>

double Likelihood::cluster_loglikelihood(int cluster_index) const {
    data.print();

    int K = data.get_K();
    double likelihood = 0;

    /* -------------------- cohesion part -------------------------- */

    int n_k = data.get_cluster_size(cluster_index);
    int pairs = 0.5 * n_k * (n_k - 1);

    likelihood = lgamma(params.delta1) * (-pairs);

    // multiplication of D_ij
    auto cls_ass_k = data.get_cluster_assignments(cluster_index);

    double prod = 1;
    double sum = 0;
    for (int i = 0; i < n_k; ++i)
        for (int j = i + 1; j < n_k; ++j) {
            prod *= data.get_distance(cls_ass_k(i), cls_ass_k(j));
            sum += data.get_distance(cls_ass_k(i), cls_ass_k(j));
        }

    likelihood += log(prod) * (params.delta1 - 1);

    // second fraction
    likelihood += log(params.beta)*params.alpha - lgamma(params.alpha);

    // third fraction
    likelihood += lgamma(pairs * params.delta1 + params.alpha);
    likelihood += log(params.beta + sum) *(-(pairs * params.delta1 + params.alpha));

    //std::cout << "[DEBUG] cohesion part: " << likelihood << std::endl;


    /* -------------------- Repulsion part -------------------------- */
    for (int t = 0; t < data.get_K(); ++t) {
        if (t == cluster_index)
            continue;

        double n_t = data.get_cluster_size(t);
        auto cls_ass_t = data.get_cluster_assignments(t);

        prod = 1;
        sum = 0;
        for (int i = 0; i < n_k; ++i)
            for (int j = 0; j < n_t; ++j) {
                prod *= data.get_distance(cls_ass_k(i), cls_ass_t(j));
                sum += data.get_distance(cls_ass_k(i), cls_ass_t(j));
            }
        
        likelihood += log(prod) * (params.delta2 - 1);
        likelihood += lgamma(params.delta2) * (-n_k*n_t);

        likelihood += log(params.gamma) * params.zeta - lgamma(params.zeta);

        likelihood += lgamma(n_k*n_t * params.delta2 + params.zeta);
        likelihood += log(params.gamma + sum) * (-(n_k*n_t * params.delta2 + params.zeta));

    }

    return likelihood;
}