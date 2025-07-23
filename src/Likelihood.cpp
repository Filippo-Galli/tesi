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

double Likelihood::point_loglikelihood_cond(int point_index, int cluster_index) const{
    /**
    @brief This function compute f(Y_l | bold_Y_c)
    */

    double loglik = 0;

    /* -------------------- cohesion part -------------------------- */

    // Useful computation
    int n_k = cluster_index != data.get_K() ? data.get_cluster_size(cluster_index) : 1;
    int pairs = 0.5*n_k*(n_k - 1);
    int pairs_i = 0.5*n_k*(n_k + 1);

    double sum = 0;
    auto cls_ass_k= cluster_index != data.get_K() ? data.get_cluster_assignments(cluster_index) : Eigen::VectorXi::Zero(1);
    if(cluster_index != data.get_K()){
        for (int i = 0; i < n_k; ++i)
            for (int j = i + 1; j < n_k; ++j) {
                sum += data.get_distance(cls_ass_k(i), cls_ass_k(j));
            }
    } 

    double sum_i = sum; // cluster sum of distances + the distances of point_index
    double prod_i = 1;

    if(cluster_index != data.get_K())
        for (int i = 0; i < n_k; ++i){
            sum_i += data.get_distance(point_index, cls_ass_k(i));
            prod_i *= data.get_distance(point_index, cls_ass_k(i));
        }

    // First fraction 
    loglik += (-n_k)*lgamma(params.delta1);
    loglik += (params.delta1 - 1)*log(prod_i);

    // Second fraction
    loglik += lgamma(pairs_i*params.delta1 + params.alpha) - lgamma(pairs*params.delta1 + params.alpha);

    // Third fraction 
    loglik += (pairs*params.delta1 + params.alpha)*log(params.beta + sum);
    loglik -= (pairs_i*params.delta1 + params.alpha)*log(params.beta + sum_i);

    /* -------------------- repulsive part -------------------------- */

    for(int t = 0; t < data.get_K(); ++t){
        if (t == cluster_index) continue;

        int n_t = data.get_cluster_size(t);

        auto cls_ass_t= data.get_cluster_assignments(t);

        sum = 0;
        if(cluster_index != data.get_K()){
            for (int i = 0; i < n_k; ++i)
                for (int j = 0; j < n_t; ++j) {
                    sum += data.get_distance(cls_ass_k(i), cls_ass_t(j));
                }
        }
        
        sum_i = sum; // cluster sum of distances + the distances of point_index
        prod_i = 1;

        for (int i = 0; i < n_t; ++i){
            sum_i += data.get_distance(point_index, cls_ass_t(i));
            prod_i *= data.get_distance(point_index, cls_ass_t(i));
        }

        // First fraction 
        loglik += (-n_t)*lgamma(params.delta2);
        loglik += (params.delta2 - 1)*log(prod_i);

        // Second fraction
        loglik += lgamma((n_k+1)*n_t*params.delta2 + params.zeta) - lgamma(n_k*n_t*params.delta2 + params.zeta);

        // Third fraction 
        loglik += (n_k*n_t*params.delta2 + params.zeta)*log(params.gamma + sum);
        loglik -= ((n_k+1)*n_t*params.delta2 + params.zeta)*log(params.gamma + sum_i);
    }

    return loglik;
}

