#include <RcppEigen.h>
#include "utils/Params.hpp"
#include "utils/Covariates.hpp"
#include "utils/Data.hpp"
#include "utils/Data_wClusterInfo.hpp"
#include "utils/ClusterInfo.hpp"
#include "utils/Likelihood.hpp"
#include "likelihoods/Natarajan_likelihood.hpp"
#include "utils/Process.hpp"
#include "processes/DP.hpp"
#include "processes/DPW.hpp"
#include "processes/NGGP.hpp"
#include "processes/NGGPW.hpp"
#include "processes/NGGPWx.hpp"
#include "processes/NGGPWxCache.hpp"
#include "processes/caches/Covariate_cache.hpp"
#include "samplers/U_sampler/U_sampler.hpp"
#include "samplers/U_sampler/RWMH.hpp"
#include "samplers/U_sampler/MALA.hpp"
#include "utils/Sampler.hpp"
#include "samplers/neal.hpp"
#include "samplers/neal_ZDNAM.hpp"
#include "samplers/splitmerge.hpp"
#include "samplers/splitmerge_SAMS.hpp"
#include "samplers/splitmerge_LSS.hpp"
#include "samplers/splitmerge_LSS_SDDS.hpp"

// Helper functions to handle inheritance with XPtrs
Data *get_data_ptr(SEXP sexp) {
    // Try Data
    try {
        Rcpp::XPtr<Data> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    // Try Data_wClusterInfo
    try {
        Rcpp::XPtr<Data_wClusterInfo> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    Rcpp::stop("Expected external pointer to Data or Data_wClusterInfo");
    return nullptr;
}

ClusterInfo *get_cluster_info_ptr(SEXP sexp) {
    // Try ClusterInfo
    try {
        Rcpp::XPtr<ClusterInfo> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    // Try Covariate_cache
    try {
        Rcpp::XPtr<Covariate_cache> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    Rcpp::stop("Expected external pointer to ClusterInfo or Covariate_cache");
    return nullptr;
}

// [[Rcpp::depends(RcppEigen)]]

// Factory functions for base classes
// [[Rcpp::export]]
Rcpp::XPtr<Params> create_Params(double delta1, double alpha, double beta, double delta2, double gamma, double zeta,
                                 int BI, int NI, double a, double sigma, double tau, Eigen::MatrixXd D) {
    return Rcpp::XPtr<Params>(new Params(delta1, alpha, beta, delta2, gamma, zeta, BI, NI, a, sigma, tau, D), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Covariates> create_Covariates(Eigen::MatrixXi W, double spatial_coefficient, Eigen::VectorXd ages, double B,
                                         double m, double v, bool fixed_v, double nu, double S0) {
    return Rcpp::XPtr<Covariates>(new Covariates(W, spatial_coefficient, ages, B, m, v, fixed_v, nu, S0), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Data> create_Data(Rcpp::XPtr<Params> params, Eigen::VectorXi initial_allocations) {
    return Rcpp::XPtr<Data>(new Data(*params, initial_allocations), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Covariate_cache> create_Covariate_cache(Rcpp::XPtr<Covariates> covariates,
                                                   Eigen::VectorXi &initial_allocations) {
    return Rcpp::XPtr<Covariate_cache>(new Covariate_cache(*covariates, initial_allocations), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Data_wClusterInfo> create_Data_wClusterInfo(Rcpp::XPtr<Params> params, SEXP cluster_info_sexp,
                                                       Eigen::VectorXi initial_allocations) {
    ClusterInfo *cluster_info = get_cluster_info_ptr(cluster_info_sexp);
    return Rcpp::XPtr<Data_wClusterInfo>(new Data_wClusterInfo(*params, *cluster_info, initial_allocations), true);
}

// Factory functions for likelihoods
// [[Rcpp::export]]
Rcpp::XPtr<Natarajan_likelihood> create_Natarajan_likelihood(SEXP data_sexp, Rcpp::XPtr<Params> params) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Natarajan_likelihood>(new Natarajan_likelihood(*data, *params), true);
}

// Factory functions for U samplers
// [[Rcpp::export]]
Rcpp::XPtr<RWMH> create_RWMH(Rcpp::XPtr<Params> params, SEXP data_sexp, bool use_V, double proposal_sd,
                             bool tuning_enabled) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<RWMH>(new RWMH(*params, *data, use_V, proposal_sd, tuning_enabled), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<MALA> create_MALA(Rcpp::XPtr<Params> params, SEXP data_sexp, bool use_V, double proposal_sd,
                             bool tuning_enabled) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<MALA>(new MALA(*params, *data, use_V, proposal_sd, tuning_enabled), true);
}

// Factory functions for processes
// [[Rcpp::export]]
Rcpp::XPtr<DP> create_DP(SEXP data_sexp, Rcpp::XPtr<Params> params) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<DP>(new DP(*data, *params), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<DPW> create_DPW(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<DPW>(new DPW(*data, *params, *covariates), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGP> create_NGGP(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<U_sampler> u_sampler) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<NGGP>(new NGGP(*data, *params, *u_sampler), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGPW> create_NGGPW(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates,
                               Rcpp::XPtr<U_sampler> u_sampler) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<NGGPW>(new NGGPW(*data, *params, *covariates, *u_sampler), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGPWx> create_NGGPWx(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates,
                                 Rcpp::XPtr<U_sampler> u_sampler) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<NGGPWx>(new NGGPWx(*data, *params, *covariates, *u_sampler), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGPWxCache> create_NGGPWxCache(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates,
                                           Rcpp::XPtr<U_sampler> u_sampler, SEXP cov_cache_sexp) {
    Data *data = get_data_ptr(data_sexp);
    // We know it's a Covariate_cache, but let's use the helper or just cast if we are sure.
    // The helper returns ClusterInfo*, we need Covariate_cache*.
    // So we should probably just use XPtr<Covariate_cache> here if we are sure.
    // But wait, if the user passed a generic ClusterInfo (if possible), this would fail.
    // But NGGPWxCache REQUIRES Covariate_cache.
    // So we should keep XPtr<Covariate_cache> in the signature?
    // YES. Because NGGPWxCache constructor takes Covariate_cache&.
    // But wait, create_Covariate_cache returns XPtr<Covariate_cache>.
    // So passing it here works fine.
    // The only issue is if we changed create_Covariate_cache to return ClusterInfo.
    // We didn't.
    // So we only need to change data argument.

    // Wait, I changed create_Data_wClusterInfo to take SEXP.
    // But create_NGGPWxCache takes Covariate_cache.
    // So I only need to change data argument here.

    // BUT, I should check if I need to change cov_cache argument too?
    // No, create_Covariate_cache returns Covariate_cache.

    return Rcpp::XPtr<NGGPWxCache>(
        new NGGPWxCache(*data, *params, *covariates, *u_sampler, *Rcpp::XPtr<Covariate_cache>(cov_cache_sexp)), true);
}

// Factory functions for samplers
// [[Rcpp::export]]
Rcpp::XPtr<Neal3> create_Neal3(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Likelihood> likelihood,
                               Rcpp::XPtr<Process> process) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Neal3>(new Neal3(*data, *params, *likelihood, *process), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Neal3ZDNAM> create_Neal3ZDNAM(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Likelihood> likelihood,
                                         Rcpp::XPtr<Process> process) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Neal3ZDNAM>(new Neal3ZDNAM(*data, *params, *likelihood, *process), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge> create_SplitMerge(SEXP data_sexp, Rcpp::XPtr<Params> params, Rcpp::XPtr<Likelihood> likelihood,
                                         Rcpp::XPtr<Process> process, bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge>(new SplitMerge(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_SAMS> create_SplitMerge_SAMS(SEXP data_sexp, Rcpp::XPtr<Params> params,
                                                   Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process,
                                                   bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge_SAMS>(new SplitMerge_SAMS(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_LSS> create_SplitMerge_LSS(SEXP data_sexp, Rcpp::XPtr<Params> params,
                                                 Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process,
                                                 bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge_LSS>(new SplitMerge_LSS(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_LSS_SDDS> create_SplitMerge_LSS_SDDS(SEXP data_sexp, Rcpp::XPtr<Params> params,
                                                           Rcpp::XPtr<Likelihood> likelihood,
                                                           Rcpp::XPtr<Process> process, bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge_LSS_SDDS>(new SplitMerge_LSS_SDDS(*data, *params, *likelihood, *process, shuffle),
                                           true);
}

// Wrapper functions for methods
// [[Rcpp::export]]
void process_update_params(Rcpp::XPtr<Process> process) { process->update_params(); }

// [[Rcpp::export]]
void sampler_step(Rcpp::XPtr<Sampler> sampler) { sampler->step(); }

// [[Rcpp::export]]
Eigen::VectorXi data_get_allocations(SEXP data_sexp) {
    Data *data = get_data_ptr(data_sexp);
    return data->get_allocations();
}

// [[Rcpp::export]]
int data_get_K(SEXP data_sexp) {
    Data *data = get_data_ptr(data_sexp);
    return data->get_K();
}

// [[Rcpp::export]]
double u_sampler_get_U(Rcpp::XPtr<U_sampler> u_sampler) { return u_sampler->get_U(); }

// [[Rcpp::export]]
double u_sampler_get_acceptance_rate(Rcpp::XPtr<U_sampler> u_sampler) { return u_sampler->get_acceptance_rate(); }

// [[Rcpp::export]]
int params_get_BI(Rcpp::XPtr<Params> params) { return params->BI; }

// [[Rcpp::export]]
int params_get_NI(Rcpp::XPtr<Params> params) { return params->NI; }

// [[Rcpp::export]]
double params_get_a(Rcpp::XPtr<Params> params) { return params->a; }

// [[Rcpp::export]]
double params_get_sigma(Rcpp::XPtr<Params> params) { return params->sigma; }

// [[Rcpp::export]]
double params_get_tau(Rcpp::XPtr<Params> params) { return params->tau; }

// [[Rcpp::export]]
void cluster_info_set_allocation(Rcpp::XPtr<ClusterInfo> cluster_info, int index, int cluster, int old_cluster) {
    cluster_info->set_allocation(index, cluster, old_cluster);
}

// [[Rcpp::export]]
void cluster_info_recompute(Rcpp::XPtr<ClusterInfo> cluster_info, int K, Eigen::VectorXi allocations) {
    cluster_info->recompute(K, allocations);
}

// [[Rcpp::export]]
void data_wclusterinfo_set_allocation(Rcpp::XPtr<Data_wClusterInfo> data, int index, int cluster) {
    data->set_allocation(index, cluster);
}

// [[Rcpp::export]]
void data_wclusterinfo_set_allocations(Rcpp::XPtr<Data_wClusterInfo> data, Eigen::VectorXi new_allocations) {
    data->set_allocations(new_allocations);
}