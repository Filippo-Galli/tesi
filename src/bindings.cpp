#include <RcppEigen.h>
#include "utils/Params.hpp"
#include "utils/Covariates.hpp"
#include "utils/Data.hpp"
#include "utils/Likelihood.hpp"
#include "likelihoods/Natarajan_likelihood.hpp"
#include "utils/Process.hpp"
#include "processes/DP.hpp"
#include "processes/DPW.hpp"
#include "processes/NGGP.hpp"
#include "processes/NGGPW.hpp"
#include "processes/NGGPWx.hpp"
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

// [[Rcpp::depends(RcppEigen)]]

// Factory functions for base classes
// [[Rcpp::export]]
Rcpp::XPtr<Params> create_Params(double delta1, double alpha, double beta, double delta2, 
                                  double gamma, double zeta, int BI, int NI, double a, 
                                  double sigma, double tau, Eigen::MatrixXd D) {
    return Rcpp::XPtr<Params>(new Params(delta1, alpha, beta, delta2, gamma, zeta, 
                                         BI, NI, a, sigma, tau, D), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Covariates> create_Covariates(Eigen::MatrixXi W, double spatial_coefficient, 
                                          Eigen::VectorXi ages, double B, double m, double v,
                                          bool fixed_v, double nu, double S0) {
    return Rcpp::XPtr<Covariates>(new Covariates(W, spatial_coefficient, ages, B, m, v, 
                                                 fixed_v, nu, S0), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Data> create_Data(Rcpp::XPtr<Params> params, Eigen::VectorXi initial_allocations) {
    return Rcpp::XPtr<Data>(new Data(*params, initial_allocations), true);
}

// Factory functions for likelihoods
// [[Rcpp::export]]
Rcpp::XPtr<Natarajan_likelihood> create_Natarajan_likelihood(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params) {
    return Rcpp::XPtr<Natarajan_likelihood>(new Natarajan_likelihood(*data, *params), true);
}

// Factory functions for U samplers
// [[Rcpp::export]]
Rcpp::XPtr<RWMH> create_RWMH(Rcpp::XPtr<Params> params, Rcpp::XPtr<Data> data, bool use_V, 
                             double proposal_sd, bool tuning_enabled) {
    return Rcpp::XPtr<RWMH>(new RWMH(*params, *data, use_V, proposal_sd, tuning_enabled), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<MALA> create_MALA(Rcpp::XPtr<Params> params, Rcpp::XPtr<Data> data, bool use_V, 
                             double proposal_sd, bool tuning_enabled) {
    return Rcpp::XPtr<MALA>(new MALA(*params, *data, use_V, proposal_sd, tuning_enabled), true);
}

// Factory functions for processes
// [[Rcpp::export]]
Rcpp::XPtr<DP> create_DP(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params) {
    return Rcpp::XPtr<DP>(new DP(*data, *params), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<DPW> create_DPW(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates) {
    return Rcpp::XPtr<DPW>(new DPW(*data, *params, *covariates), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGP> create_NGGP(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params, Rcpp::XPtr<U_sampler> u_sampler) {
    return Rcpp::XPtr<NGGP>(new NGGP(*data, *params, *u_sampler), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGPW> create_NGGPW(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates,
                               Rcpp::XPtr<U_sampler> u_sampler) {
    return Rcpp::XPtr<NGGPW>(new NGGPW(*data, *params, *covariates, *u_sampler), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGPWx> create_NGGPWx(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params, Rcpp::XPtr<Covariates> covariates,
                                 Rcpp::XPtr<U_sampler> u_sampler) {
    return Rcpp::XPtr<NGGPWx>(new NGGPWx(*data, *params, *covariates, *u_sampler), true);
}

// Factory functions for samplers
// [[Rcpp::export]]
Rcpp::XPtr<Neal3> create_Neal3(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params, Rcpp::XPtr<Likelihood> likelihood,
                               Rcpp::XPtr<Process> process) {
    return Rcpp::XPtr<Neal3>(new Neal3(*data, *params, *likelihood, *process), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Neal3ZDNAM> create_Neal3ZDNAM(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params,
                                         Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process) {
    return Rcpp::XPtr<Neal3ZDNAM>(new Neal3ZDNAM(*data, *params, *likelihood, *process), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge> create_SplitMerge(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params,
                                         Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process, bool shuffle) {
    return Rcpp::XPtr<SplitMerge>(new SplitMerge(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_SAMS> create_SplitMerge_SAMS(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params,
                                                   Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process,
                                                   bool shuffle) {
    return Rcpp::XPtr<SplitMerge_SAMS>(new SplitMerge_SAMS(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_LSS> create_SplitMerge_LSS(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params,
                                                 Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process,
                                                 bool shuffle) {
    return Rcpp::XPtr<SplitMerge_LSS>(new SplitMerge_LSS(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_LSS_SDDS> create_SplitMerge_LSS_SDDS(Rcpp::XPtr<Data> data, Rcpp::XPtr<Params> params,
                                                           Rcpp::XPtr<Likelihood> likelihood,
                                                           Rcpp::XPtr<Process> process, bool shuffle) {
    return Rcpp::XPtr<SplitMerge_LSS_SDDS>(new SplitMerge_LSS_SDDS(*data, *params, *likelihood, *process, shuffle),
                                           true);
}

// Wrapper functions for methods
// [[Rcpp::export]]
void process_update_params(Rcpp::XPtr<Process> process) { 
    process->update_params(); 
}

// [[Rcpp::export]]
void sampler_step(Rcpp::XPtr<Sampler> sampler) { 
    sampler->step(); 
}

// [[Rcpp::export]]
Eigen::VectorXi data_get_allocations(Rcpp::XPtr<Data> data) { 
    return data->get_allocations(); 
}

// [[Rcpp::export]]
int data_get_K(Rcpp::XPtr<Data> data) { 
    return data->get_K(); 
}

// [[Rcpp::export]]
double u_sampler_get_U(Rcpp::XPtr<U_sampler> u_sampler) { 
    return u_sampler->get_U(); 
}

// [[Rcpp::export]]
double u_sampler_get_acceptance_rate(Rcpp::XPtr<U_sampler> u_sampler) { 
    return u_sampler->get_acceptance_rate(); 
}

// [[Rcpp::export]]
int params_get_BI(Rcpp::XPtr<Params> params) { 
    return params->BI; 
}

// [[Rcpp::export]]
int params_get_NI(Rcpp::XPtr<Params> params) { 
    return params->NI; 
}

// [[Rcpp::export]]
double params_get_a(Rcpp::XPtr<Params> params) {
    return params->a; 
}

// [[Rcpp::export]]
double params_get_sigma(Rcpp::XPtr<Params> params) {
    return params->sigma; 
}

// [[Rcpp::export]]
double params_get_tau(Rcpp::XPtr<Params> params) {
    return params->tau; 
}