/**
 * @file Likelihood.cpp
 * @brief Implementation of the Likelihood class for distance-based clustering
 */

#include "Likelihood.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "cmath"
#include <Rcpp.h>
#include <cmath>
#include <math.h>
#include <omp.h>

double Likelihood::cluster_loglikelihood(int cluster_index) const {
  auto cls_ass_k = data.get_cluster_assignments(cluster_index);
  return cluster_loglikelihood(cluster_index, cls_ass_k);
}

double Likelihood::cluster_loglikelihood(int cluster_index,
                                  const Eigen::VectorXi &cls_ass_k) const {
  int K = data.get_K(); // number of clusters
  double rep = 0, coh = 0;
  int n_k = data.get_cluster_size(cluster_index); // cluster size of cluster_index
  int pairs = 0;
  bool cohesion_part = false;

  double log_prod = 0;
  double sum = 0;

  if (n_k == 0) {
    // Silently handle empty cluster case during MCMC sampling
    return 0;
  }
  
  if (n_k != 1) {
    pairs = n_k * (n_k - 1) / 2;
    cohesion_part = true;
  }

  /* -------------------- Repulsion part -------------------------- */
  // Compute repulsive forces between cluster_index and all other clusters
  for (int t = 0; t < data.get_K(); ++t) {
    if (t == cluster_index)
      continue;

    int n_t = data.get_cluster_size(t); // cluster size of t
    const Eigen::VectorXi &cls_ass_t = data.get_cluster_assignments(t);

    log_prod = 0;
    sum = 0;
    
    // Calculate all pairwise distances between clusters
    for (int i = 0; i < n_k; ++i) {
      for (int j = 0; j < n_t; ++j) {
        double dist = data.get_distance(cls_ass_k(i), cls_ass_t(j));
        // Add small epsilon to prevent log(0) = -inf
        log_prod += log(dist + 1e-10);
        sum += dist;
      }
    }

    // Gamma likelihood components for repulsion
    rep += log_prod * (params.delta2 - 1);              // Distance product term
    rep -= lgamma(params.delta2) * (n_k * n_t);         // Normalization

    rep += log_gamma_zeta;                              // Prior term

    rep += lgamma(n_k * n_t * params.delta2 + params.zeta);  // Posterior term
    rep -= log(params.gamma + sum) * (n_k * n_t * params.delta2 + params.zeta);
  }

  /* -------------------- Cohesion part -------------------------- */
  if (!cohesion_part)
    return rep;

  // Reset for cohesion calculation
  log_prod = 0;
  sum = 0;

  // Calculate all pairwise distances within the cluster
  for (int i = 0; i < n_k; ++i) {
    for (int j = i + 1; j < n_k; ++j) {
      double dist = data.get_distance(cls_ass_k(i), cls_ass_k(j));
      // Add small epsilon to prevent log(0) = -inf
      log_prod += log(dist + 1e-10);
      sum += dist;
    }
  }

  // Gamma likelihood components for cohesion
  coh += log_prod * (params.delta1 - 1);         // Distance product term
  coh -= lgamma(params.delta1) * (pairs);        // Normalization

  coh += log_beta_alpha;                         // Prior term

  coh += lgamma(pairs * params.delta1 + params.alpha);  // Posterior term
  coh -= log(params.beta + sum) * (pairs * params.delta1 + params.alpha);

  return rep + coh;
}

double Likelihood::point_loglikelihood_cond(int point_index,
                                            int cluster_index) const {
  // Get cluster assignments for the target cluster
  const auto cls_ass_k = data.get_cluster_assignments(cluster_index);

  // Check if this is a new cluster (index == K) or existing cluster
  int n_k = (cluster_index != data.get_K())
                ? data.get_cluster_size(cluster_index)
                : 0;

  /* -------------------- Cohesion part -------------------------- */
  double coeh = compute_cohesion(point_index, cluster_index, cls_ass_k, n_k);

  /* -------------------- Repulsion part -------------------------- */
  double rep = compute_repulsion(point_index, cluster_index, cls_ass_k, n_k);

  return coeh + rep;
}

double Likelihood::compute_cohesion(int point_index, int cluster_index,
                                    const Eigen::VectorXi cls_ass_k,
                                    int n_k) const {
  double loglik = 0;

  // Early return for empty clusters (new cluster case)
  if (n_k == 0) {
    return 0.0;
  }

  double sum_i = 0, log_prod_i = 0;

  // Compute distances from point_index to all points in the cluster
  for (int i = 0; i < n_k; ++i) {
    double dist = data.get_distance(point_index, cls_ass_k(i));
    sum_i += dist;
    log_prod_i += log(dist + 1e-10);
  }

  // Posterior parameters for gamma distribution
  double alpha_mh = params.alpha + params.delta1 * n_k;
  double beta_mh = params.beta + sum_i;

  // Cohesion likelihood using gamma distribution
  loglik += (-n_k) * lgamma(params.delta1);   // Product normalization
  loglik += (params.delta1 - 1) * log_prod_i; // Distance product term

  // Normalization constant
  loglik += lgamma(alpha_mh);        // Γ(α_mh)
  loglik += log_beta_alpha;          // β^α / Γ(α)
  loglik -= alpha_mh * log(beta_mh); // β_mh^α_mh
  
  return loglik;
}

double Likelihood::compute_repulsion(int point_index, int cluster_index,
                                     const Eigen::VectorXi cls_ass_k,
                                     int n_k) const {
  double loglik = 0;

  // Parallel computation of repulsion from all other clusters
  #pragma omp parallel for reduction(+:loglik)
  for (int t = 0; t < data.get_K(); ++t) {
    if (t == cluster_index)
      continue;

    int n_t = data.get_cluster_size(t);
    auto cls_ass_t = data.get_cluster_assignments(t);

    // Skip empty clusters
    if (n_t == 0) {
      continue;
    }

    double sum_i = 0, log_point_prod = 0;

    // Compute distances from point_index to all points in cluster t
    for (int j = 0; j < n_t; ++j) {
      double point_dist = data.get_distance(point_index, cls_ass_t(j));
      sum_i += point_dist;
      log_point_prod += log(point_dist + 1e-10); // Add epsilon to prevent log(0)
    }

    // Posterior parameters for gamma distribution
    double zeta_mt = params.zeta + params.delta2 * n_t;
    double gamma_mt = params.gamma + sum_i;

    // Repulsion likelihood using gamma distribution
    loglik -= n_t * lgamma_delta2;                  // Product normalization
    loglik += (params.delta2 - 1) * log_point_prod; // Distance product term

    // Normalization constant
    loglik += lgamma(zeta_mt);         // Γ(ζ_mt)
    loglik += log_gamma_zeta;          // γ^ζ / Γ(ζ)
    loglik -= zeta_mt * log(gamma_mt); // γ_mt^ζ_mt
  }

  return loglik;
}