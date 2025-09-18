#include "Likelihood.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "cmath"
#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <math.h>

double Likelihood::cluster_loglikelihood(int cluster_index) const {

  auto cls_ass_k = data.get_cluster_assignments(cluster_index);

  return cluster_loglikelihood(cluster_index, cls_ass_k);
}

double
Likelihood::cluster_loglikelihood(int cluster_index,
                                  const Eigen::VectorXi &cls_ass_k) const {

  int K = data.get_K(); // number of cluster
  double rep = 0, coh = 0;
  int n_k =
      data.get_cluster_size(cluster_index); // cluster size of cluster_index
  int pairs = 0;
  bool cohesion_part = false;

  double log_prod = 0;
  double sum = 0;

  if (n_k == 0) {
    std::cerr << "[ERROR] Cluster " << cluster_index << " is empty."
              << std::endl;
    return 0; // Handle empty cluster case
  }
  if (n_k != 1) {
    pairs = n_k * (n_k - 1) / 2;
    cohesion_part = true;
  }

  /* -------------------- Repulsion part -------------------------- */
  for (int t = 0; t < data.get_K(); ++t) {
    if (t == cluster_index)
      continue;

    int n_t = data.get_cluster_size(t); // cluster size of t
    const Eigen::VectorXi &cls_ass_t = data.get_cluster_assignments(t);

    log_prod = 0;
    sum = 0;
    for (int i = 0; i < n_k; ++i) {
      for (int j = 0; j < n_t; ++j) {
        double dist = data.get_distance(cls_ass_k(i), cls_ass_t(j));
        // Add small epsilon to prevent log(0) = -inf
        log_prod += log(dist + 1e-10);
        sum += dist;
      }
    }

    rep += log_prod * (params.delta2 - 1);
    rep -= lgamma(params.delta2) * (n_k * n_t);

    rep += log_gamma_zeta;

    rep += lgamma(n_k * n_t * params.delta2 + params.zeta);
    rep -=
        log(params.gamma + sum) * (n_k * n_t * params.delta2 + params.zeta);
  }

  Rcpp::Rcout << "[DEBUG] repulsive part: " << rep << std::endl;
  
  /* -------------------- cohesion part -------------------------- */

  if (!cohesion_part)
    return rep;

  // Reset log_prod and sum for cohesion calculation
  log_prod = 0;
  sum = 0;

  for (int i = 0; i < n_k; ++i) {
    for (int j = i + 1; j < n_k; ++j) {
      double dist = data.get_distance(cls_ass_k(i), cls_ass_k(j));
      // Add small epsilon to prevent log(0) = -inf
      log_prod += log(dist + 1e-10);
      sum += dist;
    }
  }

  coh += log_prod * (params.delta1 - 1);
  coh -= lgamma(params.delta1) * (pairs);

  // Second fraction
  coh += log_beta_alpha;

  // Third fraction
  coh += lgamma(pairs * params.delta1 + params.alpha);
  coh -=
      log(params.beta + sum) * (pairs * params.delta1 + params.alpha);
  Rcpp::Rcout << "[DEBUG] cohesion part: " << coh << std::endl;

  return rep + coh;
}

double Likelihood::point_loglikelihood_cond(int point_index,
                                            int cluster_index) const {
  // Cluster assignment of the cluster_index
  const auto cls_ass_k = data.get_cluster_assignments(cluster_index);
  // Check if the actual cluster is a new one or not
  int n_k = (cluster_index != data.get_K())
                ? data.get_cluster_size(cluster_index)
                : 0;
  /* -------------------- cohesion part -------------------------- */
  //double coeh = 0;
  double coeh = compute_cohesion(point_index, cluster_index, cls_ass_k, n_k);

  /* -------------------- repulsion part -------------------------- */
  //double rep = 0;
  double rep = compute_repulsion(point_index, cluster_index, cls_ass_k, n_k);

  Rcpp::Rcout << "\t [DEBUG] cluster: " << cluster_index << " point_index: " << point_index << " coeh: " << coeh << " rep: " << rep << std::endl;

  return coeh + rep;
}

double Likelihood::compute_cohesion(int point_index, int cluster_index,
                                    const Eigen::VectorXi cls_ass_k,
                                    int n_k) const {

  double loglik = 0;

  // Early return for empty clusters
  if (n_k == 0) {
    // Rcpp::Rcout << "[WARNING - cohesion] cluster " << cluster_index << " is empty" << std::endl;
    return 0.0;
  }

  double sum_i = 0, log_prod_i = 0;

  // Compute existing cluster sum and prod
  for (int i = 0; i < n_k; ++i) {
    double dist = data.get_distance(point_index, cls_ass_k(i));
    sum_i += dist; // distance of point_index from any other observation in the clsuter
    log_prod_i += log(dist + 1e-10);
  }

  // Cohesion likelihood computation
  loglik += (-n_k) * lgamma_delta1;
  loglik += (params.delta1 - 1) * log_prod_i;
  loglik += log_beta_alpha;
  loglik += lgamma(n_k * params.delta1 + params.alpha);
  loglik -= (n_k * params.delta1 + params.alpha) * log(params.beta + sum_i);

  return loglik;
}

double Likelihood::compute_repulsion(int point_index, int cluster_index,
                                     const Eigen::VectorXi cls_ass_k,
                                     int n_k) const {

  double loglik = 0;

  for (int t = 0; t < data.get_K(); ++t) {
    if (t == cluster_index)
      continue;

    int n_t = data.get_cluster_size(t);
    auto cls_ass_t = data.get_cluster_assignments(t);

    // Skip empty clusters
    if (n_t == 0) {
      continue;
    }

    // Compute inter-cluster sum and point distances in single pass
    double sum_i = 0, log_point_prod = 0;

    for (int j = 0; j < n_t; ++j) {
      double point_dist = data.get_distance(point_index, cls_ass_t(j));   // distance of i from j-element of cluster t
      sum_i += point_dist; // sum_i = sum of distances of i from all points in cluster t
      
      log_point_prod += log(point_dist + 1e-10); // Add small epsilon to prevent log(0) = -inf
    }

    loglik += (-n_t) * lgamma_delta2;
    loglik += (params.delta2 - 1) * log_point_prod;
    loglik += log_gamma_zeta;
    loglik += lgamma(n_t*params.delta2 + params.zeta);
    loglik -= (n_t*params.delta2 + params.zeta) * log(params.gamma + sum_i);
  }

  return loglik;
}