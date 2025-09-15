#include "Likelihood.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "cmath"
#include <cmath>
#include <iostream>
#include <math.h>
#include <Rcpp.h>

double Likelihood::cluster_loglikelihood(int cluster_index) const {
  int K = data.get_K();
  double likelihood = 0;
  int n_k = data.get_cluster_size(cluster_index);
  int pairs = 0;
  bool cohesion_part = false;
  // Product of D_ij - work directly in log space
  auto cls_ass_k = data.get_cluster_assignments(cluster_index);

  double log_prod = 0; // Fixed: work in log space from start
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

    int n_t = data.get_cluster_size(t); // Fixed: should be int
    auto cls_ass_t = data.get_cluster_assignments(t);

    log_prod = 0; // Fixed: work in log space from start
    sum = 0;
    for (int i = 0; i < n_k; ++i) {
      for (int j = 0; j < n_t; ++j) {
        double dist = data.get_distance(cls_ass_k(i), cls_ass_t(j));
        log_prod += log(dist); // Fixed: accumulate log directly
        sum += dist;
      }
    }

    likelihood += log_prod * (params.delta2 - 1);
    likelihood += lgamma(params.delta2) * (-n_k * n_t);

    likelihood += log(params.gamma) * params.zeta - lgamma(params.zeta);

    likelihood += lgamma(n_k * n_t * params.delta2 + params.zeta);
    likelihood +=
        log(params.gamma + sum) * (-(n_k * n_t * params.delta2 + params.zeta));
  }
  // std::cout << "[DEBUG] repulsive part: " << likelihood << std::endl;

  /* -------------------- cohesion part -------------------------- */

  if (!cohesion_part)
    return likelihood;

  likelihood += lgamma(params.delta1) * (-pairs);

  // Reset log_prod and sum for cohesion calculation
  log_prod = 0;
  sum = 0;

  for (int i = 0; i < n_k; ++i) {
    for (int j = i + 1; j < n_k; ++j) {
      double dist = data.get_distance(cls_ass_k(i), cls_ass_k(j));
      log_prod += log(dist); // Fixed: accumulate log directly
      sum += dist;
    }
  }

  likelihood += log_prod * (params.delta1 - 1);

  // Second fraction
  likelihood += log(params.beta) * params.alpha - lgamma(params.alpha);

  // Third fraction
  likelihood += lgamma(pairs * params.delta1 + params.alpha);
  likelihood +=
      log(params.beta + sum) * (-(pairs * params.delta1 + params.alpha));

  // std::cout << "[DEBUG] cohesion part: " << likelihood << std::endl;

  return likelihood;
}

double Likelihood::point_loglikelihood_cond(int point_index,
                                            int cluster_index) const {

  double coeh = 0, rep = 0;
  // Cluster assignment of the cluster_index
  const auto cls_ass_k = data.get_cluster_assignments(cluster_index);
  // Check if the actual cluster is a new one or not
  int n_k = (cluster_index != data.get_K())
                ? data.get_cluster_size(cluster_index)
                : 0;

  /* -------------------- cohesion part -------------------------- */
  coeh += compute_cohesion(point_index, cluster_index, cls_ass_k, n_k);
  
  /* -------------------- repulsion part -------------------------- */
  rep += compute_repulsion(point_index, cluster_index, cls_ass_k, n_k);
  
  Rcpp::Rcout << "\t [DEBUG] cluster: " << cluster_index << " coeh: " << coeh << " rep: " << rep << std::endl;

  return coeh + rep;
}

double Likelihood::compute_cohesion(int point_index, int cluster_index,
                                    const Eigen::VectorXi cls_ass_k,
                                    int n_k) const {

  double loglik = 0;

  // Early return for empty clusters
  if (n_k == 0) {
    return 0.0;
  }

  double sum = 0, sum_i = 0, log_prod_i = 0;

  // pairs of distances into the cluster without point_index or with
  int pairs = n_k * (n_k - 1) / 2;
  int pairs_i = n_k * (n_k + 1) / 2;

  // Compute existing cluster sum and prod
  for (int i = 0; i < n_k; ++i) {
    // sum of distances inside the cluster without i
    for (int j = i + 1; j < n_k; ++j) {
      sum += data.get_distance(cls_ass_k(i), cls_ass_k(j));
    }

    double dist = data.get_distance(point_index, cls_ass_k(i));
    sum_i += dist; // distance of i from any other observation in the cluster
    log_prod_i += log(dist);
  }

  sum_i += sum; // sum of distances inside the cluster with i

  // Cohesion likelihood computation
  loglik += (-n_k) * lgamma_delta1;
  loglik += (params.delta1 - 1) * log_prod_i;
  loglik += lgamma(pairs_i * params.delta1 + params.alpha) - lgamma(pairs * params.delta1 + params.alpha);
  loglik += (pairs * params.delta1 + params.alpha) * log(params.beta + sum);
  loglik -= (pairs_i * params.delta1 + params.alpha) * log(params.beta + sum_i);

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
    double sum = 0, sum_i = 0, log_point_prod = 0;

    for (int j = 0; j < n_t; ++j) {
      double point_dist = data.get_distance(
          point_index,
          cls_ass_t(j));       // distance of i from j-element of cluster t
      sum_i += point_dist; // sum_i = sum of distances of i from all
                               // points in cluster t
      log_point_prod += log(point_dist);

      // sum of distances between cluster k and cluster t
      for (int i = 0; i < n_k; ++i) {
        double inter_dist = data.get_distance(cls_ass_k(i), cls_ass_t(j));
        sum += inter_dist; // sum of distances of cluster k from all
                                 // points in cluster t
      }
    }

    // Pre-compute common terms
    sum_i = sum + sum_i;
    double nk_nt_delta2 = n_k * n_t * params.delta2;
    double nk_i_nt_delta2 = (n_k + 1) * n_t * params.delta2;

    loglik += (-n_t) * lgamma_delta2;
    loglik += (params.delta2 - 1) * log_point_prod;
    loglik += lgamma(nk_i_nt_delta2 + params.zeta) - lgamma(nk_nt_delta2 + params.zeta);
    loglik += (nk_nt_delta2 + params.zeta) * log(params.gamma + sum);
    loglik -= (nk_i_nt_delta2 + params.zeta) * log(params.gamma + sum_i);
  }

  return loglik;
}