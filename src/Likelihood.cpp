#include "Likelihood.hpp"
#include "cmath"
#include <cmath>
#include <iostream>
#include <math.h>

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

  double loglik = 0;

  // Pre-compute gamma values that are reused
  const double lgamma_delta1 = lgamma(params.delta1);
  const double lgamma_delta2 = lgamma(params.delta2);
  const double log_gamma_zeta =
      log(params.gamma) * params.zeta - lgamma(params.zeta);

  int n_k = (cluster_index != data.get_K())
                ? data.get_cluster_size(cluster_index)
                : 0;

  double sum = 0, sum_i = 0, prod_i = 1;

  /* -------------------- New cluster likelihood -------------------------- */
  if (cluster_index == data.get_K()) {
    for (int t = 0; t < data.get_K(); ++t) {
      int n_t = data.get_cluster_size(t);
      auto cls_ass_t = data.get_cluster_assignments(t);

      double point_sum = 0, point_prod = 1;
      for (int j = 0; j < n_t; ++j) {
        double dist = data.get_distance(point_index, cls_ass_t(j));
        point_sum += dist;
        point_prod *= dist;
      }

      loglik += (-n_t) * lgamma_delta2;
      loglik += (params.delta2 - 1) * log(point_prod);
      loglik += log_gamma_zeta;
      loglik += lgamma(n_t * params.delta2 + params.zeta);
      loglik -=
          (n_t * params.delta2 + params.zeta) * log(params.gamma + point_sum);
    }
    return loglik;
  }

  /* -------------------- cohesion part -------------------------- */
  int pairs = n_k * (n_k - 1) / 2;
  int pairs_i = n_k * (n_k + 1) / 2;

  auto cls_ass_k = data.get_cluster_assignments(cluster_index);

  // Compute existing cluster sum once
  for (int i = 0; i < n_k; ++i) {
    for (int j = i + 1; j < n_k; ++j) {
      sum += data.get_distance(cls_ass_k(i), cls_ass_k(j));
    }
    // Compute point-to-cluster distances in same loop
    double dist = data.get_distance(point_index, cls_ass_k(i));
    sum_i += dist;
    prod_i *= dist;
  }

  sum_i += sum;

  // Cohesion likelihood computation
  loglik += (-n_k) * lgamma_delta1;
  loglik += (params.delta1 - 1) * log(prod_i);
  loglik += log(params.beta) * params.alpha - lgamma(params.alpha);
  loglik += lgamma(pairs_i * params.delta1 + params.alpha) -
            lgamma(pairs * params.delta1 + params.alpha);
  loglik += (pairs * params.delta1 + params.alpha) * log(params.beta + sum);
  loglik -= (pairs_i * params.delta1 + params.alpha) * log(params.beta + sum_i);

  /* -------------------- repulsion part -------------------------- */
  for (int t = 0; t < data.get_K(); ++t) {
    if (t == cluster_index)
      continue;

    int n_t = data.get_cluster_size(t);
    auto cls_ass_t = data.get_cluster_assignments(t);

    // Compute inter-cluster sum and point distances in single pass
    double inter_sum = 0, point_sum = 0, point_prod = 1;

    for (int j = 0; j < n_t; ++j) {
      double point_dist = data.get_distance(point_index, cls_ass_t(j));
      point_sum += point_dist;
      point_prod *= point_dist;

      for (int i = 0; i < n_k; ++i) {
        inter_sum += data.get_distance(cls_ass_k(i), cls_ass_t(j));
      }
    }

    // Pre-compute common terms
    double total_sum = inter_sum + point_sum;
    double nk_nt_delta2 = n_k * n_t * params.delta2;
    double nk1_nt_delta2 = (n_k + 1) * n_t * params.delta2;

    loglik += (-n_t) * lgamma_delta2;
    loglik += (params.delta2 - 1) * log(point_prod);
    loglik += log_gamma_zeta;
    loglik += lgamma(nk1_nt_delta2 + params.zeta) -
              lgamma(nk_nt_delta2 + params.zeta);
    loglik += (nk_nt_delta2 + params.zeta) * log(params.gamma + inter_sum);
    loglik -= (nk1_nt_delta2 + params.zeta) * log(params.gamma + total_sum);
  }

  return loglik;
}
