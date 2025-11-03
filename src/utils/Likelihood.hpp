/**
 * @file Likelihood.hpp
 * @brief Likelihood computation for clustering with cohesion and repulsion
 */

#pragma once

#include "Data.hpp"
#include "Params.hpp"

/**
 * @class Likelihood
 * @brief Computes log-likelihood for clusters based on distance-based cohesion
 * and repulsion
 *
 * This class implements a likelihood model that encourages points within
 * clusters to be close (cohesion) while clusters are pushed apart (repulsion).
 * The model uses gamma priors on distance distributions.
 */
class Likelihood {
private:
  const Data &data; ///< Reference to data containing distances and allocations
  const Params &params; ///< Model parameters

  // Precomputed values for efficiency
  const double lgamma_delta1 = lgamma(params.delta1); ///< log(Γ(δ₁)) - cached for cohesion calculations
  const double log_beta_alpha = log(params.beta) * params.alpha - lgamma(params.alpha); ///< log(β^α / Γ(α)) - cached for cohesion normalization
  const double lgamma_delta2 = lgamma(params.delta2); ///< log(Γ(δ₂)) - cached for repulsion calculations
  const double log_gamma_zeta = log(params.gamma) * params.zeta - lgamma(params.zeta); ///< log(γ^ζ / Γ(ζ)) - cached for repulsion normalization

  /**
   * @brief Computes the cohesion component of the likelihood for a point
   * @details Calculates how well a point fits within a cluster based on
   *          intra-cluster distances using a gamma likelihood
   * @param point_index Index of the point being evaluated
   * @param cluster_index Index of the target cluster
   * @param cls_ass_k Vector of point indices in the cluster
   * @param n_k Number of points in the cluster
   * @return Log-likelihood contribution from cohesion
   */
  double compute_cohesion(int point_index, int cluster_index,
                          const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k,
                          int n_k) const;

  /**
   * @brief Computes the repulsion component of the likelihood for a point
   * @details Calculates the repulsive effect between a point and all other
   * clusters using inter-cluster distances with a gamma likelihood
   * @param point_index Index of the point being evaluated
   * @param cluster_index Index of the cluster containing the point
   * @param cls_ass_k Vector of point indices in the cluster
   * @param n_k Number of points in the cluster
   * @return Log-likelihood contribution from repulsion
   */
  double compute_repulsion(int point_index, int cluster_index,
                           const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k,
                           int n_k) const;

public:
  /**
   * @brief Constructs a Likelihood object
   * @param data Reference to Data object with distances and allocations
   * @param param Reference to model parameters
   */
  Likelihood(const Data &data, const Params &param)
      : data(data), params(param) {}

  /**
   * @brief Calculates the full log-likelihood of a cluster
   * @details Computes both cohesion (within-cluster) and repulsion
   * (between-cluster) components for all points in the specified cluster
   * @param cluster_index Index of the cluster to evaluate
   * @return Total log-likelihood of the cluster
   */
  double cluster_loglikelihood(int cluster_index) const;

  /**
   * @brief Calculates the full log-likelihood of a cluster with explicit
   * assignments
   * @details Overloaded version that accepts cluster assignments directly
   * @param cluster_index Index of the cluster to evaluate
   * @param cls_ass_k Vector of point indices in the cluster
   * @return Total log-likelihood of the cluster
   */
  double cluster_loglikelihood(int cluster_index, const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const;

  /**
   * @brief Calculates the conditional log-likelihood of assigning a point to a
   * cluster
   * @details Computes how likely a point is to belong to a cluster, considering
   *          both its cohesion with the cluster and repulsion from other
   * clusters
   * @param point_index Index of the point to evaluate
   * @param cluster_index Index of the target cluster (can be K for new cluster)
   * @return Conditional log-likelihood of the point given the cluster
   */
  double point_loglikelihood_cond(int point_index, int cluster_index) const;

  /**
   * @brief Calculates the conditional log-likelihood with explicit cluster
   * assignments
   * @details Overloaded version that accepts cluster assignments directly
   * @param point_index Index of the point to evaluate
   * @param cluster_index Index of the target cluster
   * @param cls_ass_k Vector of point indices in the cluster
   * @return Conditional log-likelihood of the point given the cluster
   */
  double point_loglikelihood_cond(int point_index, int cluster_index, const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const;
};