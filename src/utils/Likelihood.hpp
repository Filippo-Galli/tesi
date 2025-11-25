/**
 * @file Likelihood.hpp
 * @brief Base class for likelihood computation in clustering models
 */

#pragma once

#include "Data.hpp"
#include "Params.hpp"

/**
* @brief Abstract class with likelihood desiderable interfaces
*/

class Likelihood {
protected:
  const Data &data; ///< Reference to Data object with distances and allocations
  const Params &params; ///< Reference to model parameters

public:
  
  Likelihood(const Data &data, const Params &param) : data(data), params(param) {}

  /**
  * @brief Computes the log-likelihood for a cluster
  * @param cluster_index Index of the cluster to evaluate
  * @return Total log-likelihood of the cluster
  * @note this is useful for split-merge algorithms
  */
  virtual double cluster_loglikelihood(int cluster_index) const = 0;

  /**
   * @brief Computes the log-likelihood for a cluster with given assignments
   * @param cluster_index Index of the cluster to evaluate
   * @param cls_ass_k Vector of point indices in the cluster
   * @return Total log-likelihood of the cluster
   * @note this is useful for split-merge algorithms
  */
  virtual double cluster_loglikelihood(int cluster_index, const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const = 0;


  /**
  * @brief Conditional likelihood of a point if set in a particular cluster
  * @param point_index point of which compute the likelihood 
  * @param cluster_index cluster of which computes likelihood
  * @return Conditional likelihood
  * @note useful for gibbs sampling
  */
  virtual double point_loglikelihood_cond(int point_index, int cluster_index) const = 0;

  virtual ~Likelihood() = default;
};