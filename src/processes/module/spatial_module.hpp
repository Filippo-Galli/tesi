#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Params.hpp"
#include "Eigen/Dense"
#include <functional>
#include <utility>

/**
 * @class SpatialModule
 * @brief Module providing spatial methods for processes utilizing spatial
 * information.
 *
 * This class offers utility functions to compute neighbor counts based on an
 * adjacency matrix W, facilitating the incorporation of spatial dependencies
 * in clustering algorithms.
 */
class SpatialModule {
protected:
  /**
   * @name Spatial Methods
   * @{
   */

  const Params &params_module;
  const Data &data_module;
  std::function<const Eigen::VectorXi &()> old_allocations_provider;

  /**
   * @brief Returns the number of neighbors for a given observation in a
   * specific cluster.
   *
   * This method counts the neighbors of an observation based on the adjacency
   * matrix W, considering only neighbors that belong to the specified cluster.
   * @param obs_idx The index of the observation.
   * @param cls_idx The index of the cluster to consider for neighbor counting.
   * @return The number of neighbors for the observation in the specified
   * cluster.
   */
  int get_neighbors_obs(int obs_idx, int cls_idx) const;

  /**
   * @brief Returns the number of neighbors for a given observation
   * regardless of cluster membership.
   * This method counts the total number of neighbors of an observation
   * based on the adjacency matrix W.
   * @param obs_idx The index of the observation.
   * @return The total number of neighbors for the observation for all clusters.
   */
  Eigen::VectorXi get_neighbors_obs(int obs_idx) const;

  /**
   * @brief Returns the total number of neighbors for all observations in a
   * given cluster.
   *
   * This method computes the sum of all neighbor connections within a cluster,
   * which is used in the spatial component of the prior calculations.
   * @param cls_idx The index of the cluster.
   * @param old_allo If true, uses the old allocations for neighbor counting;
   * otherwise, uses current allocations (default: false).
   * @return The total number of neighbors for the cluster.
   */
  int get_neighbors_cls(int cls_idx, bool old_allo = false) const;
  /** @} */

public:
  SpatialModule(
      const Params &params_, const Data &data_,
      std::function<const Eigen::VectorXi &()> old_alloc_provider = {})
      : params_module(params_), data_module(data_),
        old_allocations_provider(std::move(old_alloc_provider)) {}
};