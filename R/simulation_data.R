##
## @file simulation_data.R
## @brief MCMC Clustering Analysis with Spatial Dependencies
##
## @details This script performs Bayesian clustering analysis using MCMC methods
## with spatial dependencies. It loads simulated data, sets hyperparameters,
## and runs the MCMC algorithm to identify clusters.
##
## @author Filippo Galli
## @date 2025

## Import utility functions
## @details Sources additional helper functions from utils.R
source("R/utils.R")

## Set random seed for reproducibility
set.seed(44)

# ==============================================================================
# Data Loading
# ==============================================================================

## Path to simulation data folder
## @details Contains simulated data with Natarajan model parameters
sigma <- .18
d <- 10
namefile <- paste0("Natarajan_", sigma, "sigma_", d, "d")
folder <- paste0("simulation_data/", namefile)

## Load simulated datasets
## @param all_data Matrix of observed data points
## @param ground_truth Vector of true cluster assignments
## @param dist_matrix Distance matrix between observations
all_data <- readRDS(file = paste0(folder, "/all_data.rds"))
ground_truth <- readRDS(file = paste0(folder, "/ground_truth.rds"))
dist_matrix <- readRDS(file = paste0(folder, "/dist_matrix.rds"))

# ==============================================================================
# Spatial Adjacency Matrix
# ==============================================================================

## Retrieve spatial adjacency matrix W from distance matrix
## @details Computes binary adjacency matrix based on distance threshold
W <- retrieve_W(dist_matrix)

# ==============================================================================
# C++ Integration
# ==============================================================================

## Load C++ implementation of MCMC algorithm
## @details Sources compiled C++ code for efficient computation
sourceCpp("src/launcher.cpp")

# ==============================================================================
# Hyperparameter Configuration
# ==============================================================================

## Set hyperparameters for the Dirichlet Process model
##
## @param all_data The complete dataset
## @param dist_matrix Distance matrix between observations
## @param k_elbow Number of clusters suggested by elbow method
## @param ground_truth True cluster labels for validation
## @param plot_clustering Whether to plot clustering results
## @param plot_distribution Whether to plot parameter distributions
## @return List containing hyperparameters and initial cluster assignments
## @details Uncomment plot_k_means() to visualize elbow plot for cluster selection
# plot_k_means(dist_matrix, max_k = 15)
hyperparams <- set_hyperparameters(all_data, dist_matrix,
  k_elbow = 3,
  ground_truth = ground_truth, plot_clustering = FALSE, plot_distribution = FALSE
)

# ==============================================================================
# Parameter Object Initialization
# ==============================================================================

## Create parameter object with Dirichlet Process hyperparameters
##
## @param delta1 Shape parameter for the cohesion term
## @param alpha Shape parameter for lambdas of cohesion term
## @param beta Rate parameter for lambdas of cohesion term
## @param delta2 Shape parameter for the repulsion term
## @param gamma Shape parameter for thetas of repulsion term
## @param zeta Scale parameter for thetas of repulsion term
## @param BI Burn-in iterations (2000)
## @param NI Total number of iterations (10000)
## @param a Total mass parameter (1)
## @param sigma Parameter of NGGP prior (0.4)
## @param tau Parameter of NGGP prior (1)
## @param coeff Coefficient for spatial dependence (1)
## @param W Spatial adjacency matrix
## @return Params object containing all model parameters
param <- new(
  Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
  2000, 5000, 0.1, # BI, NI, a,
  0.7, 1, 1, # sigma, tau, coeff spatial dep 
  W # Spatial adjacency matrix
)

# ==============================================================================
# Initial Cluster Allocation
# ==============================================================================

## Initialize cluster allocations
## @details Multiple initialization strategies available:
##
## Option 1: All-in-one allocation (all points in single cluster)
## @code
## hyperparams$initial_clusters <- rep(0, nrow(dist_matrix))
## @endcode
#
# hyperparams$initial_clusters <- rep(0, nrow(dist_matrix))

## Option 2: Sequential allocation (each point in separate cluster)
## @code
## hyperparams$initial_clusters <- seq(0, nrow(dist_matrix) - 1)
## @endcode
#
# hyperparams$initial_clusters <- seq(0, nrow(dist_matrix) - 1)

## Option 3: k-means allocation with different number of centers
## @code
## hyperparams$initial_clusters <- kmeans(all_data, centers = 2, nstart = 25)$cluster - 1
## @endcode
#
# hyperparams$initial_clusters <- kmeans(all_data,
#   centers = 2,
#   nstart = 25
# )$cluster - 1

## Option 4 (Selected): k-means allocation from hyperparameter computation
## @note Cluster indices are converted to 0-based indexing for C++ compatibility
hyperparams$initial_clusters <- hyperparams$initial_clusters - 1

## Display initial cluster allocation summary
print("Initial cluster allocation:")
print(table(hyperparams$initial_clusters))

# ==============================================================================
# MCMC Execution
# ==============================================================================

## Run MCMC algorithm with spatial dependencies
##
## @param dist_matrix Distance matrix between observations
## @param param Parameter object containing all hyperparameters
## @param initial_clusters Initial cluster assignment vector
## @return mcmc_result Object containing posterior samples and cluster assignments
## @details Output is captured to log file for monitoring convergence
log_file <- "mcmc_log.txt"
if (file.exists(log_file)) {
  file.remove(log_file) # Remove previous log file
}

## Execute MCMC and capture console output
results <- capture.output(
  {
    mcmc_result <- mcmc(dist_matrix, param, hyperparams$initial_clusters)
  },
  file = log_file
)

# ==============================================================================
# Save Results (Optional)
# ==============================================================================

## Save results with custom naming convention
##
## @param folder Output directory path
## @param param Parameter object used in analysis
## @param name Custom name for saved output files
## @details Uncomment to save results with specific naming scheme
process <- "NGGPW" # Dirichlet Process
method <- "Neal" # MCMC method used
initialization <- "kmeans" # Initialization strategy)
filename <- paste0(process, "_", method, "_", initialization, "_", sigma, "sigma_", d, "d")
save_with_name(folder, param, filename)

# ==============================================================================
# Visualization
# ==============================================================================

## Plot MCMC results and compare with ground truth
##
## @param mcmc_result Results object from MCMC algorithm
## @param ground_truth True cluster labels as factor
## @param BI Number of burn-in iterations to exclude
## @details Generates diagnostic plots for convergence and cluster assignments
#plot_mcmc_results(mcmc_result, as.factor(ground_truth), BI = param$BI)

# ==============================================================================
# Plot U Trace
# ==============================================================================

## Plot trace of auxiliary variable U over MCMC iterations

# plot(mcmc_result$U, type = "l", xlab = "Iteration", ylab = "U")
# abline(h = mean(mcmc_result$U), col = "red", lty = 2)
# legend("topright", legend = c("Mean U"), col = c("red"), lty = 2)
# titolo <- paste0("Trace of U over MCMC iterations (mean U = ", round(mean(mcmc_result$U), 3), ")")
# title(main = titolo)
