source("R/utils.R")
source("R/utils_plot.R")

## Set random seed for reproducibility
set.seed(44)

##############################################################################
# Data Loading ====
##############################################################################

## Load simulated data
# sigma <- .2
# d <- 50
# namefile <- paste0("Natarajan_", sigma, "sigma_", d, "d")
# folder <- paste0("simulation_data/", namefile)
# all_data <- readRDS(file = paste0(folder, "/all_data.rds"))
# ground_truth <- readRDS(file = paste0(folder, "/ground_truth.rds"))
# dist_matrix <- readRDS(file = paste0(folder, "/dist_matrix.rds"))

## Load real data
files_folder <- "real_data/CA"

files <- list.files(files_folder)
file_chosen <- files[5]
dist_matrix <- readRDS(file = paste0(files_folder, "/", file_chosen))
plot_distance(dist_matrix)

if (min(dist_matrix) < 0) {
  dist_matrix <- dist_matrix + abs(min(dist_matrix))
}
diag(dist_matrix) <- 0

##############################################################################
# Spatial Adjacency Matrix ====
##############################################################################

## Retrieve spatial adjacency matrix W from distance matrix
# W <- retrieve_W(dist_matrix)
W <- readRDS(file = paste0(files_folder, "/adj_matrix.rds"))

# Check is W is symmetric
if (!isSymmetric(W)) {
  warning("W is not symmetric!")
}
##############################################################################
# C++ Integration ====
##############################################################################

## Load C++ implementation of MCMC algorithm
sourceCpp("src/launcher.cpp")
cat("✅ C++ code compiled successfully!\n\n")

##############################################################################
# Hyperparameter Configuration ====
##############################################################################

# Plot k-means elbow method to help set hyperparameters
plot_k_means(dist_matrix, max_k = 10)

# Set hyperparameters based on distance matrix
hyperparams <- set_hyperparameters(dist_matrix,
  k_elbow = 3, plot_clustering = FALSE, plot_distribution = FALSE
)

##############################################################################
# Parameter Object Initialization ====
##############################################################################

param <- new(
  Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
  10000, 10000, 1, # BI, NI, a,
  0.1, 1, 1, # sigma, tau, coeff spatial dependence
  dist_matrix, W # distance matrix, Spatial adjacency matrix
)

##############################################################################
# Initial Cluster Allocation ====
##############################################################################

hyperparams$initial_clusters <- hyperparams$initial_clusters - 1

## Display initial cluster allocation summary
print("Initial cluster allocation:")
print(table(hyperparams$initial_clusters))

##############################################################################
# MCMC Execution ====
##############################################################################

log_file <- "mcmc_log.txt"
if (file.exists(log_file)) {
  file.remove(log_file) # Remove previous log file
}

## Execute MCMC and capture console output
results <- capture.output(
  {
    mcmc_result <- mcmc(param, hyperparams$initial_clusters)
  },
  file = log_file
)

##############################################################################
# Save Results (Optional) ====
##############################################################################
# file_chosen <- sub("\\.rds$", "", file_chosen)
# files_folder <- gsub("/", "_", files_folder)
# data_type <- paste0(files_folder, "_", file_chosen) # "simulation_data" or "real_data_{distance_used}"
# process <- "NGGPW" # Process type: "DP", "DPW", "NGGP", "NGGPW"
# method <- "LSS25+Gibbs1" # MCMC method used
# initialization <- "kmeans" # Initialization strategy
# filename <- paste0(data_type, "_", process, "_", method, "_", initialization, "_")
# save_with_name(folder, param, filename)

##############################################################################
# Visualization (Optional) ====
##############################################################################

plot_post_distr(mcmc_result, BI = param$BI)
plot_trace_cls(mcmc_result, BI = param$BI)
plot_post_sim_matrix(mcmc_result, BI = param$BI)
plot_trace_U(mcmc_result, BI = param$BI)
plot_acf_U(mcmc_result, BI = param$BI)
plot_cls_est(mcmc_result, BI = param$BI)
# plot_stats(mcmc_result, ground_truth, BI = param$BI)

puma_ids <- sf::st_read("input/CA/counties-pumas/counties-pumas.shp", quiet = TRUE)[["PUMA"]]
plot_map_prior_mean(unit_ids = puma_ids, puma_dir = "input/CA/counties-pumas")
plot_map_cls(
  results = mcmc_result,
  BI = param$BI,
  unit_ids = puma_ids,
  puma_dir = "input/CA/counties-pumas"
)

plot_hist_cls(
  results = mcmc_result,
  BI = param$BI,
  # input_dir = "input/CA/",
)
