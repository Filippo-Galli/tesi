source("R/utils.R")
source("R/utils_plot.R")
source("R/mcmc_loop.R")

## Set random seed for reproducibility
set.seed(44)

##############################################################################
# Data Loading ====
##############################################################################

## Load simulated data
# sigma <- .2
# d <- 50
# file_chosen <- paste0("Natarajan_", sigma, "sigma_", d, "d")
# files_folder <- paste0("simulation_data/", file_chosen)
# all_data <- readRDS(file = paste0(files_folder, "/all_data.rds"))
# ground_truth <- readRDS(file = paste0(files_folder, "/ground_truth.rds"))
# dist_matrix <- readRDS(file = paste0(files_folder, "/dist_matrix.rds"))

## Load real data
files_folder <- "real_data/LA"
files <- list.files(files_folder)
file_chosen <- files[3]
dist_matrix <- readRDS(file = paste0(files_folder, "/", file_chosen))
puma_age <- readRDS(file = paste0(files_folder, "/puma_agep_std_mean.rds"))
# plot_distance(dist_matrix)

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
# sourceCpp("src/launcher.cpp", rebuild = TRUE, cacheDir = "~/my_rcpp_cache") # useful for perf
sourceCpp("src/bindings.cpp")
cat("✅ C++ code compiled successfully!\n\n")

##############################################################################
# Hyperparameter Configuration ====
##############################################################################

# Plot k-means elbow method to help set hyperparameters
# plot_k_means(dist_matrix, max_k = 10)

# Set hyperparameters based on distance matrix and save it for future use
# hyperparams <- set_hyperparameters(dist_matrix,
#   k_elbow = 3, plot_clustering = FALSE, plot_distribution = FALSE
# )
# saveRDS(hyperparams, file = paste0(files_folder, "/hyperparameters_", sub("\\.rds$", "", file_chosen), ".rds"))
hyperparams <- readRDS(file = paste0(files_folder, "/hyperparameters_", sub("\\.rds$", "", file_chosen), ".rds"))

##############################################################################
# Parameter Object Initialization ====
##############################################################################

# Create Params using factory function instead of module constructor
param <- create_Params(
    hyperparams$delta1,
    hyperparams$alpha,
    hyperparams$beta,
    hyperparams$delta2,
    hyperparams$gamma,
    hyperparams$zeta,
    10000, # BI
    10000, # NI
    1, # a
    0.1, # sigma
    1, # tau
    dist_matrix # D (distance matrix)
)

##############################################################################
# Covariates Object Initialization ====
##############################################################################

B <- 10 * var(puma_age$Mean_AGEP_std) # prior variance
m <- 0 # prior mean
v <- 0.5 * var(puma_age$Mean_AGEP_std) # known variance

# Ensure W is integer matrix
W <- matrix(as.integer(W), nrow = nrow(W), ncol = ncol(W))

# Create Covariates using factory function instead of module constructor
covariates <- create_Covariates(
    W, # Spatial adjacency matrix (must be integer)
    1, # spatial_coefficient
    as.integer(puma_age$Mean_AGEP_std), # ages vector (must be integer)
    B, m, v, # covariate prior parameters
    TRUE, 1, 1 # fixed_v, nu, S0
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

mcmc_result <- run_mcmc(param, covariates, hyperparams$initial_clusters)
elapsed_time <- mcmc_result$elapsed_time

##############################################################################
# Save Results (Optional) ====
##############################################################################
file_chosen <- sub("\\.rds$", "", file_chosen)
files_folder_clean <- gsub("/", "_", files_folder)
data_type <- paste0(files_folder_clean, "_", sub("^distance_", "", file_chosen))
process <- "TEST-FIXEDV-NGGPWx" # Process type: "DP", "DPW", "NGGP", "NGGPW", NGGPWx
method <- "LSS_SDDS25+Gibbs1" # MCMC method used
initialization <- "kmeans" # Initialization strategy
filename <- paste0(data_type, "_", process, "_", method, "_", initialization, "_")
save_with_name("results/", param, filename)

##############################################################################
# Visualization (Optional) ====
##############################################################################

# plot_post_distr(mcmc_result, BI = mcmc_result$BI)
# plot_trace_cls(mcmc_result, BI = mcmc_result$BI)
# plot_post_sim_matrix(mcmc_result, BI = mcmc_result$BI)
# plot_trace_U(mcmc_result, BI = mcmc_result$BI)
# plot_acf_U(mcmc_result, BI = mcmc_result$BI)
# plot_cls_est(mcmc_result, BI = mcmc_result$BI)
# plot_stats(mcmc_result, ground_truth, BI = mcmc_result$BI)

# puma_ids <- sf::st_read("input/LA/counties-pumas/counties-pumas.shp", quiet = TRUE)[["PUMA"]]
# plot_map_prior_mean(unit_ids = puma_ids, puma_dir = "input/CA/counties-pumas")
# plot_map_cls(
#   results = mcmc_result,
#   BI = mcmc_result$BI,
#   unit_ids = puma_ids,
#   puma_dir = "input/LA/counties-pumas"
# )

# plot_hist_cls(
#   results = mcmc_result,
#   BI = mcmc_result$BI,
#   # input_dir = "input/CA/",
# )
