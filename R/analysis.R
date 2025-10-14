source("R/utils.R")

list.files("results/")

# Load the results and ground truth for later analysis
folder <- "results/DP_splitmerge50_Neal1_kmeans_0.25sigma_10d_BI2000_NI10000_a1_sigma0.4_tau1/"
filename_results <- "simulation_results.rds"
filename_gt <- "simulation_ground_truth.rds"
filename_dist <- "simulation_distance_matrix.rds"
filename_data <- "simulation_data.rds"
filename_initial_params <- "simulation_initial_params.rds"

filename_results <- paste0(folder, filename_results)
filename_gt <- paste0(folder, filename_gt)
filename_dist <- paste0(folder, filename_dist)
filename_data <- paste0(folder, filename_data)
filename_initial_params <- paste0(folder, filename_initial_params)

results <- readRDS(file = filename_results)
ground_truth <- readRDS(file = filename_gt)
dist_matrix <- readRDS(file = filename_dist)
all_data <- readRDS(file = filename_data)
param <- readRDS(file = filename_initial_params)

# Create plot directory if it doesn't exist
folder <- paste0(folder, "plots/")
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}

#################################################################
############################## Plot Data ########################
#################################################################

#plot_data(all_data, ground_truth, save = TRUE, folder)
distance_plot(all_data, ground_truth, save = TRUE, folder)

#################################################################
################ Plot Simulation Results ########################
#################################################################

cat("\nFinal hyperparameters:\n")
cat("delta1 =", param$delta1, "\n")
cat("alpha =", param$alpha, "\n")
cat("beta =", param$beta, "\n")
cat("delta2 =", param$delta2, "\n")
cat("gamma =", param$gamma, "\n")
cat("zeta =", param$zeta, "\n")
cat("Initial clusters =", param$initial_cluster, "\n")

plot_mcmc_results(results, ground_truth, BI = 2000, save = TRUE, folder = folder)
