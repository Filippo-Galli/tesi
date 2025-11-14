source("R/utils_plot.R")

files <- list.files("results/")
files

##############################################################################
# Load Results ====
##############################################################################
folder <- paste0("results/", files[10], "/")
filename_results <- "simulation_results.rds"
#filename_gt <- "simulation_ground_truth.rds"
filename_dist <- "simulation_distance_matrix.rds"
#filename_data <- "simulation_data.rds"
filename_initial_params <- "simulation_initial_params.rds"

filename_results <- paste0(folder, filename_results)
#filename_gt <- paste0(folder, filename_gt)
filename_dist <- paste0(folder, filename_dist)
#filename_data <- paste0(folder, filename_data)
filename_initial_params <- paste0(folder, filename_initial_params)

results <- readRDS(file = filename_results)
#ground_truth <- readRDS(file = filename_gt)
dist_matrix <- readRDS(file = filename_dist)
#all_data <- readRDS(file = filename_data)
param <- readRDS(file = filename_initial_params)

# Create plot directory if it doesn't exist
folder <- paste0(folder, "VI_plots/")
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}

#################################################################
# Plot Data ====
#################################################################

#plot_data(all_data, ground_truth, save = FALSE, folder) # Uncomment to plot data
#plot_distance(dist_matrix, save = FALSE, folder) # Plot distance matrix

#################################################################
# Plot Simulation Results ====
#################################################################

cat("\nFinal hyperparameters:\n")
cat("delta1 =", param$delta1, "\n")
cat("alpha =", param$alpha, "\n")
cat("beta =", param$beta, "\n")
cat("delta2 =", param$delta2, "\n")
cat("gamma =", param$gamma, "\n")
cat("zeta =", param$zeta, "\n")
cat("Initial clusters =", param$initial_cluster, "\n")

##############################################################################
# Visualization ====
##############################################################################

plot_post_distr(results, BI = 10000, save = TRUE, folder = folder)
plot_trace_cls(results, BI = 10000, save = TRUE, folder = folder)
plot_post_sim_matrix(results, BI = 10000, save = TRUE, folder = folder)
plot_trace_U(results, BI = 10000, save = TRUE, folder = folder)
plot_acf_U(results, BI = 10000, save = TRUE, folder = folder)

puma_ids <- sf::st_read("input/West_Midwest/counties-pumas/counties-pumas.shp", quiet = TRUE)[["PUMA"]]
#plot_map_prior_mean(unit_ids = puma_ids, puma_dir = "input/West_Midwest/counties-pumas", save = TRUE, folder = folder)
plot_map_cls(
  results = results,
  BI = 10000,
  unit_ids = puma_ids, 
  puma_dir = "input/West_Midwest/counties-pumas", 
  save = TRUE, folder = folder
)

plot_hist_cls(
  results = results,
  BI = 10000, 
  save = TRUE, folder = folder,
  input_dir = "input/West_Midwest/",
)
