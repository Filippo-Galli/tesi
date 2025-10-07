# Import utils.R for additional functions
source("R/utils.R")

set.seed(44)

# Read simulated data
folder <- "simulation_data/Natarajan_0.2sigma_50d"

all_data <- readRDS(file = paste0(folder, "/all_data.rds"))
ground_truth <- readRDS(file = paste0(folder, "/ground_truth.rds"))
dist_matrix <- readRDS(file = paste0(folder, "/dist_matrix.rds"))

# Retrieve adjacency matrix W from distance matrix
W <- retrieve_W(dist_matrix)

# Load the C++ code
sourceCpp("src/launcher.cpp")

# Set hyperparameters with ground truth and plotting
# plot_k_means(dist_matrix, max_k = 15)
hyperparams <- set_hyperparameters(all_data, dist_matrix,
  k_elbow = 3,
  ground_truth = ground_truth, plot_clustering = FALSE, plot_distribution = FALSE
)

# Create parameter object with computed hyperparameters
param <- new(
  Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
  2000, 10000, 1, # BI, NI, a,
  0.00001, 1, 1, # sigma, tau, coeff spatial dep
  W # Spatial adjacency matrix
)

# Initialize allocations

## All-in-one allocation
#hyperparams$initial_clusters <- rep(0, nrow(dist_matrix)) # All points in one cluster

## Sequential allocation
# hyperparams$initial_clusters <- seq(0, nrow(dist_matrix) - 1)

## k-means allocation different from the one used for hyperparameters
# hyperparams$initial_clusters <- kmeans(all_data,
#   centers = 2,
#   nstart = 25
# )$cluster - 1

## k-means allocation used for hyperparameters
hyperparams$initial_clusters <- hyperparams$initial_clusters - 1

print("Initial cluster allocation:")
print(table(hyperparams$initial_clusters))

# Run MCMC with computed hyperparameters
log_file <- "mcmc_log.txt"
if (file.exists(log_file)) {
  file.remove(log_file) # remove previous log file
}

results <- capture.output(
  {
    mcmc_result <- mcmc(dist_matrix, param, hyperparams$initial_clusters)
  },
  file = log_file
)

## Save into files - initialization name
#save_with_name(folder, param, "DP_Neal2W1_SMW1_kmeans_025_10")

plot_mcmc_results(mcmc_result, as.factor(ground_truth), BI = param$BI)

# Plot U over iterations
plot(mcmc_result$U, type = 'l', main = "U over iterations", ylab = "U", xlab = "Iteration")
abline(h = mean(mcmc_result$U), col = "red", lty = 2)
legend("topright", legend = paste("Mean U =", round(mean(mcmc_result$U), 3)), col = "red", lty = 2)
