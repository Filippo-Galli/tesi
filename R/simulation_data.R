# Import utils.R for additional functions
source("R/utils.R")

# Generate data from 4 Gaussian distributions
set.seed(44)

### ----------------------- Natarajan Data Generation ----------------- ###
# data_generation <- generate_mixture_data(N = 100, sigma = 0.25, d = 10)
# data_generation <- generate_mixture_data(N = 100, sigma = 0.2, d = 50)
data_generation <- generate_mixture_data(N = 100, sigma = 0.18, d = 10)

all_data <- data_generation$points
ground_truth <- data_generation$clusts

# Distance plot
# distance_plot(all_data, ground_truth)

### ----------------------- Gaussian Data Generation ----------------- ###

# n_points <- 30
# centers <- list(
#   c(10, 10),  # Cluster 1: top-right
#   c(6, 6),  # Cluster 2: bottom-right
#   c(0, 0),  # Cluster 3: bottom-left
#   c(3, 3)   # Cluster 4: top-left
# )

# # Generate data for each cluster
# data1 <- mvrnorm(n_points, mu = centers[[1]], Sigma = 3*diag(2))
# data2 <- mvrnorm(n_points, mu = centers[[2]], Sigma = 3*diag(2))
# data3 <- mvrnorm(n_points, mu = centers[[3]], Sigma = 3*diag(2))
# data4 <- mvrnorm(n_points, mu = centers[[4]], Sigma = 3*diag(2))

# # Combine all data points + labels
# ground_truth <- c(
#   rep(0, n_points),
#   rep(1, n_points),
#   rep(2, n_points),
#   rep(3, n_points)
# )
# all_data <- rbind(data1, data2, data3, data4)

### ----------------------- Gamma Data Generation ----------------- ###
# n_points <- 30
# # Gamma Generated data
# data1 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) + 8, ncol = 2)
# data2 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) + 1, ncol = 2)
# data3 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) - 1, ncol = 2)
# data4 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) + 4, ncol = 2)

# # Combine all data points + labels
# ground_truth <- c(
#   rep(0, n_points),
#   rep(1, n_points),
#   rep(2, n_points),
#   rep(3, n_points)
# )
# all_data <- rbind(data1, data2, data3, data4)

# Create distance matrix (n x n)
dist_matrix <- as.matrix(dist(all_data))
# Normalize distance to be between 2 and 4
lower_bound <- 0
upper_bound <- 10000
# dist_matrix <- lower_bound + (upper_bound - lower_bound)*(dist_matrix - min(dist_matrix)) / (max(dist_matrix) - min(dist_matrix))
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
  1000, 5000, 1, 1.0, 1.0, 1, W
) # BI, NI, a, sigma, tau, coeff spatial dep, W

# Initialize allocations

## All-in-one allocation
# hyperparams$initial_clusters <- rep(0, nrow(dist_matrix)) # All points in one cluster

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
# save_with_name(folder, param, "kmeans_Natarajan_02-50")

plot_mcmc_results(mcmc_result, as.factor(ground_truth), BI = param$BI)
