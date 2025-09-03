library(Rcpp)
library(ggplot2)
library(MASS) # For fitdistr function
library(dplyr)
library(tidyr)
library(spam) # For comp.psm
library(fields) # For minVI
library(viridisLite) # For color scales
library(RColorBrewer) # For color palettes
library(pheatmap) # For heatmaps

# Function to set hyperparameters based on elbow method and k-means
plot_k_means <- function(dist_matrix, max_k = 10) {
  # Step 1: Compute elbow method for K selection
  wss <- numeric(max_k)
  # For distance matrices, we need to use multidimensional scaling to get coordinates
  mds_result <- cmdscale(dist_matrix, k = 2) # 2D embedding

  for (k in 1:max_k) {
    kmeans_result <- kmeans(mds_result, centers = k, nstart = 25)
    wss[k] <- kmeans_result$tot.withinss
  }

  # Plot the elbow curve
  elbow_data <- data.frame(K = 1:max_k, WSS = wss)
  elbow_plot <- ggplot(elbow_data, aes(x = K, y = WSS)) +
    geom_line() +
    geom_point() +
    labs(
      title = "Elbow Method for Optimal K",
      x = "Number of Clusters (K)",
      y = "Within-Cluster Sum of Squares"
    ) +
    theme_minimal()
  print(elbow_plot)
}

# Function to set hyperparameters using elbow method and k-means clustering
set_hyperparameters <- function(dist_matrix, k_elbow, plot_clustering = FALSE) {
  # For distance matrices, we need to use multidimensional scaling to get coordinates
  mds_result <- cmdscale(dist_matrix, k = 2) # 2D embedding

  # Step 2: Use k-means with K_elbow to get initial clustering
  kmeans_result <- kmeans(mds_result, centers = k_elbow, nstart = 25)
  initial_clusters <- kmeans_result$cluster

  # Step 3: Split pairwise distances into within-cluster (A) and inter-cluster (B)
  n <- nrow(dist_matrix)
  cat("Processing", n, "data points...\n")

  # Much more efficient vectorized approach
  # Create matrices of cluster labels
  cluster_i <- matrix(rep(initial_clusters, n), nrow = n, byrow = FALSE)
  cluster_j <- matrix(rep(initial_clusters, n), nrow = n, byrow = TRUE)

  # Find same cluster pairs (excluding diagonal)
  same_cluster <- (cluster_i == cluster_j) & (row(dist_matrix) < col(dist_matrix))
  diff_cluster <- (cluster_i != cluster_j) & (row(dist_matrix) < col(dist_matrix))

  # Extract distances
  within_cluster_distances <- dist_matrix[same_cluster]
  inter_cluster_distances <- dist_matrix[diff_cluster]

  cat("Vectorized processing completed!\n")

  cat("Within-cluster distances (A):", length(within_cluster_distances), "values\n")
  cat("Inter-cluster distances (B):", length(inter_cluster_distances), "values\n")

  # Plot the initial clustering from k-means (optional)
  if (plot_clustering) {
    cluster_data <- data.frame(
      x = mds_result[, 1],
      y = mds_result[, 2],
      cluster = as.factor(initial_clusters)
    )

    cluster_plot <- ggplot(cluster_data, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 3) +
      labs(
        title = paste("Initial K-means Clustering (K =", k_elbow, ")"),
        x = "MDS Dimension 1",
        y = "MDS Dimension 2"
      ) +
      theme_minimal()

    print(cluster_plot)
  }

  # Step 4: Fit Gamma distribution to within-cluster distances (A)
  # Using method of moments as initial estimates, then MLE
  if (length(within_cluster_distances) > 0) {
    # Remove zeros if any (distances to self)
    A <- within_cluster_distances[within_cluster_distances > 0]

    if (length(A) > 1) {
      # Fit gamma distribution using MLE with error handling
      tryCatch(
        {
          # Add small epsilon to avoid zero values that cause issues
          A_safe <- pmax(A, 1e-10)

          # Try method of moments first for better initial estimates
          mean_A <- mean(A_safe)
          var_A <- var(A_safe)

          # Method of moments estimates
          shape_est <- mean_A^2 / var_A
          rate_est <- mean_A / var_A

          # Use these as starting values for MLE
          gamma_fit_A <- fitdistr(A_safe, "gamma",
            start = list(shape = shape_est, rate = rate_est),
            lower = c(0.01, 0.01)
          )

          delta1 <- gamma_fit_A$estimate["shape"]

          # delta_1 should be less than 1 for shape parameter
          if (delta1 > 1) {
            delta1 <- 0.99 # Adjust shape parameter if needed
          }

          # Step 5: Set alpha and beta
          n_A <- length(A)
          alpha <- delta1 * n_A
          beta <- sum(A)

          cat("For within-cluster distances:\n")
          cat("  delta1 (shape) =", delta1, "\n")
          cat("  alpha =", alpha, "\n")
          cat("  beta =", beta, "\n")
        },
        error = function(e) {
          cat("Warning: Could not fit gamma distribution to within-cluster distances\n")
          cat("Error:", e$message, "\n")
          delta1 <<- 0.5
          alpha <<- 2
          beta <<- 2
        }
      )
    } else {
      # Fallback values
      delta1 <- 0.5
      alpha <- 2
      beta <- 2
      cat("Warning: Not enough within-cluster distances, using default values\n")
    }
  } else {
    # Fallback values
    delta1 <- 0.5
    alpha <- 2
    beta <- 2
    cat("Warning: No within-cluster distances found, using default values\n")
  }

  # Step 6: Repeat for inter-cluster distances (B)
  if (length(inter_cluster_distances) > 0) {
    B <- inter_cluster_distances[inter_cluster_distances > 0]

    if (length(B) > 1) {
      # Fit gamma distribution using MLE with error handling
      tryCatch(
        {
          # Add small epsilon to avoid zero values that cause issues
          B_safe <- pmax(B, 1e-10)

          # Try method of moments first for better initial estimates
          mean_B <- mean(B_safe)
          var_B <- var(B_safe)

          # Method of moments estimates
          shape_est <- mean_B^2 / var_B
          rate_est <- mean_B / var_B

          # Use these as starting values for MLE
          gamma_fit_B <- fitdistr(B_safe, "gamma",
            start = list(shape = shape_est, rate = rate_est),
            lower = c(0.01, 0.01)
          )

          delta2 <- gamma_fit_B$estimate["shape"]

          # delta_2 should be greater than 1 for shape parameter
          if (delta2 < 1) {
            delta2 <- 1.01 # Adjust shape parameter if needed
          }

          # Set zeta and gamma (note: gamma is a parameter name, using different variable)
          n_B <- length(B)
          zeta <- delta2 * n_B
          gamma_param <- sum(B) # renamed to avoid conflict with gamma() function

          cat("For inter-cluster distances:\n")
          cat("  delta2 (shape) =", delta2, "\n")
          cat("  zeta =", zeta, "\n")
          cat("  gamma =", gamma_param, "\n")
        },
        error = function(e) {
          cat("Warning: Could not fit gamma distribution to inter-cluster distances\n")
          cat("Error:", e$message, "\n")
          delta2 <<- 2
          zeta <<- 2
          gamma_param <<- 2
        }
      )
    } else {
      # Fallback values
      delta2 <- 2
      zeta <- 2
      gamma_param <- 2
      cat("Warning: Not enough inter-cluster distances, using default values\n")
    }
  } else {
    # Fallback values
    delta2 <- 2
    zeta <- 2
    gamma_param <- 2
    cat("Warning: No inter-cluster distances found, using default values\n")
  }

  # Return the computed hyperparameters
  return(list(
    delta1 = delta1,
    alpha = alpha,
    beta = beta,
    delta2 = delta2,
    zeta = zeta,
    gamma = gamma_param,
    k_elbow = k_elbow,
    initial_clusters = initial_clusters
  ))
}

save_with_name <- function(folder, params, initialization) {
  ## Name creation
  folder <- "results/"
  # Nomenclature: initialization + BI + NI + a + sigma + tau
  subfolder <- paste(initialization, "_BI", params$BI, "_NI", params$NI,
    "_a", params$a, "_sigma", params$sigma,
    "_tau", params$tau, "/",
    sep = ""
  )
  folder <- paste0(folder, subfolder)

  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }

  filename_results <- "simulation_results.rds"
  filename_gt <- "simulation_ground_truth.rds"
  filename_data <- "simulation_data.rds"
  filename_dist <- "simulation_distance_matrix.rds"
  filename_initial_params <- "simulation_initial_params.rds"

  filename_results <- paste0(folder, filename_results)
  filename_gt <- paste0(folder, filename_gt)
  filename_data <- paste0(folder, filename_data)
  filename_dist <- paste0(folder, filename_dist)
  filename_initial_params <- paste0(folder, filename_initial_params)

  saveRDS(results, file = filename_results)
  saveRDS(ground_truth, file = filename_gt)
  saveRDS(all_data, file = filename_data)
  saveRDS(dist_matrix, file = filename_dist)
  saveRDS(hyperparams, file = filename_initial_params)
}

# Generate data from 3 Gaussian distributions -> around 300 points
set.seed(42)
n_points <- 30
data1 <- matrix(rnorm(n_points * 2, mean = c(10, 10), sd = 1), ncol = 2)
data2 <- matrix(rnorm(n_points * 2, mean = c(3, 3), sd = 1), ncol = 2)
data3 <- matrix(rnorm(n_points * 2, mean = c(1, 1), sd = 1), ncol = 2)
data4 <- matrix(rnorm(n_points * 2, mean = c(6, 6), sd = 1), ncol = 2)

# Combine all data points + labels
ground_truth <- c(
  rep(0, n_points),
  rep(1, n_points),
  rep(2, n_points),
  rep(3, n_points)
)
all_data <- rbind(data1, data2, data3, data4)

# Create distance matrix (n x n)
dist_matrix <- as.matrix(dist(all_data))

# Create cluster labels for plotting
cluster_labels <- c(
  rep("Cluster 1", n_points),
  rep("Cluster 2", n_points),
  rep("Cluster 3", n_points),
  rep("Cluster 4", n_points)
)

# Load the C++ code
sourceCpp("src/launcher.cpp")

# Set hyperparameters using the NORMALIZED distances
#plot_k_means(dist_matrix, max_k = 8)
hyperparams <- set_hyperparameters(dist_matrix, k_elbow = 4)

cat("Final hyperparameters:\n")
cat("delta1:", hyperparams$delta1, "\n")
cat("alpha:", hyperparams$alpha, "\n")
cat("beta:", hyperparams$beta, "\n")
cat("delta2:", hyperparams$delta2, "\n")
cat("zeta:", hyperparams$zeta, "\n")
cat("gamma:", hyperparams$gamma, "\n")

# Create parameter object with computed hyperparameters
param <- new(
  Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
  200, 1000, 4, 1.0, 1.0 # Increased a from 4 to 20
) # BI, NI, a, sigma, tau

# Initialize allocations

## All-in-one allocation
#hyperparams$initial_clusters <- rep(0, nrow(dist_matrix)) # All points in one cluster

## Sequential allocation
hyperparams$initial_clusters <- seq(0, nrow(dist_matrix) - 1)

## k-means allocation different from the one used for hyperparameters
# hyperparams$initial_clusters <- kmeans(all_data,
#   centers = 4,
#   nstart = 25
# )$cluster - 1

## k-means allocation used for hyperparameters
# hyperparams$initial_clusters <- hyperparams$initial_clusters - 1

print(table(hyperparams$initial_clusters))

# Run MCMC with computed hyperparameters
results <- mcmc(dist_matrix, param, hyperparams$initial_clusters)

## Save into files
# save_with_name(folder, param, "allInOne")
