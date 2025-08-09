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
library(mcclust.ext) # For MCMC clustering functions


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
          gamma_fit_A <- fitdistr(A, "gamma")
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
          gamma_fit_B <- fitdistr(B, "gamma")
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

# Function to plot MCMC analysis results
plot_mcmc_results <- function(results, true_labels) {
  ### First plot - Posterior distribution of the number of clusters
  post_k <- table(unlist(results$K)) / length(unlist(results$K))
  df <- data.frame(
    cluster_found = as.numeric(names(post_k)),
    rel_freq = as.numeric(post_k)
  )

  p1 <- ggplot(data = df, aes(x = factor(cluster_found), y = rel_freq)) +
    geom_col() +
    labs(
      x = "Cluster Found",
      y = "Relative Frequency",
    ) +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      text = element_text(size = 15),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_line(color = "grey95")
    ) +
    scale_x_discrete(drop = FALSE) # Ensures all cluster_found values are shown
  print(p1)

  ### Second plot - Trace of number of clusters
  # Ensure K values exist and handle any missing values
  k_values <- unlist(results$K)
  # Remove any NULL or NA values and create corresponding iteration indices
  valid_indices <- which(!is.null(results$K) & !is.na(unlist(results$K)))
  k_values <- k_values[!is.na(k_values)]

  k_df <- data.frame(
    Iteration = valid_indices[seq_along(k_values)],
    NumClusters = k_values
  )

  k_df_long <- k_df %>%
    pivot_longer(
      cols = starts_with("NumClusters"),
      names_to = "variable",
      values_to = "value"
    )

  p2 <- ggplot(k_df_long, aes(x = Iteration, y = value)) +
    geom_line() +
    labs(
      x = "Iteration",
      y = "Number of clusters",
    ) +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      text = element_text(size = 15),
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_line(color = "grey95")
    )
  print(p2)

  ### Third plot - Posterior Similarity Matrix
  # Compute posterior similarity matrix
  n <- nrow(dist_matrix)
  similarity_matrix <- matrix(0, nrow = n, ncol = n)

  # For each MCMC iteration, add to similarity matrix
  for (iter in seq_along(results$allocations)) {
    allocation <- results$allocations[[iter]]
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (allocation[i] == allocation[j]) {
          similarity_matrix[i, j] <- similarity_matrix[i, j] + 1
          similarity_matrix[j, i] <- similarity_matrix[j, i] + 1
        }
      }
    }
  }

  # Normalize by number of iterations
  similarity_matrix <- similarity_matrix / length(results$allocations)

  # Set diagonal to 1 (each point is always similar to itself)
  diag(similarity_matrix) <- 1

  # Create heatmap of posterior similarity matrix
  similarity_df <- expand.grid(i = 1:n, j = 1:n)
  similarity_df$similarity <- as.vector(similarity_matrix)

  p3 <- ggplot(similarity_df, aes(x = i, y = j, fill = similarity)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkblue") +
    labs(
      x = "Data Point Index",
      y = "Data Point Index",
      fill = "Posterior\nSimilarity",
      title = "Posterior Similarity Matrix"
    ) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      text = element_text(size = 12),
      panel.background = element_blank()
    ) +
    coord_fixed()
  print(p3)

  ### Fourth plot - Posterior similarity matrix
  # Check if allocations exist, otherwise skip this plot
  if (!is.null(results$allocations) && length(results$allocations) > 0) {
    C <- matrix(unlist(lapply(results$allocations, function(x) x + 1)),
      nrow = length(results$allocations),
      ncol = length(true_labels),
      byrow = TRUE
    )

    required_packages <- c("spam", "fields", "viridisLite", "RColorBrewer", "pheatmap", "mcclust")
    for (pkg in required_packages) {
      if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
        library(pkg, character.only = TRUE)
      }
    }

    psm <- comp.psm(C)
    VI <- minVI(psm)

    cat("Cluster Sizes:\n")
    print(table(VI$cl))
    cat("\nAdjusted Rand Index:", arandi(VI$cl, true_labels), "\n")

  } else {
    cat("Warning: No allocation data found, skipping posterior similarity matrix plot\n")
  }

  ### Fifth plot - Auto-correlation plot
  # Check if loglikelihood data exists and is valid
  if (!is.null(results$loglikelihood) && length(results$loglikelihood) > 0) {
    # Remove any NA or infinite values
    logl_clean <- results$loglikelihood[is.finite(results$loglikelihood)]
    k_clean <- unlist(results$K)[is.finite(unlist(results$K))]

    if (length(logl_clean) > 1 && length(k_clean) > 1) {
      # Ensure both vectors have the same length
      min_length <- min(length(logl_clean), length(k_clean))
      mcmc_list <- list(
        ncls = k_clean[1:min_length],
        logl = logl_clean[1:min_length]
      )
      mcmc_matrix <- do.call(cbind, mcmc_list)

      # Only plot if we have valid data
      if (all(is.finite(mcmc_matrix))) {
        acf(mcmc_matrix, main = "Autocorrelation of MCMC chains")
      } else {
        cat("Warning: Skipping autocorrelation plot due to invalid data\n")
      }
    } else {
      cat("Warning: Not enough valid data for autocorrelation plot\n")
    }
  } else {
    cat("Warning: No loglikelihood data available for autocorrelation plot\n")
  }
}

# Generate data from 3 Gaussian distributions
set.seed(123)
n_points <- 100
data1 <- matrix(rnorm(n_points * 2, mean = c(10, 10), sd = 1), ncol = 2)
data2 <- matrix(rnorm(n_points * 2, mean = c(1, 1), sd = 1), ncol = 2)
data3 <- matrix(rnorm(n_points * 2, mean = c(5, 5), sd = 2), ncol = 2)

# Combine all data points + labels
ground_truth <- c(
  rep(0, n_points),
  rep(1, n_points),
  rep(2, n_points)
)
all_data <- rbind(data1, data2, data3)

# Create distance matrix (n x n)
dist_matrix <- as.matrix(dist(all_data))

# Create cluster labels for plotting
cluster_labels <- c(
  rep("Cluster 1", n_points),
  rep("Cluster 2", n_points),
  rep("Cluster 3", n_points)
)

# Create data frame for plotting
plot_data <- data.frame(
  x = all_data[, 1],
  y = all_data[, 2],
  cluster = cluster_labels
)

# Plot the clusters
ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 3) +
  labs(
    title = "Three Gaussian Clusters",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_minimal()

# Load the C++ code
sourceCpp("src/launcher.cpp")

# Set hyperparameters using the new method
plot_k_means(dist_matrix, max_k = 10)
hyperparams <- set_hyperparameters(dist_matrix, k_elbow = 3)

# Create parameter object with computed hyperparameters
param <- new(
  Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
  1000, 2000, 3, 1.0, 1.0
) # BI, NI, a, sigma, tau

cat("\nFinal hyperparameters:\n")
cat("delta1 =", param$delta1, "\n")
cat("alpha =", param$alpha, "\n")
cat("beta =", param$beta, "\n")
cat("delta2 =", param$delta2, "\n")
cat("gamma =", param$gamma, "\n")
cat("zeta =", param$zeta, "\n")

# Run MCMC with computed hyperparameters
results <- mcmc(dist_matrix, param, first_allocation = "sequential")

## Check results
plot_mcmc_results(results, ground_truth)
