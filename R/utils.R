## @file utils.R
## @brief Utility functions for Bayesian clustering analysis and MCMC visualization
## @author Filippo Galli
## @date 2025

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
library(mvtnorm)
library(gtools)
library(salso)
library(aricode)
library(reshape2)
library(cluster)

retrieve_W <- function(distance_matrix, neighbours = 8) {
  # For each element find the nearest neighbours
  n <- nrow(distance_matrix)
  W <- matrix(0, n, n) # Adjacency matrix
  for (i in 1:n) {
    # Get indices of the nearest neighbours (excluding self)
    nn_indices <- order(distance_matrix[i, ])[2:(neighbours + 1)]
    W[i, nn_indices] <- 1
  }

  # Make W only upper triangular to avoid double counting
  W <- (W + t(W)) > 0
  # W[lower.tri(W)] <- 0

  return(W)
}

generate_mixture_data <- function(N = 100, K = 10, alpha = 10, dim = K, radius = 1, sigma = 0.1, ordered = TRUE) {
  # Input validation
  if (N < 1) stop("N must be greater than 1")
  if (K < 1 || K > N) stop("K must satisfy 1 ≤ K ≤ N")
  if (alpha <= 0) stop("alpha must be positive")
  if (dim < K) stop("dim must be ≥ K")
  if (radius <= 0) stop("radius must be positive")
  if (sigma <= 0) stop("sigma must be positive")

  # Generate cluster weights from Dirichlet prior
  probs <- as.numeric(rdirichlet(1, rep(alpha, K)))

  # Generate cluster assignments
  clusts <- sample(1:K, N, replace = TRUE, prob = probs)

  # Generate cluster centres as vertices of dim-dimensional simplex
  # Center i has radius at position i, zeros elsewhere
  clust_centres <- matrix(0, K, dim)
  for (i in 1:K) {
    clust_centres[i, i] <- radius
  }

  # Covariance matrix: sigma^2 * I_dim
  Sigma <- diag(sigma^2, dim)

  # Generate points
  points <- matrix(0, N, dim)
  for (i in 1:N) {
    cluster <- clusts[i]
    points[i, ] <- rmvnorm(1, mean = clust_centres[cluster, ], sigma = Sigma)
  }

  # Optional: Order data by cluster assignments
  if (ordered) {
    # Sort by cluster assignment
    order_idx <- order(clusts)
    points <- points[order_idx, ]
    clusts <- clusts[order_idx]

    # Also return the ordering for reference
    return(list(
      points = points,
      clusts = clusts,
      clust_centres = clust_centres,
      probs = probs,
      original_order = order_idx # In case you need to map back
    ))
  } else {
    return(list(
      points = points,
      clusts = clusts,
      clust_centres = clust_centres,
      probs = probs
    ))
  }
}

save_with_name <- function(folder, params, initialization, gt = FALSE) {
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
  if (gt) {
    filename_results <- "simulation_ground_truth.rds"
  }
  filename_dist <- "simulation_distance_matrix.rds"
  filename_initial_params <- "simulation_initial_params.rds"

  filename_results <- paste0(folder, filename_results)
  if (gt) {
    filename_gt <- paste0(folder, filename_gt)
  }
  filename_dist <- paste0(folder, filename_dist)
  filename_initial_params <- paste0(folder, filename_initial_params)

  saveRDS(mcmc_result, file = filename_results)
  if (gt) {
    saveRDS(ground_truth, file = filename_gt)
  }
  saveRDS(dist_matrix, file = filename_dist)
  saveRDS(hyperparams, file = filename_initial_params)
}

set_hyperparameters <- function(dist_matrix, k_elbow, ground_truth = NULL, plot_clustering = FALSE, plot_distribution = TRUE) {

  # Use k-medoids to get initial clusters
  pam_fit <- pam(dist_matrix, k = k_elbow, diss = TRUE)
  initial_clusters <- pam_fit$clustering

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
      x = coords[, 1],
      y = coords[, 2],
      cluster = as.factor(initial_clusters)
    )

    # Add ground truth if provided
    if (!is.null(ground_truth)) {
      cluster_data$ground_truth <- as.factor(ground_truth)
    }

    cluster_plot <- ggplot(cluster_data, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 3) +
      labs(
        title = paste("Initial K-means Clustering (K =", k_elbow, ")"),
        x = "X Coordinate",
        y = "Y Coordinate"
      ) +
      theme_minimal()

    print(cluster_plot)

    # Plot ground truth if provided
    if (!is.null(ground_truth)) {
      gt_plot <- ggplot(cluster_data, aes(x = x, y = y, color = ground_truth)) +
        geom_point(size = 3) +
        labs(
          title = "Ground Truth Clustering",
          x = "X Coordinate",
          y = "Y Coordinate"
        ) +
        theme_minimal()

      print(gt_plot)
    }
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

          # delta_1 should be less than 1 for cohesion
          if (delta1 > 1) {
            print(paste("Warning: Fitted delta1 > 1, adjusting to 0.9 for cohesion. Old value: ", delta1))
            delta1 <- 0.9 # Adjust shape parameter if needed
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

          # delta_2 should be greater than 1 for repulsion
          if (delta2 < 1) {
            print(paste("Warning: Fitted delta2 < 1, adjusting to 1.5 for repulsion. Old value: ", delta2))
            delta2 <- 1.5 # Adjust shape parameter if needed
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

  if (plot_distribution == TRUE) {
    # First plot: Within-cluster distances with algorithm's gamma parameters
    if (exists("A") && length(A) > 1) {
      A_df <- data.frame(Distance = A)

      # Use the algorithm's parameters: delta1, alpha, beta
      # Convert from alpha, beta to shape, rate parameterization
      # Your algorithm uses: alpha = delta1 * n_A, beta = sum(A)
      # For plotting, we need the rate parameter that corresponds to these

      # The relationship is: alpha = shape * rate, beta = sum(data)
      # So: rate = alpha / delta1 = n_A (when alpha = delta1 * n_A)
      # But we want E[X] = shape/rate to match the empirical mean
      # Better approach: use delta1 as shape, and derive rate from the mean

      algorithm_shape_A <- delta1
      # For gamma distribution: mean = shape/rate, so rate = shape/mean
      empirical_mean_A <- mean(A)
      algorithm_rate_A <- algorithm_shape_A / empirical_mean_A

      p1 <- ggplot(A_df, aes(x = Distance)) +
        geom_histogram(aes(y = after_stat(density)),
          bins = 30, fill = "lightblue",
          color = "black", alpha = 0.7
        ) +
        stat_function(
          fun = dgamma,
          args = list(shape = algorithm_shape_A, rate = algorithm_rate_A),
          color = "red", linewidth = 1
        ) +
        labs(
          title = paste(
            "Within-Cluster Distances with Algorithm's Gamma\n",
            "δ₁ (shape) =", round(algorithm_shape_A, 3),
            ", α =", round(alpha, 3), ", β =", round(beta, 3)
          ),
          x = "Distance", y = "Density"
        ) +
        theme_minimal()
      print(p1)
    }

    # Second plot: Inter-cluster distances with algorithm's gamma parameters
    if (exists("B") && length(B) > 1) {
      B_df <- data.frame(Distance = B)

      # Use the algorithm's parameters: delta2, zeta, gamma_param
      algorithm_shape_B <- delta2
      # For gamma distribution: mean = shape/rate, so rate = shape/mean
      empirical_mean_B <- mean(B)
      algorithm_rate_B <- algorithm_shape_B / empirical_mean_B

      p2 <- ggplot(B_df, aes(x = Distance)) +
        geom_histogram(aes(y = after_stat(density)),
          bins = 30, fill = "lightgreen",
          color = "black", alpha = 0.7
        ) +
        stat_function(
          fun = dgamma,
          args = list(shape = algorithm_shape_B, rate = algorithm_rate_B),
          color = "blue", linewidth = 1
        ) +
        labs(
          title = paste(
            "Inter-Cluster Distances with Algorithm's Gamma\n",
            "δ₂ (shape) =", round(algorithm_shape_B, 3),
            ", ζ =", round(zeta, 3), ", γ =", round(gamma_param, 3)
          ),
          x = "Distance", y = "Density"
        ) +
        theme_minimal()
      print(p2)
    }
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
    initial_clusters = initial_clusters,
    ground_truth = ground_truth
  ))
}

compute_hist_distances <- function(hist1, hist2, type = "Histogram-Divergence") {
  # Extract counts from histogram objects
  counts1 <- hist1$counts
  counts2 <- hist2$counts

  # Get bin widths (important for proper CM and Wasserstein)
  breaks1 <- hist1$breaks
  breaks2 <- hist2$breaks

  # Calculate bin width (assuming equal widths)
  bin_width1 <- diff(breaks1)[1]
  bin_width2 <- diff(breaks2)[1]
  if (abs(bin_width1 - bin_width2) > 1e-10) {
    warning("Histograms have different bin widths. Results may be inaccurate.")
  }
  bin_width <- bin_width1  # Use first histogram's bin width

  # Handle different bin numbers by aligning to common bins
  max_bins <- max(length(counts1), length(counts2))
  
  # Pad shorter histogram with zeros if needed
  if (length(counts1) < max_bins) {
    counts1 <- c(counts1, rep(0, max_bins - length(counts1)))
  }
  if (length(counts2) < max_bins) {
    counts2 <- c(counts2, rep(0, max_bins - length(counts2)))
  }
  
  # Normalize to probability distributions (sum to 1)
  counts1 <- counts1 / sum(counts1)
  counts2 <- counts2 / sum(counts2)

  if (type == "Histogram-Divergence") {
    # Compute Histogram Divergence (intersection)
    divergence <- sum(pmin(counts1, counts2))
    return(divergence)
  } else if (type == "Jeff") {
    # Compute Jeffrey Divergence
    temp_1 <- counts1 > 0
    temp_2 <- counts2 > 0
    temp_3 <- temp_1 & temp_2
    kl1 <- sum(counts1[temp_3] * log(counts1[temp_3] / counts2[temp_3]))
    kl2 <- sum(counts2[temp_3] * log(counts2[temp_3] / counts1[temp_3]))

    jeffrey_divergence <- kl1 + kl2
    return(jeffrey_divergence)
  } else if (type == "chi2") {
    # Compute Chi-squared distance
    temp_2 <- counts2 > 0
    chi2_distance <- ((counts1[temp_2] - counts2[temp_2])^2) / counts2[temp_2]
    return(sum(chi2_distance))
  } else if (type == "euclidean") {
    # Compute Euclidean distance
    euclidian_distance <- sqrt(sum((counts1 - counts2)^2))
    return(euclidian_distance)
  } else if (type == "CM") {
    # Cramér-von Mises: weight by bin width
    cum1 <- cumsum(counts1)
    cum2 <- cumsum(counts2)
    cm_distance <- bin_width * sum((cum1 - cum2)^2)
    return(cm_distance)
  } else if (type == "Wasserstein") {
    # Wasserstein: weight by bin width
    cum1 <- cumsum(counts1)
    cum2 <- cumsum(counts2)
    wasserstein_distance <- bin_width * sum(abs(cum1 - cum2))
    return(wasserstein_distance)
  } else {
    stop("Unsupported distance type")
  }
}
