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


retrieve_W <- function(distance_matrix, neighbours = 8) {
  # For each element find the nearest neighbours
  n <- nrow(distance_matrix)
  W <- matrix(0, n, n) # Adjacency matrix
  for (i in 1:n) {
    # Get indices of the nearest neighbours (excluding self)
    nn_indices <- order(distance_matrix[i, ])[2:(neighbours + 1)]
    W[i, nn_indices] <- 1
  }
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
      original_order = order_idx  # In case you need to map back
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

distance_plot <- function(all_data, clusts, save = FALSE, folder = "results/plots/") {
  # Calculate distance matrix
  dist_matrix <- as.matrix(dist(all_data))

  # Get upper triangular indices to avoid duplicates
  upper_tri_indices <- which(upper.tri(dist_matrix), arr.ind = TRUE)

  # Extract distances and determine if they are intra or inter-cluster
  distances <- dist_matrix[upper_tri_indices]
  cluster_pairs <- cbind(clusts[upper_tri_indices[,1]], clusts[upper_tri_indices[,2]])

  # Classify distances as intra-cluster (same cluster) or inter-cluster (different clusters)
  intra_cluster <- distances[cluster_pairs[,1] == cluster_pairs[,2]]
  inter_cluster <- distances[cluster_pairs[,1] != cluster_pairs[,2]]

  # Create histogram with overlaid distributions
  hist(inter_cluster, breaks = 30, col = rgb(0, 0, 1, 0.7), # Orange with transparency
      main = "Histogram of Pairwise Distances", 
      xlab = "Distance", 
      ylab = "Frequency",
      xlim = range(c(intra_cluster, inter_cluster)))

  # Add intra-cluster distances on top
  hist(intra_cluster, breaks = 30, col = rgb(1, 0.5, 0, 0.7), # Orange with transparency 
      add = TRUE)

  # Add legend
  legend("topright", 
        legend = c("Intra-cluster", "Inter-cluster"), 
        fill = c(rgb(1, 0.5, 0, 0.7), rgb(0, 0, 1, 0.7)),
        bty = "n")

  # Save plot if needed
  if(save){
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
    dev.copy(png, filename = paste0(folder, "distance_histogram.png"))
    dev.off()
  }
}

# Function to plot MCMC analysis results - Similarity matrix, trace plots, autocorrelation, posterior distribution
plot_mcmc_results <- function(results, true_labels, BI, save = FALSE, folder = "results/plots/") {
  # Clean previous plots
  graphics.off()

  # Apply burn-in to all data
  k_values <- unlist(results$K)
  
  # Check if K values exist and are valid
  if (is.null(k_values) || length(k_values) == 0) {
    cat("Warning: No valid K values found, skipping cluster distribution plot\n")
    return(NULL)
  }

  # Apply burn-in period
  if (BI > 0 && length(k_values) > BI) {
    k_values <- k_values[(BI + 1):length(k_values)]
    cat("Applied burn-in: Using", length(k_values), "samples after discarding first", BI, "iterations\n")
  } else if (BI > 0) {
    cat("Warning: Burn-in period (", BI, ") is greater than or equal to total iterations (", length(k_values), ")\n")
    cat("Using all available data\n")
  }

  # Remove any NA values
  k_values <- k_values[!is.na(k_values)]

  if (length(k_values) == 0) {
    cat("Warning: All K values are NA after burn-in, skipping cluster distribution plot\n")
    return(NULL)
  }

  ### First plot - Posterior distribution of the number of clusters (after burn-in)
  post_k <- table(k_values) / length(k_values)
  df <- data.frame(
    cluster_found = as.numeric(names(post_k)),
    rel_freq = as.numeric(post_k)
  )

  p1 <- ggplot(data = df, aes(x = factor(cluster_found), y = rel_freq)) +
    geom_col() +
    labs(
      x = "Cluster Found",
      y = "Relative Frequency",
      title = paste("Posterior Distribution of Clusters (After Burn-in:", BI, ")")
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

  if (save)
    ggsave(filename = paste0(folder, "posterior_num_clusters.png"), plot = p1, width = 8, height = 6)

  ### Second plot - Trace of number of clusters
  if (is.null(results$K)) {
    cat("Warning: No valid K values for trace plot\n")
  } 
  else {  
    if (length(k_values) > 0) {
      k_df <- data.frame(
        Iteration = seq_along(k_values),
        NumClusters = k_values
      )

      p2 <- ggplot(k_df, aes(x = Iteration, y = NumClusters)) +
        geom_line() +
        labs(
          x = "Iteration",
          y = "Number of clusters",
          title = paste("Trace Plot (Burn-in:", BI, "iterations)")
        ) +
        theme(
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          text = element_text(size = 15),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "grey95"),
          panel.grid.minor = element_line(color = "grey95"),
          legend.position = "top"
        )
      
      print(p2)
      if (save)
        ggsave(filename = paste0(folder, "traceplot.png"), plot = p2, width = 8, height = 6)
    } else {
      cat("Warning: No valid K values for trace plot after cleaning\n")
    }
  }

  ### Third plot - Posterior Similarity Matrix (using post burn-in data only)
  # Check if allocations exist
  if (is.null(results$allocations) || length(results$allocations) == 0) {
    cat("Warning: No allocation data found, skipping similarity matrix plot\n")
  } else {
    # Apply burn-in to allocations
    allocations_post_burnin <- results$allocations
    if (BI > 0 && length(allocations_post_burnin) > BI) {
      allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
      cat("Using", length(allocations_post_burnin), "allocation samples after burn-in\n")
    }
    
    if (length(allocations_post_burnin) == 0) {
      cat("Warning: No allocations remaining after burn-in\n")
    } else {
      # Compute posterior similarity matrix
      n <- nrow(dist_matrix)
      similarity_matrix <- matrix(0, nrow = n, ncol = n)

      # For each MCMC iteration (post burn-in), add to similarity matrix
      for (iter in seq_along(allocations_post_burnin)) {
        allocation <- allocations_post_burnin[[iter]]
        if (!is.null(allocation) && length(allocation) == n) {
          for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
              if (allocation[i] == allocation[j]) {
                similarity_matrix[i, j] <- similarity_matrix[i, j] + 1
                similarity_matrix[j, i] <- similarity_matrix[j, i] + 1
              }
            }
          }
        }
      }

      # Normalize by number of post burn-in iterations
      similarity_matrix <- similarity_matrix / length(allocations_post_burnin)

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
          title = paste("Posterior Similarity Matrix (After Burn-in:", BI, ")")
        ) +
        theme(
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          text = element_text(size = 12),
          panel.background = element_blank()
        ) +
        coord_fixed()
      print(p3)
      if (save)
        ggsave(filename = paste0(folder, "similarity_matrix.png"), plot = p3, width = 8, height = 6)
    }
  }

  ### Fourth plot - Posterior similarity matrix analysis (using SALSO)
  # Check if allocations exist, otherwise skip this plot
  if (!is.null(results$allocations) && length(results$allocations) > 0) {
    # Apply burn-in to allocations
    allocations_post_burnin <- results$allocations
    if (BI > 0 && length(allocations_post_burnin) > BI) {
      allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
    }
    
    if (length(allocations_post_burnin) > 0) {
      # Convert allocations to matrix format for SALSO
      C <- matrix(unlist(lapply(allocations_post_burnin, function(x) x + 1)),
        nrow = length(allocations_post_burnin),
        ncol = length(true_labels),
        byrow = TRUE
      )

      required_packages <- c("salso", "fields", "viridisLite", "RColorBrewer", "pheatmap")
      for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE)) {
          install.packages(pkg)
          library(pkg, character.only = TRUE)
        }
      }

      # Compute posterior similarity matrix using SALSO
      psm <- salso::psm(C)
      
      # Get point estimate using Variation of Information (VI) loss
      point_estimate <- salso::salso(C, loss = "VI")
      
      # # Reorder based on cluster assignments for better visualization
      # cluster_order <- order(point_estimate)
      # psm_reordered <- psm[cluster_order, cluster_order]
      # true_labels_reordered <- true_labels[cluster_order]
      # point_estimate_reordered <- point_estimate[cluster_order]
      
      # # Plot the reordered similarity matrix
      # pheatmap(psm_reordered,
      #         cluster_rows = FALSE,
      #         cluster_cols = FALSE,
      #         color = colorRampPalette(c("white", "blue"))(100),
      #         main = "Posterior Similarity Matrix (SALSO - Reordered)",
      #         show_rownames = FALSE,
      #         show_colnames = FALSE,
      #         border_color = NA,              # No grid
      # )
      
      # Print results
      cat("=== SALSO Clustering Results (Post Burn-in) ===\n")
      cat("Cluster Sizes:\n")
      print(table(point_estimate))
      
      cat("\nAdjusted Rand Index:", arandi(point_estimate, true_labels), "\n")
      cat("Number of post burn-in samples:", nrow(C), "\n")
      
    } else {
      cat("Warning: No allocation data remaining after burn-in\n")
    }
  } else {
    cat("Warning: No allocation data found, skipping posterior similarity matrix analysis\n")
  }

  ### Fifth plot - Auto-correlation plot (using post burn-in data)
  # Check if loglikelihood data exists and is valid
  if (!is.null(results$loglikelihood) && length(results$loglikelihood) > 0) {
    # Apply burn-in to loglikelihood
    logl_original <- results$loglikelihood
    logl_post_burnin <- logl_original
    k_post_burnin <- k_values  # k_values already has burn-in applied
    
    if (BI > 0 && length(logl_original) > BI) {
      logl_post_burnin <- logl_original[(BI + 1):length(logl_original)]
    }
    
    # Remove any NA or infinite values
    logl_clean <- logl_post_burnin[is.finite(logl_post_burnin)]
    k_clean <- k_post_burnin[is.finite(k_post_burnin)]

    cat("Debug autocorr (post burn-in): logl_clean length =", length(logl_clean), "\n")
    cat("Debug autocorr (post burn-in): k_clean length =", length(k_clean), "\n")

    if (length(logl_clean) > 1 && length(k_clean) > 1) {
      # Ensure both vectors have the same length
      min_length <- min(length(logl_clean), length(k_clean))
      mcmc_list <- list(
        ncls = k_clean[1:min_length],
        logl = logl_clean[1:min_length]
      )
      mcmc_matrix <- do.call(cbind, mcmc_list)

      # Only plot if we have valid data and variance > 0
      if (all(is.finite(mcmc_matrix)) && ncol(mcmc_matrix) > 0 && nrow(mcmc_matrix) > 1) {
        # Check if there's variance in the data
        if (any(apply(mcmc_matrix, 2, var, na.rm = TRUE) > 0)) {
          tryCatch(
            {
              acf(mcmc_matrix, main = paste("Autocorrelation of MCMC chains (Post Burn-in:", BI, ")"))
            },
            error = function(e) {
              cat("Error in MCMC chains autocorrelation:", e$message, "\n")
            }
          )
        } else {
          cat("Warning: No variance in post burn-in MCMC data for autocorrelation\n")
        }
      } else {
        cat("Warning: Skipping autocorrelation plot due to invalid post burn-in data\n")
      }
    } else {
      cat("Warning: Not enough valid post burn-in data for autocorrelation plot\n")
    }

    # Plot loglikelihood autocorrelation separately (post burn-in)
    if (length(logl_clean) > 1 && var(logl_clean, na.rm = TRUE) > 0) {
      tryCatch(
        {
          acf(logl_clean, main = paste("Autocorrelation of Log-Likelihood (Post Burn-in:", BI, ")"))
        },
        error = function(e) {
          cat("Error in loglikelihood autocorrelation:", e$message, "\n")
        }
      )
    } else {
      cat("Warning: Not enough post burn-in loglikelihood data or no variance for autocorrelation\n")
    }
  } else {
    cat("Warning: No loglikelihood data available for autocorrelation plot\n")
  }
}

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

# Function to save results with a specific naming convention in a specified folder
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

  saveRDS(mcmc_result, file = filename_results)
  saveRDS(ground_truth, file = filename_gt)
  saveRDS(all_data, file = filename_data)
  saveRDS(dist_matrix, file = filename_dist)
  saveRDS(hyperparams, file = filename_initial_params)
}

# Function to set hyperparameters using elbow method and k-means clustering
set_hyperparameters <- function(data_coords, dist_matrix, k_elbow, ground_truth = NULL, plot_clustering = FALSE, plot_distribution = TRUE) {

  coords <- data_coords
  
  # Step 2: Use k-means with K_elbow to get initial clustering
  kmeans_result <- kmeans(coords, centers = k_elbow, nstart = 25)
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
        geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", 
                      color = "black", alpha = 0.7) +
        stat_function(fun = dgamma, 
                     args = list(shape = algorithm_shape_A, rate = algorithm_rate_A), 
                     color = "red", linewidth = 1) +
        labs(title = paste("Within-Cluster Distances with Algorithm's Gamma\n",
                          "δ₁ (shape) =", round(algorithm_shape_A, 3), 
                          ", α =", round(alpha, 3), ", β =", round(beta, 3)),
             x = "Distance", y = "Density") +
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
        geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightgreen", 
                      color = "black", alpha = 0.7) +
        stat_function(fun = dgamma, 
                     args = list(shape = algorithm_shape_B, rate = algorithm_rate_B), 
                     color = "blue", linewidth = 1) +
        labs(title = paste("Inter-Cluster Distances with Algorithm's Gamma\n",
                          "δ₂ (shape) =", round(algorithm_shape_B, 3), 
                          ", ζ =", round(zeta, 3), ", γ =", round(gamma_param, 3)),
             x = "Distance", y = "Density") +
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

# Function to plot ground truth and data
plot_data <- function(all_data, cluster_labels, save = FALSE, folder = "results/plots/") {
  ground_truth <- as.factor(cluster_labels)

  # Create data frame for plotting
  plot_data <- data.frame(
    x = all_data[, 1],
    y = all_data[, 2],
    cluster = ground_truth
  )

  # Plot the clusters
  ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 3) +
    labs(
      title = "Clusters",
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme_minimal()

  # Save the plot
  if (save)
    ggsave(filename = paste0(folder, "data_clusters.png"), width = 8, height = 6)
}