library(mcclust.ext) # For MCMC clustering functions
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation
library(tidyr)      # For data reshaping

# Function to plot MCMC analysis results
plot_mcmc_results <- function(results, true_labels) {
  # Clean previous plots
  graphics.off()

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

  acf(results$loglikelihood, main = "Autocorrelation of Log-Likelihood")
}

# Load the results and ground truth for later analysis
folder <- "results/"
filename_results <- "simulation_results.rds"
filename_gt <- "simulation_ground_truth.rds"
filename_dist <- "simulation_distance_matrix.rds"

filename_results <- paste0(folder, filename_results)
filename_gt <- paste0(folder, filename_gt)
filename_dist <- paste0(folder, filename_dist)

results <- readRDS(file = filename_results)
ground_truth <- readRDS(file = filename_gt)
dist_matrix <- readRDS(file = filename_dist)

plot_mcmc_results(results, ground_truth)
