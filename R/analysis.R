library(mcclust.ext) # For MCMC clustering functions
library(ggplot2) # For plotting
library(dplyr) # For data manipulation
library(tidyr) # For data reshaping

# Function to plot MCMC analysis results
plot_mcmc_results <- function(results, true_labels) {
  # Clean previous plots
  graphics.off()

  # Check if results object has required components
  if (is.null(results) || is.null(results$K)) {
    cat("Error: Results object is NULL or missing K component\n")
    return(NULL)
  }

  # Debug: Check structure of results object
  cat("Debugging results structure:\n")
  cat("Class of results:", class(results), "\n")
  cat("Names of results:", names(results), "\n")
  if (!is.null(results$K)) {
    cat("Length of K:", length(results$K), "\n")
  } else {
    cat("K is NULL\n")
  }
  if (!is.null(results$allocations)) {
    cat("Length of allocations:", length(results$allocations), "\n")
  } else {
    cat("allocations is NULL\n")
  }
  if (!is.null(results$loglikelihood)) {
    cat("Length of loglikelihood:", length(results$loglikelihood), "\n")
  } else {
    cat("loglikelihood is NULL\n")
  }

  ### First plot - Posterior distribution of the number of clusters
  k_values <- unlist(results$K)

  # Check if K values exist and are valid
  if (is.null(k_values) || length(k_values) == 0) {
    cat("Warning: No valid K values found, skipping cluster distribution plot\n")
    return(NULL)
  }

  # Remove any NA values
  k_values <- k_values[!is.na(k_values)]

  if (length(k_values) == 0) {
    cat("Warning: All K values are NA, skipping cluster distribution plot\n")
    return(NULL)
  }

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

  ggsave(filename = paste0(folder, "posterior_num_clusters.png"), plot = p1, width = 8, height = 6)

  ### Second plot - Trace of number of clusters
  # Ensure K values exist and handle any missing values
  if (is.null(results$K) || length(k_values) == 0) {
    cat("Warning: No valid K values for trace plot\n")
  } else {
    # Remove any NULL or NA values and create corresponding iteration indices
    valid_k <- k_values[!is.na(k_values)]

    if (length(valid_k) > 0) {
      k_df <- data.frame(
        Iteration = seq_along(valid_k),
        NumClusters = valid_k
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

      ggsave(filename = paste0(folder, "traceplot.png"), plot = p2, width = 8, height = 6)
    } else {
      cat("Warning: No valid K values for trace plot after cleaning\n")
    }
  }

  ### Third plot - Posterior Similarity Matrix
  # Check if allocations exist
  if (is.null(results$allocations) || length(results$allocations) == 0) {
    cat("Warning: No allocation data found, skipping similarity matrix plot\n")
  } else {
    # Compute posterior similarity matrix
    n <- nrow(dist_matrix)
    similarity_matrix <- matrix(0, nrow = n, ncol = n)

    # For each MCMC iteration, add to similarity matrix
    for (iter in seq_along(results$allocations)) {
      allocation <- results$allocations[[iter]]
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
    ggsave(filename = paste0(folder, "similarity_matrix.png"), plot = p3, width = 8, height = 6)
  }

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
    k_clean <- k_values[is.finite(k_values)]

    cat("Debug autocorr: logl_clean length =", length(logl_clean), "\n")
    cat("Debug autocorr: k_clean length =", length(k_clean), "\n")

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
              acf(mcmc_matrix, main = "Autocorrelation of MCMC chains")
            },
            error = function(e) {
              cat("Error in MCMC chains autocorrelation:", e$message, "\n")
            }
          )
        } else {
          cat("Warning: No variance in MCMC data for autocorrelation\n")
        }
      } else {
        cat("Warning: Skipping autocorrelation plot due to invalid data\n")
      }
    } else {
      cat("Warning: Not enough valid data for autocorrelation plot\n")
    }

    # Plot loglikelihood autocorrelation separately
    if (length(logl_clean) > 1 && var(logl_clean, na.rm = TRUE) > 0) {
      tryCatch(
        {
          acf(logl_clean, main = "Autocorrelation of Log-Likelihood")
        },
        error = function(e) {
          cat("Error in loglikelihood autocorrelation:", e$message, "\n")
        }
      )
    } else {
      cat("Warning: Not enough loglikelihood data or no variance for autocorrelation\n")
    }
  } else {
    cat("Warning: No loglikelihood data available for autocorrelation plot\n")
  }
}

list.files("results/")

# Load the results and ground truth for later analysis
folder <- "results/OneInOne_2D_BI1000_NI5000_a4_sigma1_tau1/"
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

ground_truth <- as.factor(ground_truth)

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
ggsave(filename = paste0(folder, "data_clusters.png"), width = 8, height = 6)

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

plot_mcmc_results(results, ground_truth)
