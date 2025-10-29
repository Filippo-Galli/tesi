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

plot_distance <- function(dist_matrix, cls = NULL, save = FALSE, folder = "results/plots/") {
  # Calculate distance matrix
  dist_matrix <- as.matrix(dist_matrix)

  cls_exist <- !is.null(cls)

  if(!cls_exist) {
    cls <- rep(1, nrow(dist_matrix)) # All in one cluster
  }

  # Get upper triangular indices to avoid duplicates
  upper_tri_indices <- which(upper.tri(dist_matrix), arr.ind = TRUE)

  # Extract distances and determine if they are intra or inter-cluster
  distances <- dist_matrix[upper_tri_indices]
  cluster_pairs <- cbind(cls[upper_tri_indices[, 1]], cls[upper_tri_indices[, 2]])

  # Classify distances as intra-cluster (same cluster) or inter-cluster (different clusters)
  intra_cluster <- distances[cluster_pairs[, 1] == cluster_pairs[, 2]]
  inter_cluster <- distances[cluster_pairs[, 1] != cluster_pairs[, 2]]

  # Create histogram with overlaid distributions
  hist(intra_cluster,
    breaks = 30, col = rgb(1, 0.5, 0, 0.7), # Orange with transparency
    main = "Histogram of Pairwise Distances",
    xlab = "Distance",
    ylab = "Frequency",
    xlim = range(c(intra_cluster, inter_cluster)),
    ylim = c(0, max(
      hist(intra_cluster, plot = FALSE)$counts,
      if(cls_exist) hist(inter_cluster, plot = FALSE)$counts else 0
    ) + 5)
  )

  if(cls_exist) {
    hist(inter_cluster,
      breaks = 30, col = rgb(0, 0, 1, 0.7), # Orange with transparency
      add = TRUE
    )
  }

  # Add legend
  legend("topright",
    legend = ifelse(cls_exist, c("Intra-cluster", "Inter-cluster"), "Intra-cluster"),
    fill = c(rgb(1, 0.5, 0, 0.7), rgb(0, 0, 1, 0.7)),
    bty = "n"
  )

  # Save plot if needed
  if (save) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
    dev.copy(png, filename = paste0(folder, "distance_histogram.png"), width = 2400, height = 1800, res = 300)
    dev.off()
  }
}

plot_post_distr <- function(results, BI, save = FALSE, folder = "results/plots/") {
  # Apply burn-in to all data
  k_values <- unlist(results$K)

  # Apply burn-in period
  if (BI > 0 && length(k_values) > BI) {
    k_values <- k_values[(BI + 1):length(k_values)]
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

  if (save) {
    ggsave(filename = paste0(folder, "posterior_num_clusters.png"),
           plot = p1, width = 8, height = 6)
  }
}

plot_trace_cls <- function(results, BI, save = FALSE, folder = "results/plots/") {
  k_values <- unlist(results$K)

  # Apply burn-in period
  if (BI > 0 && length(k_values) > BI) {
    k_values <- k_values[(BI + 1):length(k_values)]
  }

  ### Second plot - Trace of number of clusters
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
  if (save) {
    ggsave(filename = paste0(folder, "traceplot.png"),
           plot = p2, width = 8, height = 6)
  }
}

plot_post_sim_matrix <- function(results, BI, save = FALSE, folder = "results/plots/"){
  #### Apply burn-in to allocations
  allocations_post_burnin <- results$allocations
  if (BI > 0 && length(allocations_post_burnin) > BI) {
    allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
  }

  #### Compute posterior similarity matrix
  n <- length(allocations_post_burnin[[1]])
  similarity_matrix <- matrix(0, nrow = n, ncol = n)

  # For each MCMC iteration (post burn-in), add to similarity matrix
  for (iter in seq_along(allocations_post_burnin)) {
    allocation <- allocations_post_burnin[[iter]]
    if (!is.null(allocation) && length(allocation) == n) {
      # This creates a matrix where entry [i,j] is TRUE if allocation[i] == allocation[j]
      same_cluster <- outer(allocation, allocation, "==")
      
      # Add to similarity matrix (this is already symmetric, so no need for separate i,j loop)
      similarity_matrix <- similarity_matrix + same_cluster
    }
  }

  # Normalize by number of post burn-in iterations
  similarity_matrix <- similarity_matrix / length(allocations_post_burnin)

  # Set diagonal to 1 (each point is always similar to itself)
  diag(similarity_matrix) <- 1

  # 1. Hierarchical clustering to reorder
  hc <- hclust(as.dist(1 - similarity_matrix))
  ordered_idx <- hc$order
  sim_reordered <- similarity_matrix[ordered_idx, ordered_idx]

  # 2. Reshape for ggplot
  sim_df <- melt(sim_reordered)
  colnames(sim_df) <- c("i", "j", "similarity")

  # 3. Build heatmap with ggplot2
  p3 <- ggplot(sim_df, aes(x = j, y = i, fill = similarity)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkblue") +
    labs(
      x = "Data Point Index (clustered)",
      y = "Data Point Index (clustered)",
      fill = "Posterior\nSimilarity",
      title = "Posterior Similarity Matrix (clustered order)"
    ) +
    theme_minimal(base_size = 13) +
    coord_fixed()
  print(p3)
  if (save) {
    ggsave(filename = paste0(folder, "similarity_matrix.png"), 
           plot = p3, width = 8, height = 6)
  }
}

plot_stats <- function(results, true_labels, BI, save = FALSE, folder = "results/plots/") {
  #### Apply burn-in to allocations
  allocations_post_burnin <- results$allocations
  if (BI > 0 && length(allocations_post_burnin) > BI) {
    allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
  }
  
  #### Convert allocations to matrix format for SALSO
  C <- matrix(unlist(lapply(allocations_post_burnin, function(x) x + 1)),
    nrow = length(allocations_post_burnin),
    ncol = length(true_labels),
    byrow = TRUE
  )

  #### Get point estimate using Variation of Information (VI) loss
  point_estimate <- salso::salso(C, loss = "binder",
                                 maxNClusters = 200,
                                 maxZealousAttempts = 1000)

  #### Print results
  cat("=== SALSO Clustering Results (Post Burn-in) ===\n")
  cat("Cluster Sizes:\n")
  print(table(point_estimate))

  point_labels <- as.vector(point_estimate)

  #### Adjusted Rand Index (ARI)
  cat("\nAdjusted Rand Index:",
      arandi(point_estimate, true_labels), "\n")

  #### NMI
  cat("Normalized Mutual Information:",
      NMI(point_labels, true_labels), "\n")

  #### VI
  cat("Variation of Information:", NVI(point_labels, true_labels), "\n")

  if (save) {
    stats_file <- paste0(folder, "salso_stats.txt")
    write(paste("Adjusted Rand Index:", arandi(point_estimate, true_labels)), file = stats_file)
    write(paste("Normalized Mutual Information:", NMI(point_labels, true_labels)), file = stats_file, append = TRUE)
    write(paste("Variation of Information:", NVI(point_labels, true_labels)), file = stats_file, append = TRUE)
  }
}

plot_trace_U <- function(results, BI, save = FALSE, folder = "results/plots/") {
  U_after_burnin <- results$U
  if (BI > 0 && length(U_after_burnin) > BI) {
    U_after_burnin <- U_after_burnin[(BI + 1):length(U_after_burnin)]
  }
  plot(U_after_burnin, type = "l", xlab = "Iteration", ylab = "U", width = 2400, height = 1800, res = 300)
  abline(h = mean(U_after_burnin), col = "red", lty = 2)
  legend("topright", legend = c("Mean U"), col = c("red"), lty = 2)
  titolo <- paste0("Trace of U over MCMC iterations (mean U = ", round(mean(U_after_burnin), 3), ")")
  title(main = titolo)

  if (save) {
    dev.copy(png, filename = paste0(folder, "U_trace"), width = 2400, height = 1800, res = 300)
  }
}

plot_acf_U <- function(results, BI, save = FALSE, folder = "results/plots/") {
  
  U_after_burnin <- results$U
  if (BI > 0 && length(U_after_burnin) > BI) {
    U_after_burnin <- U_after_burnin[(BI + 1):length(U_after_burnin)]
  }

  acf(U_after_burnin, main = "ACF of U over MCMC iterations")
  if (save) {
    dev.copy(png, filename = paste0(folder, "U_acf"), width = 2400, height = 1800, res = 300)
  }
}

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

plot_data <- function(all_data, cluster_labels, save = FALSE, folder = "results/plots/") {
  ground_truth <- as.factor(cluster_labels)

  # Create data frame for plotting
  plot_data <- data.frame(
    x = all_data[, 1],
    y = all_data[, 2],
    cluster = ground_truth
  )

  # Plot the clusters
  p <- ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 3) +
    labs(
      title = "Clusters",
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme_minimal()

  print(p)

  # Save the plot
  if (save) {
    ggsave(filename = paste0(folder, "data_clusters.png"), 
           plot = p3, width = 8, height = 6)
  }
}