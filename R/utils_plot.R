## @file utils.R
## @brief Utility functions for Bayesian clustering analysis and MCMC visualization
## @author Filippo Galli
## @date 2025

suppressMessages(library(Rcpp))
suppressMessages(library(ggplot2))
suppressMessages(library(MASS)) # For fitdistr function
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(spam)) # For comp.psm
suppressMessages(library(fields)) # For minVI
suppressMessages(library(viridisLite)) # For color scales
suppressMessages(library(RColorBrewer)) # For color palettes
suppressMessages(library(pheatmap)) # For heatmaps
suppressMessages(library(mcclust.ext)) # For MCMC clustering functions
suppressMessages(library(mvtnorm))
suppressMessages(library(gtools))
suppressMessages(library(salso))
suppressMessages(library(aricode))
suppressMessages(library(reshape2))
suppressMessages(library(cluster))
suppressMessages(library(label.switching))
suppressMessages(library(coda))

## @brief Generate consistent colors for clusters
## @param n_clusters Number of clusters
## @return Named vector of colors
get_cluster_colors <- function(n_clusters) {
  # Use viridis turbo for many clusters
  colors <- viridisLite::turbo(n_clusters)

  # Return named vector with cluster numbers as names
  names(colors) <- as.character(seq_len(n_clusters))
  return(colors)
}

plot_distance <- function(dist_matrix, cls = NULL,
                          save = FALSE, folder = "results/plots/",
                          title = "Distance Histogram",
                          normalize = FALSE,
                          breaks = 30) {
  # Calculate distance matrix
  dist_matrix <- as.matrix(dist_matrix)

  cls_exist <- !is.null(cls)

  if (!cls_exist) {
    cls <- rep(1, nrow(dist_matrix)) # All in one cluster
  }

  # Get upper triangular indices to avoid duplicates
  upper_tri_indices <- which(upper.tri(dist_matrix), arr.ind = TRUE)

  # Extract distances and determine if they are intra or inter-cluster
  distances <- dist_matrix[upper_tri_indices]
  cluster_pairs <- cbind(
    cls[upper_tri_indices[, 1]],
    cls[upper_tri_indices[, 2]]
  )

  # Classify distances as intra-cluster or inter-cluster
  intra_cluster <- distances[cluster_pairs[, 1] == cluster_pairs[, 2]]
  inter_cluster <- distances[cluster_pairs[, 1] != cluster_pairs[, 2]]

  # Pre-compute histograms for ylim calculation
  h_intra <- hist(intra_cluster, breaks = breaks, plot = FALSE)
  h_inter <- if (cls_exist) hist(inter_cluster, breaks = breaks, plot = FALSE) else NULL

  # Determine y-axis label and limits based on normalization
  ylab <- if (normalize) "Density" else "Frequency"

  ylim_max <- if (normalize) {
    max(h_intra$density, if (cls_exist) h_inter$density else 0) * 1.1
  } else {
    max(h_intra$counts, if (cls_exist) h_inter$counts else 0) * 1.1
  }

  # Create histogram with overlaid distributions
  hist(intra_cluster,
    breaks = breaks,
    col = rgb(1, 0.5, 0, 0.7),
    main = "",
    xlab = "Distance",
    ylab = ylab,
    probability = normalize,
    xlim = range(c(intra_cluster, inter_cluster)),
    ylim = c(0, ylim_max)
  )

  if (cls_exist) {
    hist(inter_cluster,
      breaks = breaks,
      col = rgb(0, 0, 1, 0.7),
      probability = normalize,
      add = TRUE
    )
  }

  # Add legend
  legend("topright",
    legend = if (cls_exist) c("Intra-cluster", "Inter-cluster") else "Distance",
    fill = if (cls_exist) c(rgb(1, 0.5, 0, 0.7), rgb(0, 0, 1, 0.7)) else rgb(1, 0.5, 0, 0.7),
    bty = "n"
  )
  title(main = title)

  # Save plot if needed
  if (save) {
    if (!dir.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
    words_title <- gsub(" ", "_", title)
    dev.copy(png,
      filename = paste0(folder, words_title, ".png"),
      width = 2400, height = 1800, res = 300
    )
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
    ggsave(
      filename = paste0(folder, "posterior_num_clusters.png"),
      plot = p1, width = 8, height = 6
    )
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
    ggsave(
      filename = paste0(folder, "traceplot.png"),
      plot = p2, width = 8, height = 6
    )
  }
}

plot_post_sim_matrix <- function(results, BI, save = FALSE, folder = "results/plots/") {
  #### Apply burn-in to allocations
  allocations_post_burnin <- results$allocations
  if (BI > 0 && length(allocations_post_burnin) > BI) {
    allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
  }

  #### Compute posterior similarity matrix using salso::psm
  n <- length(allocations_post_burnin[[1]])
  n_iter <- length(allocations_post_burnin)

  cat("Data dimensions: n =", n, ", n_iter =", n_iter, "\n")

  # Convert allocations to matrix format (each row is one iteration)
  alloc_matrix <- matrix(unlist(allocations_post_burnin),
    nrow = n_iter,
    ncol = n,
    byrow = TRUE
  )

  cat("Computing posterior similarity matrix using salso::psm...\n")
  similarity_matrix <- salso::psm(alloc_matrix)

  # Hierarchical clustering
  cat("Performing hierarchical clustering...\n")
  hc <- hclust(as.dist(1 - similarity_matrix))
  ordered_idx <- hc$order
  sim_reordered <- similarity_matrix[ordered_idx, ordered_idx]

  cat("Creating plot (full", n, "x", n, "matrix)...\n")

  # Convert reordered similarity matrix to long format for ggplot2
  # Use indices directly to avoid large melt operation
  sim_df <- data.frame(
    i = rep(1:n, each = n),
    j = rep(1:n, times = n),
    similarity = as.vector(sim_reordered)
  )

  # Create ggplot with geom_raster (more efficient than geom_tile)
  p <- ggplot(sim_df, aes(x = j, y = i, fill = similarity)) +
    geom_raster() +
    scale_fill_gradient(
      low = "white", high = "darkblue",
      limits = c(0, 1),
      name = "Posterior\nSimilarity"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Data Point Index (clustered)",
      y = "Data Point Index (clustered)",
      title = "Posterior Similarity Matrix"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    coord_fixed()

  print(p)

  if (save) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    ggsave(
      filename = paste0(folder, "similarity_matrix.png"),
      plot = p,
      width = 10, height = 10, dpi = 400
    )
    cat("Plot saved to:", paste0(folder, "similarity_matrix.png"), "\n")
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
  point_estimate <- salso::salso(C,
    loss = "VI",
    maxNClusters = 200,
    maxZealousAttempts = 1000
  )

  #### Print results
  cat("=== SALSO Clustering Results (Post Burn-in) ===\n")
  cat("Cluster Sizes:\n")
  print(table(point_estimate))

  point_labels <- as.vector(point_estimate)

  #### Adjusted Rand Index (ARI)
  cat(
    "\nAdjusted Rand Index:",
    arandi(point_estimate, true_labels), "\n"
  )

  #### NMI
  cat(
    "Normalized Mutual Information:",
    NMI(point_labels, true_labels), "\n"
  )

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

  draw_plot <- function() {
    plot(U_after_burnin, type = "l", xlab = "Iteration", ylab = "U")
    abline(h = mean(U_after_burnin), col = "red", lty = 2)
    legend("topright",
      legend = sprintf("Mean U = %.3f", mean(U_after_burnin)),
      col = "red", lty = 2, bty = "n"
    )
    title(main = sprintf("Trace of U over MCMC iterations\n(mean U = %.3f)", mean(U_after_burnin)))
  }

  draw_plot()

  if (save) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    file <- file.path(folder, "U_trace.png")
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file, width = 2400, height = 1800, res = 300)
    } else {
      png(file, width = 2400, height = 1800, res = 300)
    }
    draw_plot()
    dev.off()
  }
}

plot_acf_U <- function(results, BI, save = FALSE, folder = "results/plots/") {
  U_after_burnin <- results$U
  if (BI > 0 && length(U_after_burnin) > BI) {
    U_after_burnin <- U_after_burnin[(BI + 1):length(U_after_burnin)]
  }

  draw_plot <- function() {
    acf(U_after_burnin, main = "ACF of U over MCMC iterations")
  }

  draw_plot()

  if (save) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    file <- file.path(folder, "U_acf.png")
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file, width = 2400, height = 1800, res = 300)
    } else {
      png(file, width = 2400, height = 1800, res = 300)
    }
    draw_plot()
    dev.off()
  }
}

plot_k_means <- function(dist_matrix, max_k = 10) {
  # Ensure dist_matrix is in the right format
  if (!inherits(dist_matrix, "dist")) {
    dist_matrix <- as.dist(dist_matrix)
  }

  mds_result <- cmdscale(dist_matrix, k = 2)

  # Step 1: Compute elbow method for K selection using kmeans
  wss <- numeric(max_k)

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
      title = "Elbow Method for Optimal K (k-means)",
      x = "Number of Clusters (K)",
      y = "Total Dissimilarity"
    ) +
    theme_minimal()

  print(elbow_plot)

  # Return the WSS values for further analysis
  return(invisible(elbow_data))
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
    ggsave(
      filename = paste0(folder, "data_clusters.png"),
      plot = p3, width = 8, height = 6
    )
  }
}

plot_cls_est <- function(results, BI, save = FALSE, start_time, end_time, folder = "results/plots/") {
  cat("Computing point estimate using SALSO...\n")
  #### Apply burn-in to allocations
  allocations_post_burnin <- results$allocations
  if (BI > 0 && length(allocations_post_burnin) > BI) {
    allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
  }

  #### Convert allocations to matrix format for SALSO
  C <- matrix(unlist(lapply(allocations_post_burnin, function(x) x + 1)),
    nrow = length(allocations_post_burnin),
    ncol = length(allocations_post_burnin[[1]]),
    byrow = TRUE
  )

  #### Get point estimate using Variation of Information (VI) loss
  point_estimate <- salso::salso(C,
    loss = "VI",
    maxNClusters = 200,
    maxZealousAttempts = 1000
  )

  #### Print results
  cat("=== SALSO Clustering Results (Post Burn-in) ===\n")
  cat("Cluster Sizes:\n")
  print(table(point_estimate))

  cat("=== ESS ===\n")
  ess_value <- effectiveSize(as.mcmc(unlist(results$K)))
  ess_seconds <- ess_value / elapsed_time
  cat("Effective Sample Size (ESS):", ess_value, "\n")
  cat("ESS per second:", ess_seconds, "\n")

  if (save) {
    cls_file <- paste0(folder, "salso_cluster_estimate.txt")
    write("Cluster Sizes:", file = cls_file)
    write(capture.output(table(point_estimate)), file = cls_file, append = TRUE)
    write("\n=== ESS ===", file = cls_file, append = TRUE)
    write(paste("ESS:", ess_value), file = cls_file, append = TRUE)
    write(paste("ESS per second:", ess_seconds), file = cls_file, append = TRUE)
  }
  cat("Computed point estimate using SALSO...\n")
  return(invisible(point_estimate))
}

plot_map_cls <- function(results, BI, point_estimate = NULL, save = FALSE,
                         folder = "results/plots/",
                         puma_dir = "input/counties-pumas",
                         id_col = "PUMA",
                         unit_ids = NULL,
                         simplify_geom = TRUE,
                         dTolerance = 100) {
  cat("Computing map of cluster assignments...\n")

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for plot_map_cls().")
  }

  # Only compute point_estimate if not provided
  if (is.null(point_estimate)) {
    point_estimate <- plot_cls_est(results, BI, save = FALSE)
  }

  shp <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp) == 0) {
    stop("No .shp file found in '", puma_dir, "'.")
  }

  # Assign names to point_estimate if missing
  if (is.null(names(point_estimate))) {
    candidate_ids <- unit_ids %||% results$unit_ids %||% results$puma_ids
    if (is.null(candidate_ids) || length(candidate_ids) != length(point_estimate)) {
      stop("Provide unit_ids or store results$unit_ids/results$puma_ids matching the PUMAs.")
    }
    names(point_estimate) <- candidate_ids
  }

  # Pre-compute cluster factor with proper levels (avoids recomputation)
  unique_clusters <- sort(unique(point_estimate))
  n_clusters <- length(unique_clusters)
  cluster_colors <- get_cluster_colors(n_clusters)

  # Create cluster data frame with pre-factored clusters
  cluster_df <- data.frame(
    id = names(point_estimate),
    cluster = factor(point_estimate, levels = unique_clusters),
    stringsAsFactors = FALSE
  )
  names(cluster_df)[1] <- id_col

  # Read shapefile - only columns we need
  geom <- sf::st_read(shp[1], quiet = TRUE)
  geom <- geom[, c(id_col, attr(geom, "sf_column"))]
  geom[[id_col]] <- as.character(geom[[id_col]])

  # Use match-based join (faster than merge for simple cases)
  idx <- match(geom[[id_col]], cluster_df[[id_col]])
  keep <- !is.na(idx)
  geom <- geom[keep, ]
  geom$cluster <- cluster_df$cluster[idx[keep]]

  # Simplify geometry for faster rendering (optional)
  if (simplify_geom && nrow(geom) > 50) {
    geom <- sf::st_simplify(geom, preserveTopology = TRUE, dTolerance = dTolerance)
  }

  p <- ggplot2::ggplot(geom) +
    ggplot2::geom_sf(ggplot2::aes(fill = cluster), color = "grey60", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = cluster_colors, drop = FALSE, na.value = "lightgrey") +
    ggplot2::labs(
      title = "PUMAs by Cluster Assignment",
      fill = "Cluster"
    ) +
    ggplot2::theme_minimal()

  print(p)

  if (save) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    ggplot2::ggsave(file.path(folder, "puma_clusters.png"), p, width = 10, height = 8)
  }

  invisible(p)
}

plot_map_prior_mean <- function(save = FALSE, folder = "results/plots/",
                                puma_dir = "input/counties-pumas",
                                input_dir = "input/CA/",
                                id_col = "PUMA",
                                unit_ids = NULL) {
  cat("Computing map of mean income ...\n")

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for plot_map_prior_mean().")
  }

  shp <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp) == 0) {
    stop("No .shp file found in '", puma_dir, "'.")
  }

  # Try multiple file formats
  data <- NULL

  # Try .rds first
  rds_file <- paste0(input_dir, "full_dataset.rds")
  if (file.exists(rds_file)) {
    data <- readRDS(rds_file)
  }

  # Try .csv if .rds not found
  if (is.null(data)) {
    csv_file <- paste0(input_dir, "full_dataset.csv")
    if (file.exists(csv_file)) {
      data <- read.csv(csv_file)
    }
  }

  # Try .dat with readRDS as fallback
  if (is.null(data)) {
    dat_file <- paste0(input_dir, "full_dataset.dat")
    if (file.exists(dat_file)) {
      data <- tryCatch(
        readRDS(dat_file),
        error = function(e) {
          env <- new.env()
          load(dat_file, envir = env)
          get(ls(env)[1], envir = env)
        }
      )
    }
  }

  if (is.null(data)) {
    stop("Could not find or load full_dataset in any supported format (.rds, .csv, .dat) in ", input_dir)
  }

  if (is.null(unit_ids)) {
    unit_ids <- names(data)
  }
  if (is.null(unit_ids) || length(unit_ids) != length(data)) {
    stop("Provide unit_ids matching the number of PUMAs in 'data'.")
  }

  prior_means <- vapply(data, mean, numeric(1))
  names(prior_means) <- unit_ids

  geom <- sf::st_read(shp[1], quiet = TRUE)
  prior_df <- tibble::tibble(
    !!id_col := as.character(unit_ids),
    prior_mean = prior_means
  )
  geom[[id_col]] <- as.character(geom[[id_col]])
  geom <- dplyr::left_join(geom, prior_df, by = id_col)

  p <- ggplot2::ggplot(geom) +
    ggplot2::geom_sf(aes(fill = prior_mean), color = "grey60", size = 0.2) +
    ggplot2::scale_fill_viridis_c(option = "viridis", na.value = "lightgrey") +
    ggplot2::labs(
      title = "PUMAs by Prior Mean",
      fill = "Prior Mean"
    ) +
    ggplot2::theme_minimal()

  print(p)

  if (save) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    ggplot2::ggsave(file.path(folder, "puma_prior_means.png"), p, width = 10, height = 8)
  }

  invisible(p)
}

plot_hist_cls_pumas <- function(results, BI, input_dir = "input/CA/", point_estimate = NULL, save = FALSE, folder = "results/plots/", log = TRUE) {
  cat("Computing histograms of cluster assignments...\n")

  # Try multiple file formats
  data <- NULL

  # Try .rds first
  rds_file <- paste0(input_dir, "full_dataset.rds")
  if (file.exists(rds_file)) {
    data <- readRDS(rds_file)
  }

  # Try .csv if .rds not found
  if (is.null(data)) {
    csv_file <- paste0(input_dir, "full_dataset.csv")
    if (file.exists(csv_file)) {
      data <- read.csv(csv_file)
    }
  }

  # Try .dat with readRDS as fallback
  if (is.null(data)) {
    dat_file <- paste0(input_dir, "full_dataset.dat")
    if (file.exists(dat_file)) {
      data <- tryCatch(
        readRDS(dat_file),
        error = function(e) {
          env <- new.env()
          load(dat_file, envir = env)
          get(ls(env)[1], envir = env)
        }
      )
    }
  }

  if (is.null(data)) {
    stop("Could not find or load full_dataset in any supported format (.rds, .csv, .dat) in ", input_dir)
  }

  # Transform data from log scale if needed
  if (log) {
    cat("Data is in log scale, transforming to original scale (exp)...\n")
    data <- lapply(data, exp)
  }

  if (is.null(point_estimate)) {
    point_estimate <- plot_cls_est(results, BI, save = FALSE)
  }

  unique_clusters <- sort(unique(point_estimate))
  n_clusters <- length(unique_clusters)

  cat("\n=== Cluster Diagnostics ===\n")
  cat("Number of clusters found:", n_clusters, "\n")
  cat("Cluster sizes:\n")
  print(table(point_estimate))

  draw_histograms <- function() {
    if (n_clusters == 1) {
      cat("\n⚠️  WARNING: Only 1 cluster found. Showing overall distribution.\n")
      combined_data <- unlist(data)
      hist(combined_data,
        breaks = "FD",
        main = paste("Single Cluster Distribution\n(n_pumas =", length(data), ")"),
        xlab = "Income Value",
        ylab = "Density",
        col = "steelblue",
        border = "white",
        probability = TRUE
      )
      if (length(unique(combined_data)) > 1) {
        lines(density(combined_data), col = "red", lwd = 2)
      }
      legend("topleft",
        legend = c(
          paste("Mean:", round(mean(combined_data), 2)),
          paste("SD:", round(sd(combined_data), 2)),
          paste("N obs:", length(combined_data))
        ),
        bty = "n",
        cex = 0.9
      )

      cat("\nOverall statistics:\n")
      cat("Mean:", mean(combined_data), "\n")
      cat("SD:", sd(combined_data), "\n")
      cat("Min:", min(combined_data), "\n")
      cat("Max:", max(combined_data), "\n")

      puma_means <- sapply(data, mean)
      cat("\nPUMA-level variation:\n")
      cat("Mean of means:", mean(puma_means), "\n")
      cat("SD of means:", sd(puma_means), "\n")
      cat("Range of means:", range(puma_means), "\n")
    } else {
      clusters <- point_estimate
      data_split <- split(data, as.factor(clusters))
      op <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(op), add = TRUE)

      # Get consistent colors for clusters
      cluster_colors <- get_cluster_colors(n_clusters)

      # Split into batches of max 9 clusters per plot
      max_plots_per_page <- 9
      cluster_names <- names(data_split)
      n_batches <- ceiling(n_clusters / max_plots_per_page)

      for (batch_idx in seq_len(n_batches)) {
        # Determine which clusters go in this batch
        start_idx <- (batch_idx - 1) * max_plots_per_page + 1
        end_idx <- min(batch_idx * max_plots_per_page, n_clusters)
        batch_clusters <- cluster_names[start_idx:end_idx]
        n_plots_in_batch <- length(batch_clusters)

        # Calculate layout for this batch
        n_rows <- ceiling(sqrt(n_plots_in_batch))
        n_cols <- ceiling(n_plots_in_batch / n_rows)
        graphics::par(mfrow = c(n_rows, n_cols))

        for (cl in batch_clusters) {
          cluster_data <- data_split[[cl]]
          n_pumas <- length(cluster_data)
          combined_data <- unlist(cluster_data)
          # Use consistent cluster colors with transparency
          cl_color <- adjustcolor(cluster_colors[cl], alpha.f = 0.7)
          hist(combined_data,
            breaks = 30,
            main = paste("Cluster", cl, "\n(n =", n_pumas, "PUMAs)"),
            xlab = "Income Value",
            ylab = "Density",
            col = cl_color,
            border = "white",
            probability = TRUE
          )
          if (length(unique(combined_data)) > 1) {
            lines(density(combined_data), col = "black", lwd = 2)
          }
          cat("\nCluster", cl, "statistics:\n")
          cat("  N PUMAs:", n_pumas, "\n")
          cat("  N observations:", length(combined_data), "\n")
          cat("  Mean:", mean(combined_data), "\n")
          cat("  SD:", sd(combined_data), "\n")

          legend("topleft",
            legend = c(
              paste("Mean:", round(mean(combined_data), 2)),
              paste("SD:", round(sd(combined_data), 2)),
              paste("N obs:", length(combined_data))
            ),
            bty = "n",
            cex = 0.9
          )
        }

        # If saving and multiple batches, save each batch separately
        if (save && n_batches > 1) {
          if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
          out_file <- file.path(folder, paste0("cluster_histograms_batch_", batch_idx, ".png"))
          if (requireNamespace("ragg", quietly = TRUE)) {
            ragg::agg_png(filename = out_file, width = 2400, height = 1800, res = 300)
          } else {
            png_type <- if (capabilities("cairo")) "cairo" else if (capabilities("X11")) "Xlib" else "cairo"
            grDevices::png(filename = out_file, width = 2400, height = 1800, res = 300, type = png_type)
          }

          # Redraw the plots for this batch
          graphics::par(mfrow = c(n_rows, n_cols))
          for (cl in batch_clusters) {
            cluster_data <- data_split[[cl]]
            n_pumas <- length(cluster_data)
            combined_data <- unlist(cluster_data)
            cl_color <- adjustcolor(cluster_colors[cl], alpha.f = 0.7)
            hist(combined_data,
              breaks = 30,
              main = paste("Cluster", cl, "\n(n =", n_pumas, "PUMAs)"),
              xlab = "Income Value",
              ylab = "Density",
              col = cl_color,
              border = "white",
              probability = TRUE
            )
            if (length(unique(combined_data)) > 1) {
              lines(density(combined_data), col = "black", lwd = 2)
            }
            legend("topleft",
              legend = c(
                paste("Mean:", round(mean(combined_data), 2)),
                paste("SD:", round(sd(combined_data), 2)),
                paste("N obs:", length(combined_data))
              ),
              bty = "n",
              cex = 0.9
            )
          }

          grDevices::dev.off()
          cat("\nSaved batch", batch_idx, "of", n_batches, "to", out_file, "\n")
        }
      }
    }
  }

  draw_histograms()

  # Only save here if we have <= 9 clusters (single file case)
  # Multiple batch saving is handled within draw_histograms()
  device_opened <- FALSE
  if (save && n_clusters <= 9) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    out_file <- file.path(folder, "cluster_histograms.png")
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(filename = out_file, width = 2400, height = 1800, res = 300)
    } else {
      png_type <- if (capabilities("cairo")) "cairo" else if (capabilities("X11")) "Xlib" else "cairo"
      grDevices::png(filename = out_file, width = 2400, height = 1800, res = 300, type = png_type)
    }
    device_opened <- TRUE
    on.exit(
      {
        if (device_opened) grDevices::dev.off()
      },
      add = TRUE
    )

    draw_histograms()

    grDevices::dev.off()
    device_opened <- FALSE
  }

  invisible(point_estimate)
}

plot_hist_cls_comuni <- function(results, BI, input_dir = "input/Comuni/", point_estimate = NULL, save = FALSE, folder = "results/plots/") {
  cat("Computing aggregated histograms of cluster assignments for Comuni...\n")

  # Load full_dataset.csv which contains histogram data (income brackets as columns)
  csv_file <- paste0(input_dir, "full_dataset.csv")
  if (!file.exists(csv_file)) {
    stop("Could not find full_dataset.csv in ", input_dir)
  }
  data <- read.csv(csv_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(data), "comuni\n")

  # Identify income bracket columns (columns starting with X followed by digits)
  all_cols <- colnames(data)
  income_columns <- all_cols[grep("^X[0-9]", all_cols)]
  
  # Remove X0. if X0.10000 exists (avoid double counting)
  if ("X0." %in% income_columns && "X0.10000" %in% income_columns) {
    income_columns <- income_columns[income_columns != "X0."]
  }
  
  if (length(income_columns) == 0) {
    stop("Could not detect income bracket columns in the data")
  }
  n_brackets <- length(income_columns)
  cat("Detected", n_brackets, "income brackets:", paste(income_columns, collapse = ", "), "\n")

  # Extract bin breaks from column names (lower bound of each bracket)
  bin_breaks_lower <- sapply(income_columns, function(col) {
    nums <- as.numeric(unlist(regmatches(col, gregexpr("[0-9]+", col))))
    if (length(nums) > 0) return(nums[1]) else return(NA)
  })
  
  # Create midpoints for each bracket (one per income column)
  bin_mids <- numeric(n_brackets)
  for (i in seq_len(n_brackets)) {
    lower <- bin_breaks_lower[i]
    if (i < n_brackets) {
      upper <- bin_breaks_lower[i + 1]
    } else {
      # Last bracket: use 1.5x the lower bound as approximate upper
      upper <- lower * 1.5
    }
    bin_mids[i] <- (lower + upper) / 2
  }
  
  # Create better labels for x-axis with k suffix for thousands
  bracket_labels <- sapply(seq_along(bin_breaks_lower), function(i) {
    lower <- bin_breaks_lower[i]
    if (is.na(lower)) return("?")
    
    # Format with k for thousands
    if (lower >= 1000) {
      lower_label <- paste0(round(lower/1000), "k")
    } else {
      lower_label <- as.character(lower)
    }
    
    # Add upper bound for clarity
    if (i < length(bin_breaks_lower)) {
      upper <- bin_breaks_lower[i + 1]
      if (upper >= 1000) {
        upper_label <- paste0(round(upper/1000), "k")
      } else {
        upper_label <- as.character(upper)
      }
      return(paste0(lower_label, "-", upper_label))
    } else {
      return(paste0(lower_label, "+"))
    }
  })

  if (is.null(point_estimate)) {
    point_estimate <- plot_cls_est(results, BI, save = FALSE)
  }

  unique_clusters <- sort(unique(point_estimate))
  n_clusters <- length(unique_clusters)

  cat("\n=== Cluster Diagnostics ===\n")
  cat("Number of clusters found:", n_clusters, "\n")
  cat("Cluster sizes:\n")
  print(table(point_estimate))

  # Extract income counts matrix
  income_matrix <- as.matrix(data[, income_columns])

  draw_histograms <- function() {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    # Get consistent colors for clusters
    cluster_colors <- get_cluster_colors(n_clusters)

    # Split into batches of max 9 clusters per plot
    max_plots_per_page <- 9
    n_batches <- ceiling(n_clusters / max_plots_per_page)

    for (batch_idx in seq_len(n_batches)) {
      # Determine which clusters go in this batch
      start_idx <- (batch_idx - 1) * max_plots_per_page + 1
      end_idx <- min(batch_idx * max_plots_per_page, n_clusters)
      batch_clusters <- unique_clusters[start_idx:end_idx]
      n_plots_in_batch <- length(batch_clusters)

      # Calculate layout for this batch
      n_rows <- ceiling(sqrt(n_plots_in_batch))
      n_cols <- ceiling(n_plots_in_batch / n_rows)
      # Increase bottom margin for better x-axis label spacing
      graphics::par(mfrow = c(n_rows, n_cols), mar = c(6, 4, 3, 1))

      for (cl in batch_clusters) {
        # Get comuni indices in this cluster
        cluster_idx <- which(point_estimate == cl)
        n_comuni <- length(cluster_idx)
        
        # Aggregate counts across all comuni in this cluster
        aggregated_counts <- colSums(income_matrix[cluster_idx, , drop = FALSE], na.rm = TRUE)
        total_count <- sum(aggregated_counts)
        
        # Compute density (normalized counts)
        density_vals <- aggregated_counts / total_count
        
        # Compute weighted mean income for this cluster
        mean_income <- sum(bin_mids * aggregated_counts) / total_count
        
        # Use consistent cluster colors with transparency
        cl_color <- adjustcolor(cluster_colors[as.character(cl)], alpha.f = 0.7)
        
        # Create barplot (histogram-style)
        bp <- barplot(density_vals,
          names.arg = NULL,
          main = paste("Cluster", cl, "\n(n =", n_comuni, "comuni)"),
          xlab = "",
          ylab = "Density",
          col = cl_color,
          border = "white",
          space = 0,
          xaxt = "n"  # Suppress default x-axis
        )
        
        # Add x-axis labels with better spacing
        # Show every label if <= 8 brackets, otherwise show every other label
        if (n_brackets <= 8) {
          label_indices <- seq_along(bracket_labels)
        } else {
          label_indices <- seq(1, n_brackets, by = 2)
        }
        
        axis(1, at = bp[label_indices], labels = bracket_labels[label_indices], 
             las = 2, cex.axis = 0.65, padj = 0.5)
        mtext("Income Bracket (€)", side = 1, line = 4.5, cex = 0.7)

        cat("\nCluster", cl, "statistics:\n")
        cat("  N comuni:", n_comuni, "\n")
        cat("  Total taxpayers (N):", total_count, "\n")
        cat("  Weighted mean income:", round(mean_income, 0), "\n")

        legend("topright",
          legend = c(
            paste("Mean:", format(round(mean_income, 0), big.mark = ",")),
            paste("N:", format(total_count, big.mark = ","))
          ),
          bty = "n",
          cex = 0.8
        )
      }

      # If saving and multiple batches, save each batch separately
      if (save && n_batches > 1) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        out_file <- file.path(folder, paste0("cluster_histograms_batch_", batch_idx, ".png"))
        if (requireNamespace("ragg", quietly = TRUE)) {
          ragg::agg_png(filename = out_file, width = 2400, height = 1800, res = 300)
        } else {
          png_type <- if (capabilities("cairo")) "cairo" else if (capabilities("X11")) "Xlib" else "cairo"
          grDevices::png(filename = out_file, width = 2400, height = 1800, res = 300, type = png_type)
        }

        # Redraw the plots for this batch
        graphics::par(mfrow = c(n_rows, n_cols), mar = c(6, 4, 3, 1))
        for (cl in batch_clusters) {
          cluster_idx <- which(point_estimate == cl)
          n_comuni <- length(cluster_idx)
          aggregated_counts <- colSums(income_matrix[cluster_idx, , drop = FALSE], na.rm = TRUE)
          total_count <- sum(aggregated_counts)
          density_vals <- aggregated_counts / total_count
          mean_income <- sum(bin_mids * aggregated_counts) / total_count
          cl_color <- adjustcolor(cluster_colors[as.character(cl)], alpha.f = 0.7)
          
          bp <- barplot(density_vals,
            names.arg = NULL,
            main = paste("Cluster", cl, "\n(n =", n_comuni, "comuni)"),
            xlab = "",
            ylab = "Density",
            col = cl_color,
            border = "white",
            space = 0,
            xaxt = "n"
          )
          
          # Show labels based on number of brackets
          if (n_brackets <= 8) {
            label_indices <- seq_along(bracket_labels)
          } else {
            label_indices <- seq(1, n_brackets, by = 2)
          }
          
          axis(1, at = bp[label_indices], labels = bracket_labels[label_indices], 
               las = 2, cex.axis = 0.65, padj = 0.5)
          mtext("Income Bracket (€)", side = 1, line = 4.5, cex = 0.7)
          
          legend("topright",
            legend = c(
              paste("Mean:", format(round(mean_income, 0), big.mark = ",")),
              paste("N:", format(total_count, big.mark = ","))
            ),
            bty = "n",
            cex = 0.8
          )
        }

        grDevices::dev.off()
        cat("\nSaved batch", batch_idx, "of", n_batches, "to", out_file, "\n")
      }
    }
  }

  draw_histograms()

  # Only save here if we have <= 9 clusters (single file case)
  device_opened <- FALSE
  if (save && n_clusters <= 9) {
    if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
    out_file <- file.path(folder, "cluster_histograms.png")
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(filename = out_file, width = 2400, height = 1800, res = 300)
    } else {
      png_type <- if (capabilities("cairo")) "cairo" else if (capabilities("X11")) "Xlib" else "cairo"
      grDevices::png(filename = out_file, width = 2400, height = 1800, res = 300, type = png_type)
    }
    device_opened <- TRUE
    on.exit(
      {
        if (device_opened) grDevices::dev.off()
      },
      add = TRUE
    )

    draw_histograms()

    grDevices::dev.off()
    device_opened <- FALSE
  }

  invisible(point_estimate)
}
