library(Rcpp)
library(ggplot2)
library(MASS)  # For fitdistr function
library(dplyr)
library(tidyr)

# Function to set hyperparameters based on elbow method and k-means
set_hyperparameters <- function(dist_matrix, max_k = 10, plot_results = FALSE) {
  # Step 1: Compute elbow method for K selection
  wss <- numeric(max_k)
  # For distance matrices, we need to use multidimensional scaling to get coordinates
  mds_result <- cmdscale(dist_matrix, k = 2)  # 2D embedding
  
  for (k in 1:max_k) {
    kmeans_result <- kmeans(mds_result, centers = k, nstart = 25)
    wss[k] <- kmeans_result$tot.withinss
  }
  
  # Find elbow using the "elbow method" - look for the point where improvement slows down
  # Calculate second derivatives to find the elbow
  if (max_k >= 3) {
    second_deriv <- numeric(max_k - 2)
    for (i in 2:(max_k - 1)) {
      second_deriv[i - 1] <- (wss[i + 1] - wss[i]) - (wss[i] - wss[i - 1])
    }
    k_elbow <- which.max(second_deriv) + 1
  } else {
    k_elbow <- 2
  }
  
  cat("Elbow method suggests K =", k_elbow, "\n")
  
  # Plot the elbow curve (optional)
  if (plot_results) {
    elbow_data <- data.frame(K = 1:max_k, WSS = wss)
    elbow_plot <- ggplot(elbow_data, aes(x = K, y = WSS)) +
      geom_line() +
      geom_point() +
      geom_vline(xintercept = k_elbow, color = "red", linetype = "dashed") +
      labs(title = "Elbow Method for Optimal K",
           x = "Number of Clusters (K)",
           y = "Within-Cluster Sum of Squares") +
      theme_minimal() +
      annotate("text", x = k_elbow + 0.5, y = max(wss) * 0.8, 
               label = paste("K =", k_elbow), color = "red")  
    print(elbow_plot)
  }
  
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
  if (plot_results) {
    cluster_data <- data.frame(
      x = mds_result[, 1],
      y = mds_result[, 2],
      cluster = as.factor(initial_clusters)
    )
    
    cluster_plot <- ggplot(cluster_data, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 3) +
      labs(title = paste("Initial K-means Clustering (K =", k_elbow, ")"),
           x = "MDS Dimension 1",
           y = "MDS Dimension 2") +
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
      tryCatch({
        gamma_fit_A <- fitdistr(A, "gamma")
        delta1 <- gamma_fit_A$estimate["shape"]
        
        # Step 5: Set alpha and beta
        n_A <- length(A)
        alpha <- delta1 * n_A
        beta <- sum(A)
        
        cat("For within-cluster distances:\n")
        cat("  delta1 (shape) =", delta1, "\n")
        cat("  alpha =", alpha, "\n")
        cat("  beta =", beta, "\n")
      }, error = function(e) {
        cat("Warning: Could not fit gamma distribution to within-cluster distances\n")
        cat("Error:", e$message, "\n")
        delta1 <<- 0.5
        alpha <<- 2
        beta <<- 2
      })
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
      tryCatch({
        gamma_fit_B <- fitdistr(B, "gamma")
        delta2 <- gamma_fit_B$estimate["shape"]
        
        # Set zeta and gamma (note: gamma is a parameter name, using different variable)
        n_B <- length(B)
        zeta <- delta2 * n_B
        gamma_param <- sum(B)  # renamed to avoid conflict with gamma() function
        
        cat("For inter-cluster distances:\n")
        cat("  delta2 (shape) =", delta2, "\n")
        cat("  zeta =", zeta, "\n")
        cat("  gamma =", gamma_param, "\n")
      }, error = function(e) {
        cat("Warning: Could not fit gamma distribution to inter-cluster distances\n")
        cat("Error:", e$message, "\n")
        delta2 <<- 2
        zeta <<- 2
        gamma_param <<- 2
      })
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

# Generate data from 3 Gaussian distributions
set.seed(123)
n_points <- 100
data1 <- matrix(rnorm(n_points * 2, mean = c(10, 10), sd = 1), ncol = 2)
data2 <- matrix(rnorm(n_points * 2, mean = c(1, 1), sd = 1), ncol = 2)
data3 <- matrix(rnorm(n_points * 2, mean = c(5, 5), sd = 2), ncol = 2)

# Combine all data points
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

# Set hyperparameters using the new method (disable plotting for now)
cat("Setting hyperparameters using elbow method and k-means clustering...\n")
hyperparams <- set_hyperparameters(dist_matrix,
                                   max_k = 10,
                                   plot_results = FALSE)

# Create parameter object with computed hyperparameters
param <- new(Params,
             hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
             hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
             0, 1000, 3, 1.0, 1.0)  # BI, NI, a, sigma, tau

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

### First plot - Posterior distribution of the number of clusters
post_k <- table(unlist(results$K)) / length(unlist(results$K))
df <- data.frame(cluster_found = as.numeric(names(post_k)),
                 rel_freq = as.numeric(post_k))

p1 <- ggplot(data = df, aes(x = factor(cluster_found), y = rel_freq)) +
  geom_col() +
  labs(
    x = "Cluster Found",
    y = "Relative Frequency",
  ) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_line(color = "grey95")) +
  scale_x_discrete(drop = FALSE)  # Ensures all cluster_found values are shown
print(p1)

### Second plot - Trace of number of clusters
k_df <- data.frame(
  Iteration = seq_along(results$K),
  NumClusters = unlist(results$K)
)

k_df_long <- k_df %>%
  pivot_longer(cols = starts_with("NumClusters"),
               names_to = "variable",
               values_to = "value")

p2 <- ggplot(k_df_long, aes(x = Iteration, y = value)) +
  geom_line() +
  labs(
    x = "Iteration",
    y = "Number of clusters",
  ) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15),
        panel.background = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        panel.grid.minor = element_line(color = "grey95"))
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
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text = element_text(size = 12),
        panel.background = element_blank()) +
  coord_fixed()
print(p3)

