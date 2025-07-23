library(Rcpp)
library(ggplot2)

# Generate data from 3 Gaussian distributions
set.seed(123)
n_points <- 2
data1 <- matrix(rnorm(n_points * 2, mean = c(10, 10), sd = 1), ncol = 2)
data2 <- matrix(rnorm(n_points * 2, mean = c(1, 1), sd = 1), ncol = 2)
data3 <- matrix(rnorm(n_points * 2, mean = c(20, 20), sd = 1), ncol = 2)

# Combine all data points
all_data <- rbind(data1, data2, data3)

# Create distance matrix (n x n)
dist_matrix <- as.matrix(dist(all_data))

# Create cluster labels for plotting
cluster_labels <- c(rep("Cluster 1", n_points),
                    rep("Cluster 2", n_points),
                    rep("Cluster 3", n_points))

# Create data frame for plotting
plot_data <- data.frame(
  x = all_data[, 1],
  y = all_data[, 2],
  cluster = cluster_labels
)

# Plot the clusters
ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 3) +
  labs(title = "Three Gaussian Clusters",
       x = "X coordinate",
       y = "Y coordinate") +
  theme_minimal()

# Load the C++ code
sourceCpp("src/launcher.cpp")
mcmc(dist_matrix)
