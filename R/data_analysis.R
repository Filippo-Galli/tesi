source("R/utils.R")
source("R/utils_plot.R")

##############################################################################
# Load .dat files from 'input/' directory ====
##############################################################################
load("input/full_dataset.dat")
load("input/adj_matrix.dat")

##############################################################################
# Create histogram for each PUMA ====
##############################################################################
histograms_list <- list()

for (i in seq_along(data)) {
  pumas_data <- data[[i]]
  histograms_list[[i]] <- hist(pumas_data, plot = FALSE,
                               breaks = 20, probability = TRUE)
}

##############################################################################
# Distances between histograms ====
##############################################################################
distance_hist_intersection <- matrix(0, nrow = length(histograms_list),
                                     ncol = length(histograms_list))
distance_jeff_divergences <- matrix(0, nrow = length(histograms_list),
                                    ncol = length(histograms_list))
distance_chi_squared <- matrix(0, nrow = length(histograms_list),
                               ncol = length(histograms_list))
distance_euclidean <- matrix(0, nrow = length(histograms_list),
                             ncol = length(histograms_list))
distance_cm <- matrix(0, nrow = length(histograms_list),
                      ncol = length(histograms_list))
distance_wasserstein <- matrix(0, nrow = length(histograms_list),
                               ncol = length(histograms_list))

distance_mean <- matrix(0, nrow = length(histograms_list),
                        ncol = length(histograms_list))

for (i in seq_along(histograms_list)) {
  for (j in seq_along(histograms_list)) {
    distance_hist_intersection[i, j] <-
      compute_hist_distances(histograms_list[[i]], histograms_list[[j]])

    distance_jeff_divergences[i, j] <-
      compute_hist_distances(histograms_list[[i]],
                             histograms_list[[j]],
                             type = "Jeff")
    distance_chi_squared[i, j] <-
      compute_hist_distances(histograms_list[[i]],
                             histograms_list[[j]],
                             type = "chi2")
    distance_euclidean[i, j] <-
      compute_hist_distances(histograms_list[[i]],
                             histograms_list[[j]],
                             type = "euclidean")
    distance_cm[i, j] <-
      compute_hist_distances(histograms_list[[i]],
                             histograms_list[[j]],
                             type = "CM")
    distance_wasserstein[i, j] <-
      compute_hist_distances(histograms_list[[i]],
                             histograms_list[[j]],
                             type = "Wasserstein")

    distance_mean[i, j] <- sqrt(abs(mean(data[[i]])^2 - mean(data[[j]])^2))
  }
}

##############################################################################
# Plot Distance ====
##############################################################################
plot_distance(distance_euclidean,
              title = "Euclidean Distance", save = TRUE,
              folder = "results/distance_plots/")
plot_distance(distance_hist_intersection,
              title = "Histogram Intersection Distance", save = TRUE,
              folder = "results/distance_plots/")
plot_distance(distance_jeff_divergences,
              title = "Jeffreys Divergence Distance", save = TRUE,
              folder = "results/distance_plots/")
plot_distance(distance_chi_squared,
              title = "Chi-Squared Distance", save = TRUE,
              folder = "results/distance_plots/")
plot_distance(distance_cm,
              title = "Cramer-von Mises Distance", save = TRUE,
              folder = "results/distance_plots/")
plot_distance(distance_wasserstein,
              title = "Wasserstein Distance", save = TRUE,
              folder = "results/distance_plots/")
plot_distance(distance_mean,
              title = "Mean Difference", save = TRUE,
              folder = "results/distance_plots/")

##############################################################################
# Save Data ====
##############################################################################

# # Save distance matrices
folder <- "real_data"
# saveRDS(distance_hist_intersection,
#         file = paste0(folder, "/distance_hist_intersection.rds"))
# saveRDS(distance_jeff_divergences,
#         file = paste0(folder, "/distance_jeff_divergences.rds"))
# saveRDS(distance_chi_squared,
#         file = paste0(folder, "/distance_chi_squared.rds"))
# saveRDS(distance_euclidean,
#         file = paste0(folder, "/distance_euclidean.rds"))
# saveRDS(distance_cm,
#         file = paste0(folder, "/distance_cm.rds"))
# saveRDS(distance_wasserstein,
#         file = paste0(folder, "/distance_wasserstein.rds"))
saveRDS(distance_mean,
        file = paste0(folder, "/distance_mean.rds"))

# Save adjacency matrix
saveRDS(W, file = paste0(folder, "/adj_matrix.rds"))
