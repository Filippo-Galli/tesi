source("R/utils.R")
source("R/utils_plot.R")

##############################################################################
# Load data from 'input/Comuni' directory ====
##############################################################################

# Check if shapefile exists and load it
library(sf)

# Try to load the shapefile
tryCatch(
  {
    comuni_sf <- st_read("input/Comuni/full_dataset.shp", quiet = TRUE)
    cat("Shapefile loaded successfully with", nrow(comuni_sf), "features\n")
  },
  error = function(e) {
    stop("Error loading shapefile: ", e$message)
  }
)

# Check if .dat files exist, otherwise we'll need to create them
# from the shapefile data
if (!file.exists("input/Comuni/full_dataset.dat")) {
  cat("No .dat file found. Creating from shapefile...\n")

  # Extract numeric columns from the shapefile
  # Adjust column names based on your actual data
  numeric_cols <- sapply(st_drop_geometry(comuni_sf), is.numeric)

  if (sum(numeric_cols) == 0) {
    stop("No numeric columns found in shapefile data")
  }

  cat("Found", sum(numeric_cols), "numeric columns\n")
  cat("Column names:", names(comuni_sf)[numeric_cols], "\n")

  # Convert data to list format (one list element per feature)
  # Each element contains the numeric values for that feature
  data <- lapply(1:nrow(comuni_sf), function(i) {
    as.numeric(st_drop_geometry(comuni_sf)[i, numeric_cols])
  })

  # Create adjacency matrix based on spatial neighbors
  cat("Creating adjacency matrix based on spatial neighbors...\n")
  neighbors <- st_touches(comuni_sf)
  n <- nrow(comuni_sf)
  W <- matrix(0, n, n)

  for (i in 1:n) {
    if (length(neighbors[[i]]) > 0) {
      W[i, neighbors[[i]]] <- 1
    }
  }

  # Make symmetric
  W <- (W + t(W)) > 0
  W <- W * 1 # Convert logical to numeric

  # Save the data
  save(data, file = "input/Comuni/full_dataset.dat")
  save(W, file = "input/Comuni/adj_matrix.dat")

  cat("Data files created and saved\n")
} else {
  # Load existing .dat files
  load("input/Comuni/full_dataset.dat")
  load("input/Comuni/adj_matrix.dat")
  cat("Loaded existing .dat files\n")
}

##############################################################################
# Create histogram for each Comune ====
##############################################################################
hist_list <- list()

# Determine appropriate breaks for histograms
# Find global range across all data
all_data <- unlist(data)
global_min <- min(all_data, na.rm = TRUE)
global_max <- max(all_data, na.rm = TRUE)

# Create common breaks for all histograms
n_breaks <- 30
breaks <- seq(global_min, global_max, length.out = n_breaks + 1)

cat("Creating histograms for", length(data), "comuni\n")
cat("Using", n_breaks, "bins from", global_min, "to", global_max, "\n")

for (i in seq_along(data)) {
  comuni_data <- data[[i]]
  # Remove NA values
  comuni_data <- comuni_data[!is.na(comuni_data)]

  if (length(comuni_data) > 0) {
    hist_list[[i]] <- hist(comuni_data, breaks = breaks, plot = FALSE)
  } else {
    # Create empty histogram if no data
    hist_list[[i]] <- hist(numeric(0), breaks = breaks, plot = FALSE)
  }
}

##############################################################################
# Distances between histograms ====
##############################################################################
cat("Computing distances between histograms...\n")

distance_jeff_divergences <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)
distance_cm <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)
distance_wasserstein <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)
distance_mean <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)

for (i in seq_along(hist_list)) {
  if (i %% 10 == 0) {
    cat("Processing histogram", i, "of", length(hist_list), "\n")
  }
  for (j in seq_along(hist_list)) {
    distance_jeff_divergences[i, j] <-
      compute_hist_distances(hist_list[[i]],
        hist_list[[j]],
        type = "Jeff"
      )
    distance_cm[i, j] <-
      compute_hist_distances(hist_list[[i]],
        hist_list[[j]],
        type = "CM"
      )
    distance_wasserstein[i, j] <-
      compute_hist_distances(hist_list[[i]],
        hist_list[[j]],
        type = "Wasserstein"
      )

    distance_mean[i, j] <- abs(mean(data[[i]], na.rm = TRUE) -
      mean(data[[j]], na.rm = TRUE))
  }
}

cat("Distance computation completed\n")

##############################################################################
# Plot Distance (optional) ====
##############################################################################
if (exists("plot_distance")) {
  cat("Creating distance plots...\n")
  plot_distance(distance_jeff_divergences,
    title = "Jeffreys Divergence Distance", save = TRUE,
    folder = "results/distance_plots/Comuni/"
  )
  plot_distance(distance_cm,
    title = "Cramer-von Mises Distance", save = TRUE,
    folder = "results/distance_plots/Comuni/"
  )
  plot_distance(distance_wasserstein,
    title = "Wasserstein Distance", save = TRUE,
    folder = "results/distance_plots/Comuni/"
  )
  plot_distance(distance_mean,
    title = "Mean Difference", save = TRUE,
    folder = "results/distance_plots/Comuni/"
  )
}

##############################################################################
# Save Data ====
##############################################################################
cat("Saving results...\n")

# Create output folder if it doesn't exist
folder <- "real_data/Comuni"
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}

# Save distance matrices
saveRDS(distance_jeff_divergences,
  file = paste0(folder, "/distance_jeff_divergences.rds")
)
saveRDS(distance_cm,
  file = paste0(folder, "/distance_cm.rds")
)
saveRDS(distance_wasserstein,
  file = paste0(folder, "/distance_wasserstein.rds")
)
saveRDS(distance_mean,
  file = paste0(folder, "/distance_mean.rds")
)

# Save adjacency matrix
saveRDS(W, file = paste0(folder, "/adj_matrix.rds"))

cat("All results saved to", folder, "\n")
cat("\nSummary:\n")
cat("  - Number of comuni:", length(hist_list), "\n")
cat("  - Adjacency matrix dimension:", nrow(W), "x", ncol(W), "\n")
cat("  - Number of neighbors (edges):", sum(W) / 2, "\n")
cat("  - Distance matrices saved: 7 types\n")
