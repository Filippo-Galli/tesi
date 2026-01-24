##############################################################################
# Setup ====
##############################################################################

source("R/utils_plot.R")
#library(flexclust)
library(cluster)
# library(NbClust)  # optional, uncomment for nbclust()

##############################################################################
# Data Loading ====
##############################################################################

files_folder <- "real_data/LA"
files <- list.files(files_folder)
file_chosen <- files[3]
dist_matrix <- readRDS(file = file.path(files_folder, file_chosen))
puma_age <- readRDS(file = file.path(files_folder, "puma_agep_std_mean.rds"))
puma_sex <- readRDS(file = file.path(files_folder, "puma_sex_mode.rds"))

##############################################################################
# k Selection: Gap Statistic (fixed for dist) ====
##############################################################################

# Convert dist to full symmetric dissimilarity matrix (n x n)
n <- attr(dist_matrix, "Size")
diss_matrix <- as.matrix(dist_matrix)

# Now clusGap works (treats as features, but gap computed correctly)
set.seed(456)
gap_res <- clusGap(diss_matrix, FUN = pam, K.max = 8, B = 20,  # B=20 for speed; use 500+
                   diss = inherits(diss_matrix, "dist"))  # diss=FALSE since matrix
plot(gap_res)
print(gap_res, method = "firstmax")

# Alternative: NbClust (direct dist support, many indices)
# nb <- NbClust(dist_matrix, diss = NULL, method = "ward.D2")  # or "complete", etc.

##############################################################################
# PAM Clustering ====
##############################################################################

# Best k
k <- 7

if (!inherits(dist_matrix, "dist")) {
  dist_matrix <- as.dist(dist_matrix)
}

set.seed(123)
pam_cl <- pam(dist_matrix, k = k, diss = TRUE)
cluster_id <- pam_cl$clustering

##############################################################################
# Save Results ====
##############################################################################

results <- list(
  clusters     = cluster_id,
  pam_model    = pam_cl,
  medoids      = medoids,
  sil_width    = pam_cl$silinfo$avg.width,
  dist_file    = file_chosen,
  puma_age     = puma_age,
  puma_sex     = puma_sex
)

file_chosen <- sub("\\.rds$", "", file_chosen)
files_folder_clean <- gsub("/", "_", files_folder)
output_folder <- paste0("results/", files_folder_clean, "_k-pam_k", k)
if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}
saveRDS(results, paste0(output_folder, "/results.rds"))
cat("Saved to:", output_folder, "\n")

##############################################################################
# Diagnostics ====
##############################################################################

print(pam_cl)
plot(silhouette(pam_cl))
table(cluster_id)
medoids <- pam_cl$medoids
cat("Medoid indices:", medoids, "\n")

plot_hist_cls_pumas(results = NULL, BI = 0, point_estimate = cluster_id, save = TRUE, folder = output_folder )

id_col <- "PUMA"
states <- "LA"
puma_ids <- sf::st_read(paste0("input/", states, "/counties-pumas/counties-pumas.shp"), quiet = TRUE)[[id_col]]

plot_map_cls(
  results = NULL,
  BI = 0,
  point_estimate = cluster_id,
  unit_ids = puma_ids,
  puma_dir = paste0("input/", states, "/counties-pumas"),
  id_col = id_col,
  save = TRUE, 
  folder = output_folder
)
