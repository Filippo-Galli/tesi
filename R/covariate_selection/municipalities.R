source("R/utils_plot.R")
source("R/utils.R")

# Load required packages
library(dplyr)

##############################################################################
# CONFIGURATION ====
##############################################################################

# Configuration - modify these as needed
CONFIG <- list(
    input_dir = "input/municipalities",
    output_dir = "real_data/municipalities",
    results_dir = "results",
    # Choose format: "rds" (list format, recommended) or "csv" (long format)
    covariates_format = "csv",
    # Minimum proportion of non-NA values to keep a column (0 = keep all, 1 = no NA allowed)
    na_threshold = 1
)

# Create output directory if it doesn't exist
if (!dir.exists(CONFIG$output_dir)) {
    dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
}

##############################################################################
# LOAD DATA ====
##############################################################################

cat("=== Loading Covariates Data ===\n")

# Load covariates - supports both RDS (list format, recommended) and CSV (long format)
if (CONFIG$covariates_format == "rds") {
    covariates_file <- file.path(CONFIG$input_dir, "full_dataset_covariates.rds")
    if (!file.exists(covariates_file)) {
        stop(sprintf("File not found: %s. Try setting covariates_format = 'csv' or run municipalities.R first.", covariates_file))
    }

    # Load as list (one dataframe per MUN)
    full_dataset_covariates <- readRDS(covariates_file)
    cat(sprintf("Loaded RDS format: %d MUNs\n", length(full_dataset_covariates)))

    # Combine into single dataframe for processing, preserving COD_MUN
    all_covariates <- do.call(rbind, lapply(names(full_dataset_covariates), function(MUN_id) {
        df <- full_dataset_covariates[[MUN_id]]
        df$COD_MUN <- MUN_id
        return(df)
    }))
} else {
    covariates_file <- file.path(CONFIG$input_dir, "full_dataset_covariates.csv")
    if (!file.exists(covariates_file)) {
        stop(sprintf("File not found: %s", covariates_file))
    }

    # Load long format CSV
    all_covariates <- read.csv(covariates_file, stringsAsFactors = FALSE)
    cat(sprintf("Loaded CSV format: %d rows, %d columns\n", nrow(all_covariates), ncol(all_covariates)))

    # Convert to list format for compatibility with existing code
    MUN_ids <- unique(all_covariates$COD_MUN)
    full_dataset_covariates <- lapply(MUN_ids, function(id) {
        subset(all_covariates, COD_MUN == id)
    })
    names(full_dataset_covariates) <- MUN_ids
}

# Ensure COD_MUN exists (for merging with spatial data later)
if (!"COD_MUN" %in% colnames(all_covariates)) {
    if ("STATE_MUN" %in% colnames(all_covariates)) {
        all_covariates$COD_MUN <- all_covariates$STATE_MUN
    } else {
        stop("Neither COD_MUN nor STATE_MUN column found in data")
    }
}

print("=== Checking Municipality Population Totals ===")
problem <- FALSE

# Group by COD_MUN and check totals for each municipality
mun_groups <- all_covariates %>%
    group_by(COD_MUN) %>%
    summarise(
        totale_femmine = sum(TOTALE_FEMMINE, na.rm = TRUE),
        totale_maschi = sum(TOTALE_MASCHI, na.rm = TRUE),
        totale = sum(TOTALE, na.rm = TRUE),
        .groups = "drop"
    )

# Check if totals match for each municipality
inconsistent_muns <- mun_groups %>%
    filter(totale_femmine + totale_maschi != totale)

if (nrow(inconsistent_muns) > 0) {
    cat(sprintf("Found %d municipalities with inconsistent totals:\n", nrow(inconsistent_muns)))
    print(inconsistent_muns)
    problem <- TRUE
} else {
    cat("All municipality totals are consistent.\n")
}

##############################################################################
# Inspect DATA ====
##############################################################################

cat("=== Analysing NA Covariates Data ===\n")

# Get initial column info
initial_cols <- ncol(all_covariates)
cat(sprintf("Initial columns: %d\n", initial_cols))

# 1. Handle NA values based on threshold
if (CONFIG$na_threshold > 0) {
    na_proportions <- sapply(all_covariates, function(col) mean(is.na(col)))
    cols_to_keep <- names(na_proportions)[na_proportions <= (1 - CONFIG$na_threshold)]

    removed_cols <- setdiff(names(all_covariates), cols_to_keep)
    if (length(removed_cols) > 0) {
        cat(sprintf(
            "Removing %d columns with >%.0f%% NA values\n",
            length(removed_cols), (1 - CONFIG$na_threshold) * 100
        ))
    }

    all_covariates <- all_covariates[, cols_to_keep, drop = FALSE]
}

final_cols <- ncol(all_covariates)
cat(sprintf("Final columns after NA filtering: %d (removed %d)\n", final_cols, initial_cols - final_cols))
print("Remaining columns:")
print(colnames(all_covariates))

##############################################################################
# Percentage of women (continuous) ====
##############################################################################

cat("=== Calculating Percentage of Women ===\n")

mun_sex_stats <- all_covariates %>%
    group_by(COD_MUN) %>%
    summarise(
        Women_perc = mean(TOTALE_FEMMINE / TOTALE, na.rm = TRUE),
        number_women = sum(TOTALE_FEMMINE, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(COD_MUN)

saveRDS(mun_sex_stats, file = file.path(CONFIG$output_dir, "mun_sex_stats.rds"))
write.csv(mun_sex_stats,
    file = file.path(CONFIG$output_dir, "mun_sex_stats.csv"),
    row.names = FALSE
)
cat(sprintf("  Saved sex statistics for %d municipalities\n", nrow(mun_sex_stats)))

##############################################################################
# Women Marital Status (categorical) ====
##############################################################################

cat("=== Calculating Marital Status Proportions ===\n")
mun_marital_stats <- all_covariates %>%
    group_by(COD_MUN) %>%
    summarise(
        Married_perc = mean(CONIUGATE, na.rm = TRUE),
        Single_perc = mean(NUBILI, na.rm = TRUE),
        Divorced_perc = mean(DIVORZIATE, na.rm = TRUE),
        Widowed_perc = mean(VEDOVE, na.rm = TRUE),
        number_women = sum(TOTALE_FEMMINE, na.rm = TRUE),
        majority_status = case_when(
            Married_perc >= max(Single_perc, Divorced_perc, Widowed_perc) ~ "Married",
            Single_perc >= max(Married_perc, Divorced_perc, Widowed_perc) ~ "Single",
            Divorced_perc >= max(Married_perc, Single_perc, Widowed_perc) ~ "Divorced",
            Widowed_perc >= max(Married_perc, Single_perc, Divorced_perc) ~ "Widowed",
            TRUE ~ NA_character_
        ),
        .groups = "drop"
    ) %>%
    arrange(COD_MUN)

head(mun_marital_stats)
saveRDS(mun_marital_stats, file = file.path(CONFIG$output_dir, "mun_marital_stats.rds"))
write.csv(mun_marital_stats,
    file = file.path(CONFIG$output_dir, "mun_marital_stats.csv"),
    row.names = FALSE
)
cat(sprintf("  Saved marital status statistics for %d municipalities\n", nrow(mun_marital_stats)))
