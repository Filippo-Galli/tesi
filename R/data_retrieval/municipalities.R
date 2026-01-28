# Script to download and process Italian municipality data
# Simplified version of the original mechanism in spatialCMC by TeoGiane
# Author: Filippo Galli
# Modified to match naming conventions with ACF.R

# Load required libraries ====
library(dplyr)
library(sf)

# Create necessary directories ====
dir.create("raw", showWarnings = FALSE, recursive = TRUE)
dir.create("input/municipalities", showWarnings = FALSE, recursive = TRUE)

cat("Starting data download and processing...\n\n")

# Step 1: Download raw data ====
cat("Step 1: Downloading raw data files...\n")

download.file(
  url = "https://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip",
  destfile = "raw/Limiti01012020.zip",
  mode = "wb"
)

download.file(
  url = "https://www.istat.it/wp-content/uploads/2020/03/Dataset-decessi-comunali-giornalieri-e-tracciato-record_5marzo.zip",
  destfile = "raw/Dataset-decessi-comunali.zip",
  mode = "wb"
)

download.file(
  url = "https://demo.istat.it/data/posas/POSAS_2021_it_Comuni.zip",
  destfile = "raw/POSAS_2021.zip",
  mode = "wb"
)

download.file(
  url = "https://www1.finanze.gov.it/finanze/analisi_stat/public/v_4_0_0/contenuti/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.zip",
  destfile = "raw/Redditi_2020.zip",
  mode = "wb"
)

cat("Downloads complete!\n\n")

# Step 2: Unzip files ====
cat("Step 2: Extracting files...\n")

unzip("raw/Limiti01012020.zip", exdir = "raw")
unzip("raw/Dataset-decessi-comunali.zip", exdir = "raw")
unzip("raw/POSAS_2021.zip", exdir = "raw")
unzip("raw/Redditi_2020.zip", exdir = "raw")

cat("Extraction complete!\n\n")

# Step 3: Import shapefiles ====
cat("Step 3: Importing geographic data (shapefiles)...\n")

mun_sf <- read_sf("raw/Limiti01012020/Com01012020/Com01012020_WGS84.shp")

cat("Imported", nrow(mun_sf), "municipalities\n\n")

# Step 4: Merge Monteciccardo into Pesaro ====
cat("Step 4: Merging Monteciccardo into Pesaro municipality...\n")

merged_mun <- mun_sf %>%
  filter(COMUNE %in% c("Monteciccardo", "Pesaro")) %>%
  summarise(
    COD_RIP = 3, COD_REG = 11, COD_PROV = 41, COD_CM = 0, COD_UTS = 41,
    PRO_COM = 41044, PRO_COM_T = "041044", COMUNE = "Pesaro", COMUNE_A = NA, CC_UTS = 1,
    SHAPE_LENG = as.numeric(st_length(st_boundary(st_union(geometry)))),
    SHAPE_AREA = as.numeric(st_area(st_union(geometry))),
    geometry = st_union(geometry)
  )

mun_sf <- mun_sf %>%
  filter(!COMUNE %in% c("Monteciccardo", "Pesaro")) %>%
  bind_rows(merged_mun)

cat("Municipalities after merge:", nrow(mun_sf), "\n\n")

# Step 5: Import and clean income data ====
cat("Step 5: Processing income data...\n")

# Read column names first
income_cols <- names(read.csv(
  "raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.csv",
  sep = ";", header = TRUE, nrows = 0
))

# Read actual data
income_df_raw <- read.csv(
  "raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.csv",
  sep = ";", skip = 1, header = FALSE
) %>%
  select(-c(ncol(.) - 1, ncol(.))) %>%  # Remove last 2 empty columns
  `colnames<-`(income_cols) %>%
  filter(Regione != "Mancante/errata")

# Extract income distribution by brackets
income_distribution <- income_df_raw %>%
  select(
    Codice.Istat.Comune, Denominazione.Comune,
    matches("Reddito.complessivo.*Frequenza")
  ) %>%
  rename_with(
    ~ gsub("Reddito.complessivo.", "", .x) %>%
      gsub("\\.+Frequenza", "", .) %>%
      gsub("\\.", " ", .) %>%
      gsub("da ([0-9]+) a ([0-9]+)", "\\1-\\2", .) %>%
      gsub("oltre ([0-9]+)", "\\1+", .) %>%
      gsub("minore o uguale a zero", "0-", .) %>%
      gsub(" euro", "", .),
    matches("Frequenza")
  ) %>%
  mutate(across(matches("^[0-9]|^0-"), ~ replace(., is.na(.), 0))) %>%
  rename(COD_MUN = Codice.Istat.Comune, NAME_MUN = Denominazione.Comune)

cat("Processed income data for", nrow(income_distribution), "municipalities\n\n")

# Step 6: Prepare geometry reference ====
cat("Step 6: Preparing geometry reference...\n")

geometry_sf <- mun_sf %>%
  select(PRO_COM, COD_REG, COD_PROV, COMUNE, geometry) %>%
  rename(COD_MUN = PRO_COM, NAME_MUN = COMUNE)

cat("Geometry prepared\n\n")

# Step 7: Compute adjacency matrix ====
cat("Step 7: Computing spatial adjacency matrix...\n")

# Get neighbors using queen contiguity
neighbors <- st_touches(geometry_sf, sparse = FALSE)
adj_matrix <- matrix(as.integer(neighbors), nrow = nrow(geometry_sf))
rownames(adj_matrix) <- geometry_sf$COD_MUN
colnames(adj_matrix) <- geometry_sf$COD_MUN

# Check symmetry
if (!isSymmetric(adj_matrix)) {
  warning("Adjacency matrix is not symmetric!")
} else {
  cat("Adjacency matrix is symmetric ✓\n")
}

cat("Adjacency matrix:", nrow(adj_matrix), "x", ncol(adj_matrix), "\n\n")

# Step 8: Save outputs ====
cat("Step 8: Saving final datasets...\n")

# Save income distribution data
write.csv(income_distribution, "input/municipalities/full_dataset.csv", row.names = FALSE)
cat("  ✓ Income data saved\n")

# Save geometry as RDS
saveRDS(geometry_sf, "input/municipalities/geometry.rds")
cat("  ✓ Geometry RDS saved\n")

# Save geometry as shapefile in geometry folder
geometry_dir <- "input/municipalities/geometry"
dir.create(geometry_dir, showWarnings = FALSE, recursive = TRUE)
st_write(geometry_sf, file.path(geometry_dir, "municipalities.shp"), 
         append = FALSE, quiet = TRUE)
cat("  ✓ Shapefile saved\n")

# Save adjacency matrix
saveRDS(adj_matrix, "input/municipalities/adj_matrix.rds")
cat("  ✓ Adjacency matrix saved\n")

cat("\n✓ All processing complete!\n\n")

# Summary ====
cat("=== OUTPUT FILES ===\n")
cat("1. input/municipalities/full_dataset.csv - Income distribution data (", 
    nrow(income_distribution), " municipalities)\n", sep = "")
cat("2. input/municipalities/geometry.rds - Geographic boundaries (", 
    nrow(geometry_sf), " municipalities)\n", sep = "")
cat("3. input/municipalities/geometry/municipalities.shp - Shapefile boundaries\n")
cat("4. input/municipalities/adj_matrix.rds - Spatial adjacency matrix (", 
    nrow(adj_matrix), "x", ncol(adj_matrix), ")\n\n", sep = "")

cat("=== USAGE ===\n")
cat("# Load data:\n")
cat("data <- read.csv('input/municipalities/full_dataset.csv')\n")
cat("geo <- readRDS('input/municipalities/geometry.rds')\n")
cat("adj <- readRDS('input/municipalities/adj_matrix.rds')\n\n")
cat("# Or load shapefile:\n")
cat("geo <- st_read('input/municipalities/geometry/municipalities.shp')\n\n")
cat("# Join for spatial analysis:\n")
cat("full <- left_join(geo, data, by = 'COD_MUN')\n")

# Optional: Clean up raw files (uncomment if desired)
# unlink("raw", recursive = TRUE)
# cat("\nRaw files cleaned up.\n")