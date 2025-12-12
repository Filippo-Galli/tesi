source("R/utils_plot.R")
source("R/utils.R")

##############################################################################
# Load .dat files from 'input/' directory ====
##############################################################################
load("input/LA/full_dataset.dat") # it load the variable 'data'
load("input/LA/adj_matrix.dat")
full_dataset_covariates <- readRDS("input/LA/full_dataset_covariates.rds")

##############################################################################
# Covariates cleaning and preprocessing ====
##############################################################################

covariate_names <- colnames(full_dataset_covariates[[1]])

pums_variables <- list(
  
  # ========== CONTINUOUS VARIABLES (Numeric amounts, counts, time, years) ==========
  continuous = c(
    # Adjustment factors
    "ADJHSG",      # Adjustment factor for housing dollar amounts [cite: 69]
    "ADJINC",      # Adjustment factor for income amounts [cite: 75]
    
    # Income and earnings (dollar amounts)
    "FINCP",       # Family income (past 12 months) [cite: 388]
    "HINCP",       # Household income (past 12 months) [cite: 422]
    "INTP",        # Interest, dividends, net rental income [cite: 1490]
    "OIP",         # All other income (past 12 months) [cite: 1720]
    "PAP",         # Public assistance income [cite: 1742]
    "PERNP",       # Person's earnings [cite: 4014]
    "PINCP",       # Total person's income [cite: 4027]
    "RETP",        # Retirement income [cite: 1761]
    "SEMP",        # Self-employment income [cite: 1814]
    "SSIP",        # Supplementary Security Income [cite: 1846]
    "SSP",         # Social Security income [cite: 1853]
    "WAGP",        # Wages or salary income [cite: 1860]
    
    # Housing costs (dollar amounts)
    "CONP",        # Condo fee (monthly) [cite: 137]
    "ELEP",        # Electricity (monthly cost) [cite: 161]
    "FULP",        # Fuel cost (yearly) [cite: 170]
    "GASP",        # Gas (monthly cost) [cite: 182]
    "GRNTP",       # Gross rent (monthly) [cite: 398]
    "INSP",        # Fire/hazard/flood insurance (yearly) [cite: 204]
    "MHP",         # Mobile home costs (yearly) [cite: 210]
    "MRGP",        # First mortgage payment (monthly) [cite: 218]
    "RNTP",        # Monthly rent [cite: 258]
    "SMP",         # Total payment on mortgages/home equity loans [cite: 269]
    "SMOCP",       # Selected monthly owner costs [cite: 569]
    "VALP",        # Property value [cite: 323]
    "WATP",        # Water (yearly cost) [cite: 342]
    
    # Time and Years
    "JWMNP",       # Travel time to work (minutes) [cite: 1500]
    "WKHP",        # Usual hours worked per week [cite: 1866]
    "AGEP",        # Age [cite: 1365]
    "CITWP",       # Year of naturalization [cite: 1379]
    "YOEP",        # Year of entry [cite: 1875]
    
    # Counts
    "BDSP",        # Number of bedrooms [cite: 126]
    "NOC",         # Number of own children [cite: 486]
    "NRC",         # Number of related children [cite: 490]
    "NPF",         # Number of persons in family [cite: 487]
    "RMSP",        # Number of rooms [cite: 234]
    "RACNUM"       # Number of major race groups represented [cite: 4349]
  ),
  
  # ========== CATEGORICAL NOMINAL VARIABLES (Unordered categories/Codes) ==========
  categorical_nominal = c(
    # Geography & Record IDs
    "RT",          # Record type [cite: 1]
    "ST",          # State code [cite: 34]
    "PUMA",        # Public use microdata area [cite: 12]
    "DIVISION",    # Division code [cite: 11]
    "REGION",      # Region code [cite: 24]
    "POBP",        # Place of birth (Recode) [cite: 4039]
    "POWPUMA",     # Place of work PUMA [cite: 4236]
    "POWSP",       # Place of work State/Country [cite: 4242]
    "MIGPUMA",     # Migration PUMA [cite: 2989]
    "MIGSP",       # Migration recode (State/Country) [cite: 3048]
    "WAOB",        # World area of birth [cite: 4398]
    
    # Household & Family Composition
    "HHT",         # Household/family type [cite: 412]
    "HHL",         # Household language [cite: 409]
    "RELP",        # Relationship [cite: 1749]
    "SFN",         # Subfamily number [cite: 4353]
    "SFR",         # Subfamily relationship [cite: 4353]
    "MSP",         # Married status (spouse present/absent) [cite: 3139]
    "PAOC",        # Presence and age of own children [cite: 4006]
    "QTRBIR",      # Quarter of birth [cite: 4305]
    
    # Housing Characteristics
    "TYPE",        # Type of unit [cite: 94]
    "HFL",         # House heating fuel [cite: 191]
    "RESMODE",     # Response mode [cite: 568]
    "TEN",         # Tenure [cite: 298]
    "VACS",        # Vacancy status [cite: 302]
    
    # Person Demographics & Origins
    "SEX",         # Sex [cite: 1841]
    "MAR",         # Marital status [cite: 1540]
    "CIT",         # Citizenship status [cite: 1368]
    "MIG",         # Mobility status (lived here 1 year ago) [cite: 1598]
    "NATIVITY",    # Nativity [cite: 3519]
    "NOP",         # Nativity of parent [cite: 3524]
    "ANC",         # Ancestry recode [cite: 1910]
    "ANC1P",       # Recoded Detailed Ancestry first entry [cite: 1937]
    "ANC2P",       # Recoded Detailed Ancestry second entry [cite: 2136]
    "HISP",        # Recoded detailed Hispanic origin [cite: 2657]
    "RAC1P",       # Recoded detailed race code [cite: 4315]
    "RAC2P",       # Recoded detailed race code [cite: 4316]
    "RAC3P",       # Recoded detailed race code [cite: 4321]
    "LANP",        # Language spoken at home [cite: 2932]
    
    # Education & Employment Codes
    "SCH",         # School enrollment [cite: 1780]
    "COW",         # Class of worker [cite: 1387]
    "ESR",         # Employment status recode [cite: 2333]
    "INDP",        # Industry recode [cite: 2687]
    "NAICSP",      # NAICS Industry code [cite: 3139]
    "OCCP",        # Occupation recode 2010 [cite: 3538]
    "SOCP",        # SOC Occupation code [cite: 4353]
    "FOD1P",       # Field of degree first entry [cite: 2350]
    "FOD2P",       # Field of degree second entry [cite: 2523]
    "JWTR",        # Means of transportation to work [cite: 1505]
    "JWAP",        # Time of arrival at work [cite: 2903]
    "JWDP",        # Time of departure for work [cite: 2924]
    "MIL",         # Military service [cite: 1599]
    "VPS",         # Veteran period of service [cite: 4396]
    
    # Allocation Flags (Categorical)
    "FHINS3C", "FHINS4C", "FHINS5C" # Coverage edits [cite: 4403]
  ),
  
  # ========== CATEGORICAL ORDINAL VARIABLES (Ordered levels) ==========
  categorical_ordinal = c(
    # Housing & Vehicles
    "ACR",         # Lot size [cite: 105]
    "BLD",         # Units in structure [cite: 135]
    "YBL",         # When structure first built [cite: 349]
    "MV",          # When moved into house or apartment [cite: 485]
    "VEH",         # Vehicles available [cite: 330]
    "TAXP",        # Property taxes (yearly) - Ordinal codes ($1-$49, etc.) [cite: 589]
    
    # Income/Cost Ratios & Ranges
    "AGS",         # Sales of agriculture products (ranges) [cite: 107]
    "GRPIP",       # Gross rent as percentage of income [cite: 404]
    "OCPIP",       # Owner costs as percentage of income [cite: 525]
    "POVPIP",      # Income-to-poverty ratio recode [cite: 4226]
    
    # Education & Language
    "SCHG",        # Grade level attending [cite: 1781]
    "SCHL",        # Educational attainment [cite: 1782]
    "ENG",         # Ability to speak English [cite: 1416]
    
    # Household Composition (Ordered Counts/Presence)
    "NP",          # Number of person records [cite: 85]
    "R18",         # Presence of persons under 18 [cite: 536]
    "R60",         # Presence of persons 60 and over [cite: 554]
    "R65",         # Presence of persons 65 and over [cite: 555]
    "FPARC",       # Family presence and age of related children [cite: 392]
    "HUPAC",       # HH presence and age of children [cite: 446]
    "HUPAOC",      # HH presence and age of own children [cite: 457]
    "HUPARC",      # HH presence and age of related children [cite: 460]
    "GCM",         # Length of time responsible for grandchildren [cite: 1443]
    "MARHT",       # Number of times married [cite: 1561]
    
    # Work & Commute
    "WKW",         # Weeks worked (ranges) [cite: 1873]
    "WKL",         # When last worked [cite: 1872]
    "WKEXREL",     # Work experience of householder and spouse [cite: 748]
    "WORKSTAT",    # Work status of householder or spouse [cite: 762]
    "WIF",         # Workers in family [cite: 745]
    "JWRIP",       # Vehicle occupancy [cite: 1504]
    "DRIVESP",     # Number of vehicles calculated from JWRI [cite: 2312]
    "DECADE",      # Decade of entry [cite: 2255]
    
    # Disability
    "DRAT"         # Veteran service connected disability rating [cite: 1401]
  ),
  
  # ========== BINARY (YES/NO) VARIABLES ==========
  binary = c(
    # Housing Amenities
    "BATH",        # Bathtub or shower [cite: 123]
    "BUS",         # Business or medical office on property [cite: 136]
    "KIT",         # Complete kitchen facilities [cite: 461]
    "PLM",         # Complete plumbing facilities [cite: 534]
    "REFR",        # Refrigerator [cite: 233]
    "RWAT",        # Hot and cold running water [cite: 265]
    "SINK",        # Sink with a faucet [cite: 268]
    "STOV",        # Stove or range [cite: 285]
    "TEL",         # Telephone in unit [cite: 287]
    "TOIL",        # Flush toilet [cite: 301]
    
    # Household/Family Flags
    "HUGCL",       # Grandchild living in housing unit [cite: 445]
    "LNGI",        # No one speaks English well [cite: 476]
    "MULTG",       # Multigenerational household [cite: 483]
    "NPP",         # Grandparent headed with no parent [cite: 488]
    "NR",          # Presence of nonrelative [cite: 489]
    "PSF",         # Presence of subfamilies [cite: 535]
    "PARTNER",     # Unmarried partner household [cite: 531]
    "RNTM",        # Meals included in rent [cite: 235]
    "SRNT",        # Specified rent unit [cite: 587]
    "SVAL",        # Specified value owner unit [cite: 588]
    
    # Health & Disability
    "DDRS",        # Self-care difficulty [cite: 1387]
    "DEAR",        # Hearing difficulty [cite: 1387]
    "DEYE",        # Vision difficulty [cite: 1387]
    "DOUT",        # Independent living difficulty [cite: 1389]
    "DPHY",        # Ambulatory difficulty [cite: 1390]
    "DREM",        # Cognitive difficulty [cite: 1415]
    "DIS",         # Disability recode [cite: 2318]
    "DRATX",       # Veteran disability rating checkbox [cite: 1404]
    "HICOV",       # Health insurance coverage recode [cite: 2656]
    "PRIVCOV",     # Private health insurance coverage recode [cite: 4265]
    "PUBCOV",      # Public health coverage recode [cite: 4301]
    "HINS1", "HINS2", "HINS3", "HINS4", "HINS5", "HINS6", "HINS7", # Ins types [cite: 1459-1483]
    
    # Demographics & Origins Recodes (Binary)
    "FER",         # Gave birth within past 12 months [cite: 1417]
    "LANX",        # Language other than English [cite: 1506]
    "GCL",         # Grandparents living with grandchildren [cite: 1418]
    "GCR",         # Grandparents responsible for grandchildren [cite: 1444]
    "RACAIAN",     # American Indian/Alaska Native recode [cite: 4349]
    "RACASN",      # Asian recode [cite: 4349]
    "RACBLK",      # Black recode [cite: 4349]
    "RACNHPI",     # Native Hawaiian/Pacific Islander recode [cite: 4349]
    "RACSOR",      # Some other race recode [cite: 4349]
    "RACWHT",      # White recode [cite: 4351]
    "RC",          # Related child [cite: 4352]
    "OC",          # Own child [cite: 3534]
    
    # Mortgage & Financial
    "FS",          # Food stamp/SNAP recipiency [cite: 167]
    "MRGI",        # Mortgage includes insurance [cite: 217]
    "MRGT",        # Mortgage includes taxes [cite: 227]
    "MRGX",        # First mortgage status [cite: 232]
    "SMX",         # Second/junior mortgage status [cite: 586]
    
    # Work & Education Flags
    "SCIENGP",     # Science/Engineering degree flag [cite: 4352]
    "SCIENGRLP",   # Science/Engineering related degree flag [cite: 4352]
    "WRK",         # Worked last week [cite: 1874]
    "NWAB",        # Temporary absence from work [cite: 1660]
    "NWAV",        # Available for work [cite: 1661]
    "NWLA",        # On layoff from work [cite: 1708]
    "NWLK",        # Looking for work [cite: 1718]
    "NWRE",        # Informed of recall [cite: 1719]
    
    # Marital & Military Events (Binary)
    "MARHD",       # Divorced in past 12 months [cite: 1548]
    "MARHM",       # Married in past 12 months [cite: 1555]
    "MARHW",       # Widowed in past 12 months [cite: 1570]
    "MLPA", "MLPB", "MLPC", "MLPD", "MLPE", "MLPF", 
    "MLPG", "MLPH", "MLPI", "MLPJ", "MLPK" # Military service periods [cite: 1600-1659]
  ),
  
  # ========== ALLOCATION FLAGS (Binary 0/1) ==========
  allocation_flags = c(
    # Housing Flags
    "FACRP", "FAGSP", "FBATHP", "FBDSP", "FBLDP", "FBUSP", 
    "FCONP", "FELEP", "FFSP", "FFULP", "FGASP", "FHFLP",
    "FINSP", "FKITP", "FMHP", "FMRGIP", "FMRGP", "FMRGTP",
    "FMRGXP", "FMVP", "FPLMP", "FREFRP", "FRMSP", "FRNTMP",
    "FRNTP", "FRWATP", "FSINKP", "FSMP", "FSMXHP", "FSMXSP",
    "FSTOVP", "FTAXP", "FTELP", "FTENP", "FTOILP", "FVACSP",
    "FVALP", "FVEHP", "FWATP", "FYBLP",
    
    # Person Flags
    "FAGEP", "FANCP", "FCITP", "FCITWP", "FCOWP", "FDDRSP",
    "FDEARP", "FDEYEP", "FDOUTP", "FDPHYP", "FDRATP", "FDRATXP",
    "FDREMP", "FENGP", "FESRP", "FFERP", "FFODP", "FGCLP",
    "FGCMP", "FGCRP", "FHINS1P", "FHINS2P", "FHINS3P", "FHINS4P",
    "FHINS5P", "FHINS6P", "FHINS7P", "FHISP", "FINDP", "FINTP",
    "FJWDP", "FJWMNP", "FJWRIP", "FJWTRP", "FLANP", "FLANXP",
    "FMARP", "FMARHDP", "FMARHMP", "FMARHTP", "FMARHWP", "FMARHYP",
    "FMIGP", "FMIGSP", "FMILPP", "FMILSP", "FOCCP", "FOIP",
    "FPAP", "FPOBP", "FPOWSP", "FRACP", "FRELP", "FRETP",
    "FSCHGP", "FSCHLP", "FSCHP", "FSEMP", "FSEXP", "FSSIP",
    "FSSP", "FWAGP", "FWKHP", "FWKLP", "FWKWP", "FWRKP",
    "FYOEP"
  ),
  
  # ========== WEIGHTS (Continuous but special use) ==========
  weights = c(
    paste0("WGTP", 1:80),    # Housing weight replicates
    paste0("PWGTP", 1:80),   # Person weight replicates
    "WGTP",                  # Main housing weight
    "PWGTP"                  # Main person weight
  )
)

# Create a dataframe to hold all covariates
all_covariates <- data.frame()

for (i in seq_along(full_dataset_covariates)) {
  puma_covariates <- full_dataset_covariates[[i]]
  all_covariates <- rbind(all_covariates, puma_covariates)
}

# Set as factors the categorical variables for the ones existing in the dataset
categorical_vars <- unlist(c(
  pums_variables$categorical_nominal,
  pums_variables$categorical_ordinal,
  pums_variables$binary,
  pums_variables$allocation_flags
))

for (var in categorical_vars) {
  if (var %in% colnames(all_covariates)) {
    all_covariates[[var]] <- as.factor(all_covariates[[var]])
  }
}

# Number of NA for each covariate
na_counts <- sapply(all_covariates, function(col) sum(is.na(col)))
#na_counts

# Remove columns with any NA values
all_covariates <- all_covariates[, na_counts == 0]
covariate_names <- colnames(all_covariates)

###########################################################################
# Gaussianity of LPINCP ====
###########################################################################

lpincp_values <- all_covariates$PINCP

# Box-Cox transformation to approximate normality
library(MASS)
boxcox_result <- boxcox(lpincp_values ~ 1, lambda = seq(-5, 5, 0.1))
best_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
cat("Best Box-Cox lambda for LPINCP:", best_lambda, "\n") # best_lambda = 0.2

# Apply Box-Cox transformation
if (best_lambda == 0) {
  lpincp_values <- log(lpincp_values)
} else {
  lpincp_values <- (lpincp_values^best_lambda - 1) / best_lambda
}

# QQ plot
qqnorm(lpincp_values, main = "QQ Plot of LPINCP")
qqline(lpincp_values, col = "red")

# Histogram density plot
hist(lpincp_values,
  breaks = 30,
  main = "Histogram of LPINCP",
  xlab = "LPINCP Values",
  col = "lightblue",
  probability = TRUE
)
# Overlay normal distribution curve
mu <- mean(lpincp_values, na.rm = TRUE)
sigma <- sd(lpincp_values, na.rm = TRUE)
x_seq <- seq(min(lpincp_values, na.rm = TRUE), max(lpincp_values, na.rm = TRUE), length.out = 100)
y_seq <- dnorm(x_seq, mean = mu, sd = sigma)
lines(x_seq, y_seq, col = "red", lwd = 2)

# Modify LPINCP in the dataset
all_covariates$LPINCP <- lpincp_values 

###########################################################################
# Linear Regression Analysis for Covariate Importance ====
###########################################################################

# Remove identifier and weight columns that shouldn't be predictors
exclude_cols <- c(
  "LPINCP", "PINCP", # Response variable and its non-log version
  "RT", "SERIALNO", "STATE_PUMA", # Identifiers
  "ADJINC",      # Adjustment factor
  pums_variables$allocation_flags,
  grep("^PWGTP", covariate_names, value = TRUE) # Person weights
)

# Get potential predictors
predictors <- setdiff(covariate_names, exclude_cols)

# Remove columns with zero or near-zero variance (constants or single-level factors)
valid_predictors <- sapply(predictors, function(col) {
  x <- all_covariates[[col]]
  if (is.factor(x) || is.character(x)) {
    # For factors/characters: need at least 2 levels
    return(length(unique(x)) >= 2)
  } else {
    # For numeric: need non-zero variance
    return(var(x, na.rm = TRUE) > 0)
  }
})
predictors <- predictors[valid_predictors]

cat("Number of predictors after filtering:", length(predictors), "\n")

# Fit a linear model with LPINCP as the response variable
# and all other covariates as predictors
formula_str <- paste("LPINCP ~", paste(predictors, collapse = " + "))
lm_model <- lm(as.formula(formula_str), data = all_covariates)

# Summary of the linear model
summary(lm_model)

# Extract coefficients, standard errors, t-values, and p-values
coef_summary <- summary(lm_model)$coefficients
coef_df <- as.data.frame(coef_summary)
colnames(coef_df) <- c("Estimate", "Std_Error", "t_value", "p_value")
coef_df$Covariate <- rownames(coef_df)
coef_df <- coef_df[coef_df$Covariate != "(Intercept)", ]

# Sort by absolute t-value (importance)
coef_df <- coef_df[order(abs(coef_df$t_value), decreasing = TRUE), ]

# Extract variable names from covariate names (removing level suffixes)
extract_base_var <- function(covariate_name) {
  # Strategy: Remove numeric+letter suffixes that represent levels
  # But preserve the base variable name even if it's all caps
  
  # First, handle specific patterns we know about:
  # Pattern 1: Variable ends with digits only (e.g., SCHL22 -> SCHL)
  # Pattern 2: Variable ends with digits + letters (e.g., NAICSP221MP -> NAICSP, SOCP1110XX -> SOCP)
  
  # Remove suffix that starts with a digit and ends with optional letters
  # This preserves base names like SEX, MAR, DIS, etc.
  base <- sub("[0-9]+[A-Z]*$", "", covariate_name)
  
  # If we removed everything (variable was all digits), return original
  if (base == "" || nchar(base) == 0) {
    return(covariate_name)
  }
  
  return(base)
}

# Aggregate coefficient results by base variable
coef_df$Base_Var <- sapply(coef_df$Covariate, extract_base_var)

# For each base variable, summarize:
var_importance <- coef_df %>%
  group_by(Base_Var) %>%
  summarise(
    Max_t_value = max(abs(t_value)),
    Min_p_value = min(p_value),
    n_levels = n(),
    Avg_t_value = mean(abs(t_value)),
    # Add: how many levels are "significant" (though with your sample size, nearly all will be)
    n_sig_levels = sum(p_value < 0.05),
    .groups = "drop"
  ) %>%
  arrange(desc(Max_t_value))

cat("\n=== Top 40 Most Important Variables (Grouped by Base Variable) ===\n")
print(head(var_importance, 40), n = 40)

# Add a note about interpretation
cat("\n=== INTERPRETATION NOTES ===\n")
cat("1. Max_t_value: Largest |t-value| among all levels of this variable\n")
cat("2. Avg_t_value: Average |t-value| across all levels (useful for multi-level variables)\n")
cat("3. n_levels: Number of dummy variables for this categorical variable\n")
cat("4. With large sample sizes, focus on t-values (effect sizes) rather than p-values\n")
cat("5. For multi-level variables (n_levels > 1), consider both Max and Avg t-values\n")
cat("6. Partial R² shows unique variance explained (most reliable importance metric)\n")