# Bayesian Nonparametric Clustering with MCMC

A powerful framework for discovering clusters in your data without specifying the number of clusters in advance. This project combines the flexibility of Bayesian statistics with the efficiency of modern computational methods.

**📚 Full Documentation:** [API Reference & Technical Details](https://filippo-galli.github.io/tesi/)  
**🎯 What it does:** Automatically finds natural groupings in complex datasets  
**⚡ How it works:** R for easy analysis + C++ for fast computation
## ✨ What Makes This Special?

### 🎯 **Smart Clustering Methods**
You don't need to guess how many clusters exist in your data. The framework includes:
- **Basic Models (DP)**: Great for getting started with automatic clustering
- **Flexible Models (NGGP)**: For complex data with varying cluster shapes and sizes
- **Spatial Models (DPW, NGGPW)**: When your data points have locations or network connections
- **Auto-tuning**: The system helps estimate good starting parameters

� *[See technical details in documentation](https://filippo-galli.github.io/tesi/)*

### 🔄 **Multiple Analysis Methods**
Different algorithms for different needs:
- **Fast & Simple**: Neal's algorithm for quick results on well-separated clusters
- **Robust**: Split-Merge methods for tricky clustering problems
- **Optimized**: SAMS variant for large datasets

### 🌍 **Location-Aware Clustering**
Perfect for geographic data or networks:
- Automatically consider nearby points when forming clusters
- Works with any adjacency structure (neighbors, social networks, etc.)
- Adjustable spatial influence

### 📊 **Rich Visualization & Diagnostics**
Understand your results with:
- Real-time progress monitoring
- Professional plots and heatmaps
- Quality metrics (how well did it work?)
- Reproducible analysis scripts

## 🚀 Quick Start

### Option 1: Automatic Setup (Recommended)

The easiest way using [devenv](https://devenv.sh/) - it handles everything automatically:

1. **Install Nix** (one-time setup):
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf -L https://install.determinate.systems/nix | sh -s -- install
   ```

2. **Install devenv**:
   ```bash
   nix profile install --accept-flake-config github:cachix/devenv/latest
   ```

3. **Start working** (from the project folder):
   ```bash
   cd /path/to/tesi
   devenv shell
   ```

That's it! Everything you need is now installed and ready.

### Option 2: Manual Setup

If you prefer to install components yourself:

**In R:**
```r
install.packages(c(
  "Rcpp", "RcppEigen", "ggplot2", "dplyr", 
  "mcclust.ext", "salso", "MASS"
))
```

**System requirements:**
- R version 4.0 or newer
- A C++ compiler
- Eigen3 library (for linear algebra)

Most of these are pre-installed on modern systems or easily available through your package manager.

## 📖 How to Use

### Step 1: Prepare Your Data

You can work with either simulated or real data:

**Option A: Simulated Data**
Use the data retrieval scripts in `R/data_retrieval/`:
```r
source("R/data_retrieval/simulated_distances.R")
```

**Option B: Real Data**  
The framework includes support for real-world datasets. Example datasets are in the `real_data/` folder with pre-computed distance matrices and hyperparameters.

**Option C: Your Own Data**
Prepare your data as:
- An N×N distance/dissimilarity matrix
- Optional: spatial adjacency matrix W (for geographic data)
- Save as `.rds` files for easy loading

### Step 2: Run the Clustering Analysis

Analyze your data to find clusters:

```bash
Rscript R/launcher.R
```

**What happens:**
1. Loads your data (simulated or real)
2. Compiles the C++ MCMC implementation
3. Sets hyperparameters (automatically or from saved values)
4. Runs the clustering algorithm
5. Shows progress in real-time (logged to `mcmc_log.txt`)
6. Saves results in the `results/` folder

**Choosing Your Algorithm** (edit `src/launcher.cpp`):

The framework lets you mix and match different components:

```cpp
// Which clustering model? (uncomment one)
// DP process(data, param);              // Standard Dirichlet Process
// DPW process(data, param);             // DP + spatial awareness
// NGGP process(data, param);            // Normalized Generalized Gamma Process
NGGPW process(data, param);              // NGGP + spatial ✓ RECOMMENDED

// Which sampling method? (uncomment one)
// Neal3 sampler(...);                   // Neal's Algorithm 3 (fast, simple)
// Neal3_ZDNAM sampler(...);             // ZDNAM-enhanced Neal sampler
// SplitMerge sampler(...);              // Basic Split-Merge
// SplitMerge_SAMS sampler(...);            // Sequential Allocation MS 
// SplitMerge_LSS sampler(...);          // Launch State Split-Merge
SplitMerge_LSS_SDDS sampler(...);     // LSS with SDDS moves ✓ RECOMMENDED
```

**Tuning the Analysis** (in `R/launcher.R`):

Most parameters are set automatically, but you can adjust these:

```r
param <- new(Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,  # Auto-set ✓
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,  # Auto-set ✓
  
  # How long to run? (adjust these)
  10000,        # Burn-in: initial warm-up phase (larger = more reliable)
  10000,        # Sampling: how many samples to collect
  1,            # Prior mass parameter (higher = more clusters expected)
  
  # NGGP-specific parameters (for NGGP/NGGPW models)
  0.1,          # sigma: shape parameter (0.1 to 0.7)
  1,            # tau: scale parameter (usually leave at 1)
  
  # Spatial awareness (for DPW/NGGPW models)
  1,            # coefficient: spatial dependence strength (0 = ignore, 1+ = consider)
  dist_matrix,  # Distance/dissimilarity matrix
  W             # Spatial adjacency matrix
)
```

### Step 3: Visualize Results

Create plots and analyze performance:

```bash
Rscript R/mcmc_analysis.R
```

**Select your results** (the script will list available results):
```r
files <- list.files("results/")
file_chosen <- files[3]  # Choose which result to analyze
```

**You'll get:**
- **Trace plots**: Convergence diagnostics for number of clusters and log-likelihood
- **Posterior similarity matrix**: Heatmap showing clustering structure
- **Cluster assignments**: Distribution of cluster sizes and allocations
- **ACF plots**: Autocorrelation analysis for U parameter (NGGP models)
- **Geographic maps**: Spatial visualization of clusters (for real data with shapefiles)
- **Histograms**: Distribution analysis of cluster characteristics
- **All plots saved** in `results/{your_result}/VI_plots/` folder

## Output Files

### Generated Data (`simulation_data/`)
```
simulation_data/Natarajan_{σ}sigma_{d}d/
├── all_data.rds           # N×d matrix of data points
├── ground_truth.rds       # True cluster labels (length N vector)
└── dist_matrix.rds        # N×N distance matrix
```

### MCMC Results (`results/`)
Organized by configuration: `{process}_{sampler}_{init}_{params}/`

Example: `NGGPW_Neal_kmeans_0.2sigma_50d_BI2000_NI5000_a0.1_sigma0.7_tau1/`

```
results/{configuration}/
├── simulation_results.rds              # Main MCMC output
│   ├── $allocations                       # List of cluster assignments per iteration
│   ├── $K                                 # Number of clusters per iteration  
│   ├── $loglikelihood                     # Log-likelihood trace
│   └── $U                                 # NGGP U parameter trace (if applicable)
├── simulation_ground_truth.rds         # True labels (for simulated data)
├── simulation_data.rds                 # Original data matrix (for simulated data)
├── simulation_distance_matrix.rds      # Distance matrix used
├── simulation_initial_params.rds       # Params object (hyperparameters)
├── time_taken.txt                      # Elapsed time in seconds
└── VI_plots/                           # Generated visualizations
    ├── trace_K.png                        # Number of clusters over iterations
    ├── trace_loglik.png                   # Log-likelihood trace
    ├── post_sim_matrix.png                # Posterior similarity matrix heatmap
    ├── cluster_distribution.png           # Distribution of cluster sizes
    ├── cluster_estimates.png              # Point estimates of clustering
    ├── trace_U.png                        # U parameter trace (NGGP only)
    ├── acf_U.png                          # Autocorrelation of U (NGGP only)
    ├── map_clusters.png                   # Geographic visualization (real data)
    └── histogram_clusters.png             # Cluster characteristics
```

**Naming Convention:**
- **Data Type**: `simulation_data` or `real_data_{dataset}_{distance_metric}`
- **Process**: `DP`, `DPW`, `NGGP`, `NGGPW`
- **Sampler**: `Neal`, `Neal_ZDNAM`, `LSS`, `LSS_SDDS`, `SAMS`
- **Initialization**: `kmeans`, `sequential`, `single`
- **Parameters**: `BI{burn-in}_NI{iterations}_a{mass}_sigma{σ}_tau{τ}` (for NGGP variants)

## 🔧 For Developers

### Main Functions You'll Use

**In R:**

**Utilities (`R/utils.R`):**
```r
# Set up parameters automatically
set_hyperparameters(dist_matrix, k_elbow = 3)

# Retrieve spatial adjacency matrix
W <- retrieve_W(dist_matrix, k = 5)

# Check optimal number of clusters using elbow method
plot_k_means(dist_matrix, max_k = 15)

# Save results with informative naming
save_with_name(folder, param, filename)
```

**Plotting (`R/utils_plot.R`):**
```r
# Trace plots for diagnostics
plot_trace_cls(results, BI = 10000)
plot_trace_U(results, BI = 10000)

# Posterior analysis
plot_post_sim_matrix(results, BI = 10000)
plot_post_distr(results, BI = 10000)

# Cluster visualization
plot_cls_est(results, BI = 10000)
plot_map_cls(results, BI = 10000, unit_ids, puma_dir)
plot_hist_cls(results, BI = 10000)

# Autocorrelation analysis
plot_acf_U(results, BI = 10000)
```

**Main entry point (`src/launcher.cpp`):**
```cpp
// The function that does all the clustering work
Rcpp::List mcmc(distances, param, initial_allocations)
```

### Architecture Overview

The code is organized in layers:

```
R Scripts (Easy to use)
    ↓
launcher.cpp (Connects R to C++)
    ↓
Processes (DP, NGGP, etc.) ← What kind of clustering?
    ↓
Samplers (Neal, SplitMerge, etc.) ← How to find clusters?
    ↓
Data & Likelihood (Math happens here)
```
## 💡 Example Workflow

### Your First Analysis (Step by Step)

**1. Make some test data:**
```r
# Open R and run:
source("R/utils.R")
set.seed(123)

# Create 2 clear clusters
library(MASS)
cluster1 <- mvrnorm(30, mu = c(0, 0), Sigma = diag(2))
cluster2 <- mvrnorm(30, mu = c(5, 5), Sigma = diag(2))
my_data <- rbind(cluster1, cluster2)

# Save it
saveRDS(my_data, "my_test_data.rds")
```

**2. Run clustering:**
```r
# Load the framework
source("R/utils.R")
sourceCpp("src/launcher.cpp")

# Load your data
my_data <- readRDS("my_test_data.rds")
dist_matrix <- as.matrix(dist(my_data))

# Let it figure out good parameters
W <- retrieve_W(dist_matrix, k = 5)
hyperparams <- set_hyperparameters(my_data, dist_matrix, k_elbow = 2)

# Set up the analysis (using defaults)
param <- new(Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,
  2000, 5000,    # Burn-in: 2000, Sampling: 5000
  1, 0.5, 1,     # Expect ~2 clusters, moderate flexibility
  0, W           # No spatial (our data isn't geographic)
)

# Run it!
result <- mcmc(dist_matrix, param, hyperparams$initial_clusters - 1)
```

**3. See what you found:**
```r
# Make plots
plot_mcmc_results(result, BI = 2000)

# How many clusters did it find?
final_K <- result$K[2001:length(result$K)]  # After burn-in
mean(final_K)  # Should be close to 2
```

### Working with Real Data

The project supports analysis of real-world datasets:

```bash
# Run analysis on real data (e.g., municipality data)
Rscript R/launcher.R

# Analyze results with geographic visualization
Rscript R/mcmc_analysis.R
```

The framework automatically handles:
- Distance matrix loading
- Spatial adjacency matrices
- Pre-computed hyperparameters
- Geographic shapefiles for mapping

Each script is well-commented - open them to see exactly what's happening!

## ⚙️ Advanced Options

### Switching Algorithms

The default setup works well for most cases, but you can customize:

**In `src/launcher.cpp`, find these lines and uncomment what you want:**

```cpp
// Line ~110: Pick your clustering model
NGGPW process(data, params, likelihood);  // Current choice

// Line ~125: Pick your sampling method  
SplitMerge_LSS_SDDS sampler(...);  // Current choice
```

After changing, recompile in R:
```r
sourceCpp("src/launcher.cpp")
```

### Spatial/Geographic Data

If your data has locations or network structure:

```r
# Create neighbor connections (5 nearest neighbors)
W <- retrieve_W(dist_matrix, k = 5)

# Or use your own adjacency matrix
# W <- my_neighbor_matrix  # 1 if neighbors, 0 otherwise

# Tell the model to use spatial info
param <- new(Params, ...,
  coefficient = 1,  # How much to consider location (0 to 2)
  W                 # Your neighbor matrix
)
```

## Modularity & Extensibility

### Sampler or Process Addition
The codebase is designed to be modular, allowing you to easily add new processes, samplers, or likelihood functions. Each component is encapsulated in its own class, making it straightforward to extend functionality without modifying existing code. 
To add a new process or sampler, simply create a new class inheriting from the base class `src/utils/Process.hpp` or `src/utils/Sampler.hpp` and implement the required methods. This design promotes code reuse and simplifies maintenance.

### Likelihood addition 
At the moment there isn't a base class for likelihoods since at the moment I used the one specificied in the paper "Cohesion and Repulsion in Bayesian Distance Clustering" by Natarajan and De Iorio (2023). 
To add a new likelihood function, you can create a new class that implements the necessary methods for computing likelihoods based on your specific requirements following the same method that I exposed in the existing likelihood class. Following this approach will ensure consistency and integration with the existing framework without any problem.


## 📚 References & Citations

- **Neal, R. M. (2000)**: "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
- **Jain, S. & Neal, R. M. (2004)**: "A Split-Merge Markov Chain Monte Carlo Procedure for the Dirichlet Process Mixture Model"  
- **Dahl, D. B. and Newcomb, S. (2022)**: "Sequentially allocated merge-split samplers for conjugate Bayesian nonparametric models"
- **Favaro, S. & Teh, Y. W. (2013)**: "MCMC for Normalized Random Measure Mixture Models"
- **Martinez, A. F. and Mena, R. H. (2014)**: "On a Nonparametric Change Point Detection Model in Markovian Regimes"
- **Natarajan, A. and  De Iorio, M. (2023)**: "Cohesion and Repulsion in Bayesian Distance Clustering"