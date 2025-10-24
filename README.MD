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

👉 *[Algorithm comparison in docs](https://filippo-galli.github.io/tesi/)*

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

## 📁 Project Architecture

```
tesi/
├── 📂 R/                                      # R Analysis & Visualization Suite
│   ├── 🔬 sim_data_production.R                  # Data simulation (Natarajan/Gaussian/Gamma)
│   ├── 🚀 simulation_data.R                      # Main MCMC workflow orchestration  
│   ├── 📊 analysis.R                             # Post-processing and visualization
│   ├── 🛠️  utils.R                               # Utility functions and diagnostics
│   ├── � Likelihood_params_analisys.R           # Hyperparameter analysis tools
│   └── 🧪 U_cond_density.R                       # NGGP U-parameter density analysis
├── �📂 src/                                    # High-Performance C++ Core
│   ├── 📂 utils/                                  # Core Framework Components
│   │   ├── Sampler.hpp                              # Abstract MCMC sampler base class
│   │   ├── Process.hpp                              # Bayesian nonparametric process interface  
│   │   ├── Data.hpp/.cpp                            # Efficient cluster data management
│   │   ├── Likelihood.hpp/.cpp                      # Optimized likelihood computations
│   │   └── Params.hpp                               # Comprehensive parameter management
│   ├── � samplers/                               # MCMC Sampling Algorithms
│   │   ├── neal.hpp/.cpp                            # Neal's Algorithm 3 (collapsed Gibbs)
│   │   ├── neal_ZDNAM.hpp/.cpp                      # ZDNAM-enhanced Neal sampler
│   │   ├── splitmerge.hpp/.cpp                      # Split-Merge sampler  
│   │   ├── splitmerge_SAMS.hpp/.cpp                 # Sequential Allocation Merge-Split
│   │   └── � U_sampler/                            # Continuous Parameter Samplers
│   │       ├── RWMH.hpp                                 # Random Walk Metropolis-Hastings
│   │       └── MALA.hpp                                 # Metropolis-Adjusted Langevin
│   ├── 📂 processes/                              # Bayesian Nonparametric Processes
│   │   ├── DP.hpp/.cpp                              # Dirichlet Process
│   │   ├── DPW.hpp/.cpp                             # Weighted Dirichlet Process (spatial)
│   │   ├── NGGP.hpp/.cpp                            # Normalized Generalized Gamma Process
│   │   └── NGGPW.hpp/.cpp                           # Weighted NGGP (spatial)
│   └── 🔗 launcher.cpp                           # R-C++ integration interface
├── 📂 docs/                                   # Auto-Generated Documentation
│   ├── 📖 HTML files                              # Doxygen API reference
│   └── 🎨 Custom CSS                              # Documentation styling
├── 📂 doxygen-theme/                          # Doxygen Theme Components
├── 📂 simulation_data/                        # Simulated Datasets Repository
│   └── 📊 Natarajan_{σ}sigma_{d}d/              # Organized by parameters
│       ├── all_data.rds                             # Raw data points
│       ├── ground_truth.rds                         # True cluster labels
│       └── dist_matrix.rds                          # Distance matrix
├── 📂 results/                                # MCMC Analysis Outputs  
│   └── 📈 {process}_{sampler}_{init}_{params}/  # Results by configuration
│       ├── simulation_results.rds                   # MCMC output
│       ├── simulation_ground_truth.rds              # True labels
│       ├── simulation_data.rds                      # Original data
│       ├── simulation_distance_matrix.rds           # Distance matrix
│       ├── simulation_initial_params.rds            # Parameter object
│       └── plots/                                   # Generated visualizations
├── ⚙️  devenv.nix                              # Reproducible Development Environment
├── 🔧 devenv.yaml                             # Environment Configuration  
├── 🔒 devenv.lock                             # Locked Dependencies
├── 📋 Doxyfile                                # Doxygen configuration
└── 📝 README.MD                              # This comprehensive guide
```

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

### Step 1: Create Test Data

Generate some clustering data to practice with:

```bash
Rscript R/sim_data_production.R
```

**What this does:**
- Creates synthetic datasets with known cluster structures
- Saves them in the `simulation_data/` folder
- You can verify the clustering algorithms work correctly

**Customize your data** (open `R/sim_data_production.R` and edit):
```r
sigma <- 0.25  # How separated are clusters? (0.18=close, 0.25=far apart)
d <- 10        # How many dimensions? (10 or 50)
N <- 100       # How many data points?
```

The script supports different data types - see comments in the file for Gaussian or Gamma distributions.

### Step 2: Run the Clustering Analysis

Analyze your data to find clusters:

```bash
Rscript R/simulation_data.R
```

**What happens:**
1. Loads your data
2. Makes smart guesses about initial parameters (using k-means)
3. Runs the clustering algorithm
4. Shows you progress in real-time
5. Saves results for later analysis

**Choosing Your Algorithm** (edit `src/launcher.cpp`):

The framework lets you mix and match different components. Here's what's currently active:

```cpp
// Which clustering model? (uncomment one)
// DP process(data, param);              // Standard (simple, fast)
// DPW process(data, param);             // Standard + spatial awareness
// NGGP process(data, param);            // Flexible (handles complex shapes)
NGGPW process(data, param);              // Flexible + spatial ✓ CURRENTLY USED

// Which sampling method? (uncomment one)
// Neal3 sampler(...);                   // Fast (good for simple cases)
// SplitMerge sampler(...);              // Robust (better for tricky data)
SplitMerge_SAMS sampler(...);            // Optimized (best for large data) ✓ CURRENTLY USED
```

**Tuning the Analysis** (in `R/simulation_data.R`):

Most parameters are set automatically, but you can adjust these:

```r
param <- new(Params,
  hyperparams$delta1, hyperparams$alpha, hyperparams$beta,  # Auto-set ✓
  hyperparams$delta2, hyperparams$gamma, hyperparams$zeta,  # Auto-set ✓
  
  # How long to run? (adjust these)
  10000,        # Burn-in: initial warm-up phase (larger = more reliable)
  10000,        # Sampling: how many samples to collect
  
  # Clustering behavior (tune these if needed)
  2,            # How many clusters to expect? (higher = finds more clusters)
  0.1,          # Flexibility: how varied can cluster shapes be? (0.1 to 0.7)
  1,            # Advanced parameter (usually leave at 1)
  
  # Spatial awareness (for geographic data)
  1,            # How much to consider location? (0 = ignore, 1+ = consider)
  W             # Which points are neighbors?
)
```

**Quick tuning guide:**
- **Finding too many clusters?** Reduce the "expect" parameter (try 1 or 0.5)
- **Finding too few?** Increase it (try 3 or 4)
- **Need better quality?** Increase both burn-in and sampling (try 20000 each)
- **Fast test run?** Use 2000 burn-in, 5000 sampling

👉 *[Full parameter guide in documentation](https://filippo-galli.github.io/tesi/classParams.html)*

### Step 3: Visualize Results

Create plots and check how well it worked:

```bash
Rscript R/analysis.R
```

**Point it to your results** (edit the folder path in the script):
```r
folder <- "results/NGGPW_Neal_kmeans_0.2sigma_50d_BI2000_NI5000_a0.1_sigma0.7_tau1/"
```

**You'll get:**
- **Quality checks**: Did the algorithm converge? Is it stable?
- **Cluster visualizations**: Heatmaps and plots showing the discovered clusters
- **Performance scores**: How well do the clusters match reality? (for test data)
- **Beautiful plots**: Publication-ready figures saved in the `plots/` folder

**Understanding the output:**
- Stable horizontal lines in trace plots = good convergence ✓
- High ARI score (>0.8) = excellent clustering ✓
- Clear blocks in heatmap = distinct clusters ✓

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
├── simulation_ground_truth.rds         # True labels (for validation)
├── simulation_data.rds                 # Original data matrix
├── simulation_distance_matrix.rds      # Distance matrix used
├── simulation_initial_params.rds       # Params object (hyperparameters)
└── plots/                              # Generated visualizations
    ├── trace_plot_K.png
    ├── trace_plot_loglikelihood.png
    ├── posterior_similarity_matrix.png
    ├── cluster_assignments.png
    ├── distance_distribution.png
    └── ... (additional diagnostic plots)
```

**Naming Convention:**
- **Process**: `DP`, `DPW`, `NGGP`, `NGGPW`
- **Sampler**: `Neal`, `Neal_ZDNAM`, `SM` (SplitMerge), `SAMS`
- **Initialization**: `kmeans`, `sequential`, `single`
- **Parameters**: `{σ}sigma_{d}d_BI{burn-in}_NI{iterations}_a{mass}_sigma{σ}_tau{τ}`

## 🔧 For Developers

### Main Functions You'll Use

**In R (`R/utils.R`):**
```r
# Create test data
generate_mixture_data(N = 100, sigma = 0.25, d = 10)

# Set up parameters automatically
set_hyperparameters(data, dist_matrix, k_elbow = 3)

# Visualize results
plot_mcmc_results(mcmc_result, ground_truth, BI = 2000)

# Check if the right number of clusters is being found
plot_k_means(dist_matrix, max_k = 15)
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

### Working with the Example Data

The project includes pre-made examples:

```bash
# Generate example datasets
Rscript R/sim_data_production.R

# Run analysis on them
Rscript R/simulation_data.R

# Make plots
Rscript R/analysis.R
```

Each script is well-commented - open them to see exactly what's happening!

## ⚙️ Advanced Options

### Switching Algorithms

The default setup works well for most cases, but you can customize:

**In `src/launcher.cpp`, find these lines and uncomment what you want:**

```cpp
// Line ~110: Pick your clustering model
NGGPW process(data, params, likelihood);  // Current choice

// Line ~125: Pick your sampling method  
SplitMerge_SAMS sampler(...);  // Current choice
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


## 📚 References & Citations

- **Neal, R. M. (2000)**: "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
- **Jain, S. & Neal, R. M. (2004)**: "A Split-Merge Markov Chain Monte Carlo Procedure for the Dirichlet Process Mixture Model"  
- **Dahl, D. B. and Newcomb, S. (2022)**: "Sequentially allocated merge-split samplers for conjugate Bayesian nonparametric models"
- **Favaro, S. & Teh, Y. W. (2013)**: "MCMC for Normalized Random Measure Mixture Models"
- **Martinez, A. F. and Mena, R. H. (2014)**: "On a Nonparametric Change Point Detection Model in Markovian Regimes"
- **Natarajan, A. and  De Iorio, M. (2023)**: "Cohesion and Repulsion in Bayesian Distance Clustering"