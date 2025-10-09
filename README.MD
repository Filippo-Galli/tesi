# Bayesian Nonparametric Clustering with MCMC

This project implements advanced Bayesian nonparametric clustering methods using state-of-the-art MCMC techniques. It provides a comprehensive framework for clustering analysis with Dirichlet Process (DP), Normalized Generalized Gamma Process (NGGP), and their weighted spatial variants, combining R for data analysis with high-performance C++ implementations.

**📚 Documentation:** [Doxygen API Reference](file:///home/filippo/Documents/tesi/docs/index.html)  
**🔬 Research Focus:** Bayesian nonparametric models with spatial dependencies  
**⚡ Performance:** Optimized C++ backend with R interface
## ✨ Key Features

### 🎯 **Advanced Clustering Models**
- **Dirichlet Process (DP)**: Classic nonparametric clustering with automatic cluster discovery
- **Normalized Generalized Gamma Process (NGGP)**: Enhanced flexibility and control over cluster structures  
- **Weighted Variants (DPW, NGGPW)**: Spatial dependency integration for geographic/network data
- **Automatic Model Selection**: Data-driven hyperparameter estimation

### 🔄 **State-of-the-Art MCMC Samplers**
- **Neal's Algorithm 3**: Efficient collapsed Gibbs sampling for standard applications
- **Split-Merge**: Advanced joint updates for improved mixing and faster convergence
- **SAMS (Sequential Allocation)**: Optimized proposal generation for large-scale problems
- **Hybrid Approaches**: Flexible combination of sampling strategies

### 🌍 **Spatial Dependencies**
- **Adjacency Matrix Support**: Incorporate spatial/network structure into clustering
- **Geographic Clustering**: Specialized algorithms for spatial data analysis

### 📊 **Comprehensive Analysis Suite**
- **Real-time Monitoring**: Live convergence diagnostics and progress tracking
- **Advanced Visualization**: Heatmaps, trace plots, and clustering evaluation
- **Performance Metrics**: ARI, silhouette scores, and posterior analysis
- **Reproducible Workflows**: Automated analysis pipelines with devenv integration

## 📁 Project Architecture

```
tesi/
├── 📂 R/                              # R Analysis & Visualization Suite
│   ├── 🔬 sim_data_production.R          # Advanced data simulation (Natarajan method)
│   ├── 🚀 simulation_data.R              # Main MCMC workflow orchestration  
│   ├── 📊 analysis.R                     # Post-processing and visualization
│   └── 🛠️  utils.R                       # Utility functions and diagnostics
├── 📂 src/                            # High-Performance C++ Core
│   ├── 🎯 Core Framework
│   │   ├── Sampler.hpp                    # Abstract MCMC sampler base class
│   │   ├── Process.hpp                    # Bayesian nonparametric process interface  
│   │   ├── Data.hpp                       # Efficient cluster data management
│   │   ├── Likelihood.hpp                 # Optimized likelihood computations
│   │   └── Params.hpp                     # Comprehensive parameter management
│   ├── 🔄 MCMC Samplers  
│   │   ├── neal.hpp/.cpp                  # Neal's Algorithm 3 (Gibbs sampling)
│   │   ├── splitmerge.hpp/.cpp            # Split-Merge sampler  
│   │   └── splitmerge_SAMS.hpp/.cpp       # Sequential Allocation Merge-Split
│   ├── 📈 Process Implementations
│   │   ├── DP.hpp/.cpp                    # Dirichlet Process
│   │   ├── NGGP.hpp/.cpp                  # Normalized Generalized Gamma Process
│   │   ├── DPW.hpp/.cpp                   # Weighted Dirichlet Process (spatial)
│   │   └── NGGPW.hpp/.cpp                 # Weighted NGGP (spatial)
│   └── 🔗 launcher.cpp                   # R-C++ integration interface
├── 📂 docs/                           # Auto-Generated Documentation
│   ├── 📖 Doxygen HTML documentation      # Comprehensive API reference
│   └── 🎨 doxygen-theme/                 # Custom documentation styling
├── 📂 simulation_data/                # Simulated Datasets Repository
│   └── 📊 Natarajan_*sigma_*d/           # Organized by parameters (σ, dimensions)
├── 📂 results/                        # MCMC Analysis Outputs  
│   └── 📈 {algorithm}_{config}/          # Results organized by method and parameters
├── ⚙️  devenv.nix                      # Reproducible Development Environment
├── 🔧 devenv.yaml                     # Environment Configuration  
├── 🔒 devenv.lock                     # Locked Dependencies
└── 📝 README.MD                      # This comprehensive guide
```

## Installation and Setup

### Prerequisites

This project uses [devenv](https://devenv.sh/) for reproducible development environments with Nix.

1. **Install Nix** (if not already installed):
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf -L https://install.determinate.systems/nix | sh -s -- install
   ```

2. **Install devenv**:
   ```bash
   nix profile install --accept-flake-config github:cachix/devenv/latest
   ```

3. **Enter the development environment**:
   ```bash
   cd /path/to/tesi
   devenv shell
   ```

This will automatically set up:
- R with all required packages (Rcpp, RcppEigen, ggplot2, etc.)
- C++ compiler toolchain
- MCMC analysis libraries (salso, mcclust, etc.)
- VS Code integration (language server, plot viewer)

### Manual Installation (Alternative)

If not using devenv, install the following:

**R Packages:**
```r
install.packages(c(
  "Rcpp", "RcppEigen", "ggplot2", "dplyr", "tidyr",
  "spam", "fields", "viridisLite", "RColorBrewer", 
  "pheatmap", "mcclust.ext", "mvtnorm", "gtools", 
  "salso", "MASS"
))
```

**System Dependencies:**
- R (≥ 4.0)
- C++ compiler with C++23 support
- Eigen3 library
- BLAS/LAPACK libraries

## Usage Guide

### 1. Data Simulation

First, generate simulated clustering data using the Natarajan method:

```bash
cd /path/to/tesi
Rscript R/sim_data_production.R
```

This script:
- Generates mixture data with configurable parameters (σ, dimensions)
- Creates distance matrices for spatial analysis
- Saves data to `simulation_data/` directory
- Supports different data generation methods (Gaussian, Gamma, Natarajan)

**Key Parameters:**
- `sigma`: Controls cluster separation (0.18, 0.2, 0.25)
- `d`: Dimensionality (10, 50)
- `N`: Number of data points (default: 100)

### 2. MCMC Analysis

Run the main clustering analysis:

```bash
Rscript R/simulation_data.R
```

This script performs:
- **Data Loading**: Reads simulated data and distance matrices
- **Hyperparameter Estimation**: Uses k-means and elbow method for initialization
- **MCMC Sampling**: Runs Bayesian clustering with selected algorithm
- **Real-time Monitoring**: Displays progress and convergence diagnostics

**Algorithm Configuration (in `src/launcher.cpp`):**

The framework provides multiple combinations of **Process × Sampler × Spatial** configurations:

```cpp
// 🎯 PROCESS SELECTION
// DP              - Dirichlet Process (classic nonparametric prior)
// NGGP            - Normalized Generalized Gamma Process (enhanced flexibility)
// DPW/NGGPW       - Weighted variants with spatial dependencies

// 🔄 SAMPLER SELECTION  
// Neal3           - Collapsed Gibbs sampling (efficient, standard)
// SplitMerge      - Joint updates (better mixing, complex structures)
// SplitMerge_SAMS - Sequential allocation (optimized proposals)

// Example configurations:
DP process(data, param); // Basic DP process
Neal3 neal_sampler(data, params, likelihood, process); 

```

**MCMC Parameters:**
- `BI`: Burn-in iterations (default: 2000)
- `NI`: Sampling iterations (default: 10000)
- `a`: Total mass parameter (default: 1)
- `sigma`: NGGP parameter (default: 0.5)
- `tau`: NGGP parameter (default: 1.0)

### 3. Results Analysis

Analyze and visualize results:

```bash
Rscript R/analysis.R
```

This generates:
- **Convergence Diagnostics**: Trace plots, autocorrelation analysis
- **Clustering Results**: Posterior similarity matrices, cluster assignments
- **Performance Metrics**: Adjusted Rand Index (ARI), cluster distribution
- **Visualization**: Heatmaps, scatter plots, distance distributions

## Output Files

### Generated Data
- `simulation_data/Natarajan_{sigma}sigma_{d}d/`
  - `all_data.rds`: Raw data points
  - `ground_truth.rds`: True cluster labels
  - `dist_matrix.rds`: Distance matrix

### MCMC Results
- `results/{algorithm}_{initialization}_{parameters}/`
  - `simulation_results.rds`: MCMC output (allocations, K, log-likelihood)
  - `simulation_ground_truth.rds`: True labels
  - `simulation_data.rds`: Original data
  - `simulation_distance_matrix.rds`: Distance matrix
  - `plots/`: Generated visualizations

## 🔧 Core Components & API

### 📊 **Data Generation & Simulation** (`utils.R`)
```r
generate_mixture_data(N, sigma, d)     # Advanced Natarajan mixture simulation
distance_plot(data, labels)            # Intra/inter-cluster distance analysis  
```

### ⚙️ **Hyperparameter Optimization** (`utils.R`)  
```r
set_hyperparameters(data, dist_matrix) # Automatic parameter estimation via k-means
plot_k_means(data, max_k)              # Elbow method for optimal cluster number
```

### 🚀 **High-Performance MCMC Interface** (`launcher.cpp`)
```cpp
// Core classes (fully documented with Doxygen)
class Params;           // Parameter management with validation
class Data;             // Cluster data structures  
class Likelihood;       // Loglikelihood computations
class Process;          // Abstract nonparametric process interface
class Sampler;          // Abstract MCMC sampler base class

// Main interface
Rcpp::List mcmc(matrix, params, initial_clusters); // Primary MCMC function
```

### 📚 **Class Hierarchy** (See [Doxygen Docs](docs/index.html))
```
Sampler (abstract base)
├── Neal3              # Algorithm 3 implementation
├── SplitMerge         # Standard split-merge sampler  
└── SplitMerge_SAMS    # Sequential allocation variant

Process (abstract base)  
├── DP                 # Dirichlet Process
├── DPW                # Weighted Dirichlet Process
├── NGGP               # Normalized Generalized Gamma Process
└── NGGPW              # Weighted NGGP
```

## Example Workflow

1. **Generate Data**:
   ```r
   # In R/sim_data_production.R
   sigma <- 0.25
   d <- 10
   data_generation <- generate_mixture_data(N = 100, sigma = sigma, d = d)
   ```

2. **Run MCMC**:
   ```r
   # In R/simulation_data.R
   sourceCpp("src/launcher.cpp")
   hyperparams <- set_hyperparameters(all_data, dist_matrix, k_elbow = 3)
   param <- new(Params, hyperparams$delta1, hyperparams$alpha, ...)
   mcmc_result <- mcmc(dist_matrix, param, hyperparams$initial_clusters)
   ```

3. **Analyze Results**:
   ```r
   # Automatic analysis with burn-in
   plot_mcmc_results(mcmc_result, ground_truth, BI = 2000)
   ```

## ⚙️ Advanced Configuration

### 🔧 **Compilation & Build System**
```r
# Automatic compilation with dependency management
sourceCpp("src/launcher.cpp")  # Compiles entire C++ framework
```

### 🌍 **Spatial Dependency Configuration**
```cpp
// In Params constructor
Params params(
    delta1, alpha, beta,           // Process hyperparameters
    delta2, gamma, zeta,           // Likelihood hyperparameters  
    BI, NI,                        // MCMC iterations
    a, sigma, tau,                 // Process parameters
    coefficient,                   // Spatial weight (0 = no spatial, 1+ = strong spatial)
    W                             // Adjacency matrix (k-NN or custom)
);
```

### 📊 **Hyperparameter Tuning Guide**
```r
# Process parameters (control clustering behavior)
a      <- 1.0      # Total mass (higher = more clusters)
sigma  <- 0.5      # NGGP flexibility (0-1, higher = more flexible)  
tau    <- 1.0      # NGGP tail parameter

# MCMC parameters (control sampling quality)
BI     <- 2000     # Burn-in (increase if slow convergence)
NI     <- 10000    # Iterations (increase for precision)

# Spatial parameters  
coefficient <- 0.5  # Spatial strength (tune based on spatial structure)
k_neighbors <- 5    # Adjacency matrix connectivity
```

## 🛠️ Troubleshooting & FAQ

### 🔧 **Installation Issues**
```bash
# Missing dependencies  
sudo apt-get install r-base-dev libeigen3-dev liblapack-dev libblas-dev

# Rcpp compilation errors
R -e "install.packages('Rcpp', type='source')"

# devenv not working
nix --version  # Ensure Nix is installed
devenv --version  # Verify devenv installation
```

### 📞 **Getting Help**
- 📚 **Documentation**: [Doxygen API Reference](docs/index.html)
- 🔍 **Code Examples**: See `R/simulation_data.R` for complete workflows  
- 🧪 **Test Cases**: Run `R/sim_data_production.R` for validation
- 📖 **Theory**: References in code comments and Doxygen documentation

## 📚 References & Citations

- **Neal, R. M. (2000)**: "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
- **Jain, S. & Neal, R. M. (2004)**: "A Split-Merge Markov Chain Monte Carlo Procedure for the Dirichlet Process Mixture Model"  
- **Dahl, D. B. and Newcomb, S. (2022)**: "Sequentially allocated merge-split samplers for conjugate Bayesian nonparametric models"
- **Favaro, S. & Teh, Y. W. (2013)**: "MCMC for Normalized Random Measure Mixture Models"
- **Martinez, A. F. and Mena, R. H. (2014)**: "On a Nonparametric Change Point Detection Model in Markovian Regimes"
- **Natarajan, A. and  De Iorio, M. (2023)**: "Cohesion and Repulsion in Bayesian Distance Clustering"