# MCMC Clustering with Dirichlet Process and Normalized Generalized Gamma Process

This project implements Bayesian nonparametric clustering methods using Markov Chain Monte Carlo (MCMC) techniques. It combines R for data analysis and visualization with C++ for high-performance MCMC sampling, focusing on Dirichlet Process (DP) and Normalized Generalized Gamma Process (NGGP) mixture models with spatial dependency.

## Features

- **Bayesian Clustering**: Implementation of DP and NGGP mixture models
- **Spatial Dependencies**: Support for spatial clustering with adjacency matrices
- **Multiple Samplers**: Neal's Algorithm 2 and Split-Merge sampling strategies
- **Data Generation**: Simulated data generation from mixture models (Natarajan method)
- **Comprehensive Analysis**: Visualization, convergence diagnostics, and clustering evaluation
- **C++/R Integration**: High-performance C++ backend with R interface using Rcpp

## Project Structure

```
├── R/                          # R scripts for analysis and utilities
│   ├── sim_data_production.R       # Data generation and simulation
│   ├── simulation_data.R           # Main MCMC analysis script
│   ├── analysis.R                  # Results analysis and visualization
│   └── utils.R                     # Utility functions and plotting
├── src/                        # C++ source code
│   ├── launcher.cpp                # Main MCMC interface
│   ├── DP_*.cpp                    # Dirichlet Process implementations
│   ├── NGGP_*.cpp                  # NGGP implementations
│   ├── Data.hpp                    # Cluster management
│   ├── Likelihood.hpp              # Likelihood computation
│   └── Params.hpp                  # Parameter management
├── simulation_data/            # Generated simulation datasets
├── results/                    # MCMC results and analysis outputs
├── devenv.nix                  # Development environment configuration
└── README.MD                   # This file
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
- C++ compiler with C++11 support
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

**Algorithm Selection (in `src/launcher.cpp`):**
```cpp
// Uncomment the desired sampler:
// DPNeal2 sampler(data, param, likelihood);           // DP + Neal's Algorithm 2
// DPNeal2W sampler(data, param, likelihood);          // DP + Neal 2 + Spatial
// DPSplitMerge sampler(data, param, likelihood);      // DP + Split-Merge
// DPSplitMergeW sampler(data, param, likelihood);     // DP + Split-Merge + Spatial
// NGGPNeal2 sampler(data, param, likelihood);         // NGGP + Neal 2
// NGGPNeal2W sampler(data, param, likelihood);        // NGGP + Neal 2 + Spatial
// NGGPSplitMerge sampler(data, param, likelihood);    // NGGP + Split-Merge
NGGPSplitMergeW sampler(data, param, likelihood);      // NGGP + Split-Merge + Spatial (default)
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

## Key Functions

### Data Generation (`utils.R`)
- `generate_mixture_data()`: Natarajan mixture model simulation
- `distance_plot()`: Intra vs inter-cluster distance visualization

### Hyperparameter Setting (`utils.R`)
- `set_hyperparameters()`: Automatic hyperparameter estimation using k-means
- `plot_k_means()`: Elbow method for optimal K selection

### Analysis (`utils.R`)
- `plot_mcmc_results()`: Comprehensive MCMC analysis and visualization
- `plot_data()`: Basic data and clustering visualization

### MCMC Interface (`launcher.cpp`)
- `mcmc()`: Main MCMC sampling function
- `Params`: Parameter management class

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

## Configuration

### Compilation
The project automatically compiles C++ code when sourcing:
```r
sourceCpp("src/launcher.cpp")
```

### Algorithm Selection
Modify `src/launcher.cpp` to choose different sampling strategies:
- **DP vs NGGP**: Different nonparametric priors
- **Neal2 vs SplitMerge**: Different sampling algorithms  
- **W vs non-W**: With or without spatial dependencies

### Spatial Dependencies
- Adjacency matrix `W` controls spatial clustering
- `coefficient` parameter weights spatial vs non-spatial terms
- `retrieve_W()` generates adjacency from distance (k-nearest neighbors)

## Performance

- **C++ Backend**: High-performance MCMC sampling
- **Progress Monitoring**: Real-time iteration tracking
- **Memory Efficient**: Streaming results storage
- **Parallel Potential**: Ready for multi-threading extensions

## Troubleshooting

1. **Compilation Issues**: Ensure Eigen3 and R development headers are installed
2. **Package Dependencies**: Use `devenv shell` for automatic setup
3. **Memory Issues**: Reduce `NI` parameter for large datasets
4. **Convergence**: Increase burn-in `BI` if trace plots show non-stationarity