# Bayesian Nonparametric Clustering with MCMC

A flexible and efficient C++ framework for Bayesian nonparametric clustering using Markov Chain Monte Carlo (MCMC) methods, with R interfaces.

This repository provides a modular codebase designed for researchers and practitioners. While the current implementation focuses on **distance-based clustering** (Natarajan et al., 2023), the architecture allows for easy extension to other models, likelihoods, and sampling strategies.

[Link to the docs](https://filippo-galli.github.io/tesi/)

## Directory Structure

- `R`: R interface functions and scripts.
- `docs`: Documentation.
- `doxygen-theme`: Doxygen theme files for documentation generation.
- `src`: C++ source code implementing the core framework.

### R 
This folder contains some helper script to load and run the C++ code from R. The main entry point is in `R/launcher.R`.
There are also some script to do mcmc post-processing and plotting.

### docs
This folder contains the documentation files generated with Doxygen. You can find here the class diagrams and the detailed documentation of the code. To generate the documentation, you need to have Doxygen installed and run the command `doxygen Doxyfile` in the root directory.

### doxygen-theme
This folder contains the theme files for Doxygen to generate the documentation.

### src
This folder contains the C++ source code implementing the core framework. The main components are organized in subfolders:
- `processes`: Contains the implementations of the stochastic processes (e.g., Dirichlet Process, NGGP).
- `likelihoods`: Contains the implementations of the likelihood models (e.g., distance-based clustering).
- `samplers`: Contains the implementations of the MCMC samplers (e.g., Neal's Algorithm 3, Split-Merge).
- `utils`: Contains utility functions and classes used throughout the codebase and all the base classes.
The entry point for the C++ code is in `src/launcher.cpp`, which is exposed to R via Rcpp.

## 🏗️ Architecture

The framework is built around four main logical components:

1.  **`Params`**: Manages model hyperparameters and MCMC settings.
2.  **`Likelihood`**: Defines the data observation model.
3.  **`Process`**: Handles the Bayesian nonparametric prior (e.g., Dirichlet Process) and cluster assignments.
4.  **`Sampler`**: Implements the MCMC inference engine.

### Extensibility
The code is structured to facilitate the integration of new components. You can implement custom models by extending the base classes (`Likelihood`, `Process`, `Sampler`).

> **Note**: When modifying `Params` or adding new classes, ensure you update the R bindings in `src/launcher.cpp`.

## 📦 Implemented Methods

### Stochastic Processes (`Process`)
- **Dirichlet Process (DP)**
- **Normalized Generalized Gamma Process (NGGP)**
- *Modules*: Support for spatial/geographical constraints via similarity matrices.

### Likelihood Models (`Likelihood`)
- **Distance-based Clustering**: Based on "Cohesion and Repulsion in Bayesian Distance Clustering" (Natarajan et al., 2023).

### MCMC Samplers (`Sampler`)
- **Neal's Algorithm 3**: Conjugate Gibbs sampler (Neal, 2000).
- **ZDNAM**: Gibbs sampling with Zero-self Downward Nested Antithetic Modification (Neal, 2024).
- **Split-Merge (Jain & Neal)**: Standard Split-Merge MCMC (Jain & Neal, 2004).
- **SAMS**: Sequentially-Allocated Merge-Split (Dahl, 2021).
- **LSS Split-Merge**: Locality Sensitive Sampling for scalable Split-Merge (Luo et al., 2018).
- **LSS-SDDS Split-Merge**: Smart-Dumb/Dumb-Smart moves with LSS.

## 🛠️ Installation & Usage

### Prerequisites
- **R**
- **C++ Compiler** (supporting C++11 or later)
- R Packages: `Rcpp`, `RcppEigen`

### Setup
1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/repo-name.git
    cd repo-name
    ```
2.  Open the project in R/RStudio.
3.  Ensure dependencies are installed:
    ```R
    install.packages(c("Rcpp", "RcppEigen"))
    ```

### Running the Code
The entry point for the C++ code is exposed to R. Currently, the specific combination of Process, Likelihood, and Sampler is configured in `src/launcher.cpp`.

1.  Modify `src/launcher.cpp` to select your desired model components.
2.  Compile and source the code in R (e.g., using `Rcpp::sourceCpp`).
3.  Run the provided R functions to start the MCMC chain.

## 📚 References

- **Natarajan et al. (2023)**: *Cohesion and Repulsion in Bayesian Distance Clustering*.
- **Neal (2000)**: *Markov Chain Sampling Methods for Dirichlet Process Mixture Models*.
- **Neal (2024)**: *Modifying Gibbs Sampling to Avoid Self Transitions*.
- **Jain & Neal (2004)**: *A Split-Merge Markov Chain Monte Carlo Procedure for the Dirichlet Process Mixture Model*.
- **Dahl (2021)**: *Sequentially-allocated merge-split samplers for conjugate Bayesian nonparametric models*.
- **Luo et al. (2018)**: *Scaling-up Split-Merge MCMC with Locality Sensitive Sampling (LSS)*.

## 🤝 Contributing

Contributions are welcome! Feel free to submit Pull Requests or open Issues for bugs and feature requests.

---
*If you find this project useful, please leave a star! ⭐*

