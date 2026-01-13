# BNPClust: A C++ Framework for Bayesian Nonparametric Clustering with MCMC

Ispired by the need for flexible and efficient tools in Bayesian nonparametric clustering, this repository presents a C++ framework that leverages Markov Chain Monte Carlo (MCMC) methods. The framework is designed to be modular and extensible, allowing users to implement various stochastic processes, likelihood models, and sampling strategies.

[Link to the docs](https://filippo-galli.github.io/tesi/)

## Directory Structure

- `R`: R interface functions and scripts.
- `docs`: Documentation.
- `doxygen-theme`: Doxygen theme files for documentation generation.
- `src`: C++ source code implementing the core framework.

### R 
Each class of the framework has an R bindings to facilitate interaction with R users. Personally, I use `R/launcher.R` to load the data, set parameters, call the MCMC functions and save the results. The MCMC loop is inside `R/mcmc_loop.R`, here you can choose which combination of Process, Likelihood and Sampler to use in total freedom.
To analyze a specific chain I use `R/mcmc_analysis.R`, which contains functions to visualize and summarize the MCMC output.
Other folder or script are utility functions to load data, plot results, fetch and clean data, etc.

### docs
This folder contains the documentation files generated with Doxygen. You can find here the class diagrams and the detailed documentation of the code. To regenerate the documentation, you need to have Doxygen installed and run the command `doxygen Doxyfile` in the root directory.

### doxygen-theme
This folder contains the theme files for Doxygen to generate the documentation.

### src
This folder contains the C++ source code implementing the core framework. The main components are organized in subfolders:
- `processes`: Contains the implementations of the stochastic processes (e.g., Dirichlet Process, NGGP) and their modules (e.g. continuos covariate modules, spatial modules and so on).
- `likelihoods`: Contains the implementations of the likelihood models (e.g., the in Natarajan et al., 2023).
- `samplers`: Contains the implementations of the MCMC samplers (e.g., Neal's Algorithm 3, Split-Merge).
- `utils`: Contains utility functions and classes used throughout the codebase and all the base classes.

The file `src/bindings.cpp` contains the Rcpp bindings to expose the C++ classes and functions to R. This make possible to create your own MCMC scheme directly from R by combining the different components.

## 🏗️ Architecture

The framework is built around four main logical components:

1.  **`Params`**: Manages model hyperparameters and MCMC settings.
2.  **`Likelihood`**: Defines the data observation model.
3.  **`Process`**: Handles the Bayesian nonparametric prior.
4.  **`Data`**: Manages cluster assignments and data handling.
5.  **`Sampler`**: Implements the MCMC inference engine.

## 📦 Implemented Methods

### Stochastic Processes (`Process`)
Already implemented processes include:
- **Dirichlet Process (DP)**: The classic Dirichlet Process.
- **Normalized Generalized Gamma Process (NGGP)**: A more general family of which the DP is a special case.

There is also a modular structure to add covariate-dependent or spatial modules to the processes, such as:
- **Continuous Covariate Module**: Incorporates continuous covariates into the clustering process.
- **Spatial Module**: Accounts for spatial dependencies in the clustering given a matrix of neighbors.
- **Binary Covariate Module**: Incorporates binary covariates into the clustering process.

Along with these there are their respective cached version to speed up computations if needed.

### Likelihood Models (`Likelihood`)
Every likelihood components is implemented as log-likelihood function, allowing more robust computation. Currently implemented likelihoods are:
- **Distance-based Clustering**: Based on "Cohesion and Repulsion in Bayesian Distance Clustering" (Natarajan et al., 2023).
- **Gamma Likelihood**: Which is based on the previous one without the repulsion term.
- **Null Likelihood**: A placeholder likelihood that does not contribute to the posterior, useful for prior inspection.

### MCMC Samplers (`Sampler`)
The framework includes several MCMC sampling strategies:
- **Neal's Algorithm 3**: Conjugate Gibbs sampler (Neal, 2000).
- **ZDNAM**: Gibbs sampling with Zero-self Downward Nested Antithetic Modification (Neal, 2024) (WIP).
- **Split-Merge (Jain & Neal)**: Standard Split-Merge MCMC (Jain & Neal, 2004).
- **SAMS**: Sequentially-Allocated Merge-Split (Dahl, 2021).
- **LSS Split-Merge**: Locality Sensitive Sampling for scalable Split-Merge (Luo et al., 2018).
- **LSS-SDDS Split-Merge**: Split-Merge with Smart-Dumb/Dumb-Smart moves, LSS and SAMS.

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
Everything is binded to R using Rcpp in the file `src/bindings.cpp`. You can create your own MCMC scheme directly from R by combining the different components. Use the provided `R/launcher.R` as a starting point to set up and run your MCMC scheme as in `R/mcmc_loop.R`.

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

