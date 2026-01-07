library(Rcpp)
library(RcppEigen)

# Load the C++ module
sourceCpp("src/bindings.cpp")

run_mcmc <- function(params, covariates, initial_allocations = integer(0)) {
    # Ensure types are correct for C++
    initial_allocations <- as.integer(initial_allocations)

    cache <- create_Covariate_cache(covariates, initial_allocations)

    # Instantiate Data using factory function
    data <- create_Data_wClusterInfo(params, cache, initial_allocations)
    # data <- create_Data(params, initial_allocations)

    # Instantiate Likelihood using factory function
    likelihood <- create_Natarajan_likelihood(data, params)

    # Instantiate U_sampler (RWMH) using factory function
    # Constructor: Params&, Data&, bool use_V, double proposal_sd, bool tuning_enabled
    u_sampler <- create_RWMH(params, data, TRUE, 2.0, TRUE)

    # Instantiate Process (NGGPx) using modules
    # 1. Spatial module
    mod_spatial <- create_SpatialModule(covariates, data)
    # 2. Covariate module (cached)
    mod_cov <- create_CovariatesModuleCache(covariates, data, cache)
    # Combine modules into NGGPx process
    process <- create_NGGPx(data, params, u_sampler, list(mod_spatial, mod_cov))

    # Instantiate Sampler (SplitMerge_LSS_SDDS) using factory function
    # Constructor: Data&, Params&, Likelihood&, Process&, bool shuffle
    sampler <- create_SplitMerge_LSS_SDDS(data, params, likelihood, process, TRUE)

    neal3 <- create_Neal3(data, params, likelihood, process)

    # Get parameters for loop using getter functions
    BI <- params_get_BI(params)
    NI <- params_get_NI(params)
    total_iters <- BI + NI

    # Results storage
    allocations_out <- vector("list", total_iters)
    K_out <- integer(total_iters)
    U_out <- numeric(total_iters)

    cat("Starting MCMC with", NI, "iterations after", BI, "burn-in...\n")

    start_time <- Sys.time()

    for (i in 1:total_iters) {
        # Update process parameters (U)
        process_update_params(process)

        # MCMC Step
        sampler_step(sampler)

        # Neal3 Step
        if (i %% 25 == 0) {
            sampler_step(neal3)
        }
        # Store results
        allocations_out[[i]] <- data_get_allocations(data)
        K_out[i] <- data_get_K(data)
        U_out[i] <- u_sampler_get_U(u_sampler)

        # Progress
        if (i %% max(1, floor(total_iters / 20)) == 0) {
            elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            iter_per_sec <- i / elapsed
            eta <- (total_iters - i) / iter_per_sec
            cat(sprintf("Iteration %d: Clusters: %d - iter/s: %.2f eta: %.2f\n ", i, data_get_K(data), iter_per_sec, eta))
        }
    }

    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat("MCMC completed.\n")
    cat("Total time (secs):", elapsed_time, "\n")
    cat("U acceptance rate:", u_sampler_get_acceptance_rate(u_sampler) * 100, "%\n")

    return(list(
        allocations = allocations_out,
        K = K_out,
        U = U_out,
        BI = BI,
        NI = NI,
        elapsed_time = elapsed_time
    ))
}
