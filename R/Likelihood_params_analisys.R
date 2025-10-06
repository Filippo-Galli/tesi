library(ggplot2)

# Likelihood components for 2 points definition ----

cohesion_term <- function(distance, delta_1, alpha, beta) {
  # Cohesion term calculation
  term1 <- distance^(delta_1 - 1) / gamma(delta_1)
  term2 <- (beta)^(alpha) / gamma(alpha)
  term3 <- gamma(delta_1 + alpha) / (distance + beta)^(delta_1 + alpha)

  term1 * term2 * term3
}

repulsion_term <- function(distance, delta_2, zeta, gamma) {
  # Repulsion term calculation
  term1 <- distance^(delta_2 - 1) / gamma(delta_2)
  term2 <- (gamma)^(zeta) / gamma(zeta)
  term3 <- gamma(delta_2 + zeta) / (distance + gamma)^(delta_2 + zeta)

  term1 * term2 * term3
}

distances <- seq(0.01, 10, by = 0.1)

# Cohesion Analysis ----
delta_1_values <- seq(0.1, 0.9, by = 0.2)
alpha_values <- seq(1, 10, by = 1)
beta_values <- seq(1, 10, by = 1)

## Plotting Cohesion Term for different delta_1 values for all distances

cohesion_data <- expand.grid(distance = distances, delta_1 = delta_1_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    cohesion_data$delta_1, 
                    alpha = 1, beta = 1)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(delta_1))) +
  geom_line() +
  labs(x = "Distance", 
     y = "Cohesion Term",
     title = "Cohesion Term vs Distance - alpha = 1, beta = 1",
     color = "delta_1") +
  theme_minimal()

## Plotting Cohesion Term for different alpha values for all distances

cohesion_data <- expand.grid(distance = distances, alpha = alpha_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    0.9, 
                    alpha = cohesion_data$alpha, beta = 1)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(alpha))) +
  geom_line() +
  labs(x = "Distance",
     y = "Cohesion Term",
     title = "Cohesion Term vs Distance - delta_1 = 0.9, beta = 1",
     color = "alpha") +
  theme_minimal()

## Plotting Cohesion Term for different alpha values for all distances

cohesion_data <- expand.grid(distance = distances, beta = beta_values)
cohesion_data$cohesion_value <- mapply(cohesion_term, 
                    cohesion_data$distance, 
                    0.9, 
                    alpha = 1, beta = cohesion_data$beta)

ggplot(cohesion_data, aes(x = distance, y = cohesion_value, color = factor(beta))) +
  geom_line() +
  labs(x = "Distance",
       y = "Cohesion Term",
       title = "Cohesion Term vs Distance - delta_1 = 0.9, alpha = 1",
       color = "beta") +
  theme_minimal()

# Repulsion Analysis ----
delta_2_values <- seq(5, 100, by = 5)
zeta_values <- seq(1, 10, by = 1)
gamma_values <- seq(1, 10, by = 1)

## Plotting Repulsion Term for different delta_2 values for all distances
repulsion_data <- expand.grid(distance = distances, delta_2 = delta_2_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      repulsion_data$delta_2,
                      zeta = 1, gamma = 1)

ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(delta_2))) +
  geom_line() +
  labs(x = "Distance",
       y = "Repulsion Term",
       title = "Repulsion Term vs Distance - zeta = 1, gamma = 1",
       color = "delta_2") +
  theme_minimal()

## Plotting Repulsion Term for different zeta values for all distances
repulsion_data <- expand.grid(distance = distances, zeta = zeta_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      10,
                      zeta = repulsion_data$zeta, gamma = 1)  
ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(zeta))) +
  geom_line() +
  labs(x = "Distance",
        y = "Repulsion Term", 
        title = "Repulsion Term vs Distance - delta_2 = 10, gamma = 1",
        color = "zeta") +
  theme_minimal()

## Plotting Repulsion Term for different gamma values for all distances
repulsion_data <- expand.grid(distance = distances, gamma = gamma_values)
repulsion_data$repulsion_value <- mapply(repulsion_term,
                      repulsion_data$distance,
                      10,
                      zeta = 1, gamma = repulsion_data$gamma)
ggplot(repulsion_data, aes(x = distance, y = repulsion_value, color = factor(gamma))) +
  geom_line() +
  labs(x = "Distance",
       y = "Repulsion Term",
       title = "Repulsion Term vs Distance - delta_2 = 10, zeta = 1",
       color = "gamma") +
  theme_minimal()


# Plot Weidbull distribution for all the distances ----
weibull_data <- data.frame(distance = distances)
weibull_data$weibull_value <- dweibull(weibull_data$distance, shape = 3, scale = 2)
ggplot(weibull_data, aes(x = distance, y = weibull_value)) +
  geom_line(color = "blue") +
  labs(x = "Distance",
       y = "Weibull Density",
       title = "Weibull Distribution (shape=3, scale=2)") +
  theme_minimal()
