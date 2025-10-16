library(shiny)

# Define the conditional density function
conditional_density_V <- function(v, K, n, a, sigma, tau) {
  
  # Pre-compute exp(v)
  exp_v <- exp(v)
  
  # Pre-compute frequently used values
  a_over_sigma <- a / sigma
  tau_power_sigma <- tau^sigma
  
  # Compute log density components
  # log(e^{vn}) = vn
  term1 <- v * n
  
  # log((e^v + τ)^{n-a|π|}) = -(n - a*K) * log(e^v + τ)
  term2 <- -(n - a * K) * log(exp_v + tau)
  
  # -(a/σ)((e^v+τ)^σ - τ^σ)
  term3 <- -a_over_sigma * ((exp_v + tau)^sigma - tau_power_sigma)
  
  # Return density
  exp(term1 + term2 + term3)
}

ui <- fluidPage(
  titlePanel("Conditional Density V(v | ·)"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("K", 
                  "K (number of categories):", 
                  min = 1, 
                  max = 40, 
                  value = 20, 
                  step = 1),
      
      sliderInput("n", 
                  "n (sample size):", 
                  min = 10, 
                  max = 200, 
                  value = 100, 
                  step = 10),
      
      sliderInput("a", 
                  "a (parameter):", 
                  min = 0.01, 
                  max = 5, 
                  value = 0.1, 
                  step = 0.1),
      
      sliderInput("sigma", 
                  "σ (sigma):", 
                  min = 0.1, 
                  max = 1, 
                  value = 0.7, 
                  step = 0.1),
      
      sliderInput("tau", 
                  "τ (tau):", 
                  min = 0.01, 
                  max = 5, 
                  value = 1, 
                  step = 0.01),
      
      hr(),
      
      sliderInput("v_range", 
                  "v range:", 
                  min = -5, 
                  max = 20, 
                  value = c(-5, 20), 
                  step = 0.5)
    ),
    
    mainPanel(
      plotOutput("densityPlot", height = "500px"),
      hr(),
      verbatimTextOutput("paramInfo")
    )
  )
)

server <- function(input, output) {
  
  output$densityPlot <- renderPlot({
    # Generate v values
    v <- seq(input$v_range[1], input$v_range[2], length.out = 500)
    
    # Compute density for each v
    density_values <- sapply(v, function(v_i) {
      conditional_density_V(v_i, 
                           K = input$K, 
                           n = input$n, 
                           a = input$a, 
                           sigma = input$sigma, 
                           tau = input$tau)
    })
    
    # Handle numerical issues (Inf, NA)
    valid_idx <- is.finite(density_values) & density_values > 0
    
    if (sum(valid_idx) > 0) {
      # Plot
      plot(v[valid_idx], density_values[valid_idx], 
           type = "l", 
           lwd = 2, 
           col = "steelblue",
           main = "Conditional Density of V",
           xlab = "v",
           ylab = "Density",
           las = 1)
      
      grid()
      
      # Add vertical line at mode (approximate)
      mode_idx <- which.max(density_values[valid_idx])
      mode_v <- v[valid_idx][mode_idx]
      abline(v = mode_v, col = "red", lty = 2, lwd = 1.5)
      legend("topright", 
             legend = paste("Mode ≈", round(mode_v, 3)),
             col = "red", lty = 2, lwd = 1.5, bty = "n")
    } else {
      plot(0, 0, type = "n", 
           xlab = "v", ylab = "Density",
           main = "No valid density values (numerical overflow/underflow)")
      text(0, 0, "Try different parameter values", col = "red", cex = 1.5)
    }
  })
  
  output$paramInfo <- renderText({
    paste(
      "Current Parameters:",
      paste("K =", input$K),
      paste("n =", input$n),
      paste("a =", input$a),
      paste("σ =", input$sigma),
      paste("τ =", input$tau),
      paste("v range = [", input$v_range[1], ",", input$v_range[2], "]"),
      sep = "\n"
    )
  })
}

shinyApp(ui = ui, server = server)