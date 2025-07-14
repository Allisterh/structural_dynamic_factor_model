#' Example Usage of SDFM Implementation
#' 
#' This script demonstrates how to use the cleaned SDFM implementation
#' following Alessi & Kerssenfischer (2019) methodology.

# Load required functions
source("script/model_alessi.R")

# Basic usage - estimates SDFM with default parameters
results <- main_sdfm()

# Access different components
model_estimates <- results$model
impulse_responses <- results$irfs  
processed_data <- results$data

# Custom parameters
custom_results <- main_sdfm(
  r = 5,      # fewer static factors
  q = 3,      # fewer dynamic factors  
  h = 24,     # shorter IRF horizon
  nboot = 500 # fewer bootstrap replications
)

# Generate specific IRF plots
financial_vars <- list(
  c("IBRx-100" = 64),
  c("IMob" = 68),
  c("Yield - 1A" = 50),
  c("USD/BRL" = 39)
)

financial_plot <- plot_irf(results$irfs,
  response_vars = financial_vars,
  shock = 3,
  horizon = 24
)

# Save results
ggplot2::ggsave("img/financial_responses.png", financial_plot,
  width = 8, height = 6, dpi = 300
)
