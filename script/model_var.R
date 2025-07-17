# Clean environment
rm(list = ls())

# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  vars,
  readr,
  dplyr,
  tidyr,
  ggplot2
)

# ----------------------------------------------------------------------------
# 1. Load and Prepare Data
# ----------------------------------------------------------------------------

# Load main dataset
data <- readr::read_csv("data/processed/final_data.csv") |>
  # Filter for the correct period, as suggested
  dplyr::filter(ref.date >= as.Date("2013-01-01") & ref.date <= as.Date("2024-12-31"))

# Load external instrument data (assuming this is used for shock identification)
# instrument <- readr::read_csv("data/processed/instrumento.csv")

# Define the core variables for all VAR models
core_vars <- c(
  "producao_industrial" = "producao_industria_geral",
  "ipca" = "price_ipca_cheio",
  "selic" = "juros_selic"
)

# Define the asset variables to be tested, one by one
asset_vars <- c(
  "cambio" = "cambio_usd_venda_ponta",
  "ibrx100" = "asset_ibrx_100",
  "ida" = "retorno_mensal",
  "spread_credito" = "spread_credito_livre_total_pf"
  # Add other asset variables here if you want
)

# List to store the IRF results from each VAR
all_irfs <- list()

# Set fixed VAR parameters as per the research plan
p <- 6 # Fixed lag order
h <- 48 # Horizon for IRFs
ci <- 0.95 # Confidence interval

cat("Starting VAR estimation loop for each asset variable...\n")

# ----------------------------------------------------------------------------
# 2. Loop to Estimate VAR for Each Asset
# ----------------------------------------------------------------------------

for (asset_name in names(asset_vars)) {
  
  cat("\nEstimating VAR for asset:", asset_name, "\n")
  
  # Get the actual column name for the current asset
  asset_col <- asset_vars[asset_name]
  
  # Combine core variables with the current asset variable
  vars_to_select <- c("ref.date", setNames(core_vars, names(core_vars)), setNames(asset_col, asset_name))
  
  # Prepare data for this specific VAR
  var_data <- data |>
    dplyr::select(all_of(vars_to_select)) |>
    tidyr::drop_na() |>
    dplyr::select(-ref.date) |>
    # Convert to time series object
    ts(start = c(2013, 1), frequency = 12) # Adjust start date as needed
  
  # Estimate the VAR
  var_model <- vars::VAR(var_data, p = p, type = "const")
  
  # Identify structural shocks using Cholesky decomposition
  # The ordering is important: activity, prices, policy rate, asset price
  irf_result <- vars::irf(var_model,
    impulse = "selic", # Shock from the policy instrument
    n.ahead = h,
    boot = TRUE,
    ci = ci
  )
  
  # Store the results
  all_irfs[[asset_name]] <- irf_result
  
  cat("Finished estimation for:", asset_name, "\n")
}

cat("\nAll VAR models estimated. IRFs are stored in 'all_irfs'.\n")


# ----------------------------------------------------------------------------
# 3. Plot IRFs
# ----------------------------------------------------------------------------

# Now you can plot the results for each model
# Example: Plotting the results for the first asset (cambio)
plot(all_irfs$cambio)

# Example: Plotting the results for the second asset (ibrx100)
plot(all_irfs$ibrx100)

# For a more direct comparison with your SDFM, you'll want to create
# custom plots with ggplot2, combining the results.

cat("Next steps:\n")
cat("1. Review the generated plots for each asset.\n")
cat("2. Create a custom plotting script to compare VAR and SDFM IRFs side-by-side.\n")
