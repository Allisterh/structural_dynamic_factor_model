rm(list = ls())

source("R/modeling/factor_estimation.R")
source("R/modeling/impulse_responde.R")

#' Main SDFM Analysis Function
#'
#' @description
#' Complete implementation of Structural Dynamic Factor Model following 
#' Alessi & Kerssenfischer (2019) methodology for Brazilian economic data.
#'
#' @param data_path Character path to processed data file. Default: "data/processed/final_data.csv"
#' @param r Integer number of static factors. Default: 7
#' @param q Integer number of dynamic factors. Default: 5  
#' @param p Integer VAR lag order. Default: 1
#' @param h Integer IRF horizon. Default: 50
#' @param nboot Integer bootstrap replications for IRFs. Default: 800
#'
#' @return List containing SDFM estimation results and IRF analysis
#' 
#' @export
main_sdfm <- function(data_path = "data/processed/final_data.csv",
                      r = 7, q = 5, p = 1, h = 50, nboot = 800) {
  
  # Load and prepare data
  data <- readr::read_csv(data_path) |>
    dplyr::select(-ref.date) |>
    tidyr::drop_na() |>
    dplyr::select(
      dplyr::contains("consumo_"),
      dplyr::contains("vendas_"),
      dplyr::contains("veiculos_"),
      dplyr::contains("producao_"),
      capacidade_instalada_industria,
      dplyr::contains("trab_"),
      dplyr::contains("price_"),
      dplyr::contains("cambio_"),
      dplyr::contains("commodity_"),
      dplyr::contains("juros_"),
      dplyr::contains("titulo_"),
      retorno_mensal,
      dplyr::contains("spread_"),
      dplyr::contains("credito_"),
      dplyr::contains("asset_")
    ) |>
    as.matrix()
  
  colnames(data)[colnames(data) == "retorno_mensal"] <- "ida"
  
  # Estimate SDFM
  dfm_results <- estimate_dfm(data, r, q, p)
  
  # Validate results
  validation <- validate_dfm_results(dfm_results)
  if (length(validation$missing_components) > 0) {
    warning("Missing DFM components: ", paste(validation$missing_components, collapse = ", "))
  }
  
  # Compute IRFs
  irf_results <- compute_irf_dfm(dfm_results, h = h, nboot = nboot)
  
  return(list(
    model = dfm_results,
    irfs = irf_results,
    data = data
  ))
}

# Execute main analysis
sdfm_results <- main_sdfm()

# Generate IRF plots for key economic variables
response_vars <- list(
  c("USD/BRL" = 39),
  c("Spread-J" = 55),
  c("Spread-F" = 56), 
  c("IPCA" = 33),
  c("IPP" = 38),
  c("IBRx-100" = 64),
  c("IMob" = 68),
  c("IDA" = 54),
  c("Yield - 1A" = 50),
  c("Yield - 5A" = 53)
)

irf_plot <- plot_irf(sdfm_results$irfs,
  response_vars = response_vars,
  shock = 3,
  horizon = 20  
)

# Save plot
ggplot2::ggsave("img/irf_monetary_shock.png", irf_plot,
  width = 10, height = 10, dpi = 320
)
