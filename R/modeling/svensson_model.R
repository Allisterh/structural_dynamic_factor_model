#' Calculate Spot Rate Using Svensson Model
#'
#' This function implements the Svensson (1994) model to calculate spot interest rates.
#' The model extends Nelson-Siegel by adding a fourth term to better fit complex yield curves.
#'
#' @param t Numeric. Time to maturity in years.
#' @param beta0 Numeric. Long-term interest rate level parameter.
#' @param beta1 Numeric. Short-term component parameter.
#' @param beta2 Numeric. Medium-term component parameter.
#' @param beta3 Numeric. Second medium-term component parameter.
#' @param tau1 Numeric. First decay parameter (must be positive).
#' @param tau2 Numeric. Second decay parameter (must be positive).
#'
#' @return Numeric. The spot rate for the given maturity and parameters.
#'
#' @examples
#' # Calculate 5-year spot rate
#' rate <- svensson_rate(
#'   t = 5, beta0 = 0.04, beta1 = -0.02,
#'   beta2 = -0.01, beta3 = 0.005,
#'   tau1 = 1.5, tau2 = 4
#' )
#'
#' @references
#' Svensson, L. E. (1994). Estimating and Interpreting Forward Interest Rates:
#' Sweden 1992-1994. NBER Working Paper Series, No. 4871.
#'
#' @export
svensson_rate <- function(t, beta0, beta1, beta2, beta3, tau1, tau2) {
  term1 <- beta0
  term2 <- beta1 * ((1 - exp(-t / tau1)) / (t / tau1))
  term3 <- beta2 * ((1 - exp(-t / tau1)) / (t / tau1) - exp(-t / tau1))
  term4 <- beta3 * ((1 - exp(-t / tau2)) / (t / tau2) - exp(-t / tau2))

  return(term1 + term2 + term3 + term4)
}





#' Fit Svensson Model Parameters to Observed Rates
#'
#' Estimates the parameters of the Svensson model using numerical optimization.
#' Uses L-BFGS-B method to minimize the sum of squared errors between observed
#' and predicted rates.
#'
#' @param maturities Numeric vector. Time to maturity in years for observed rates.
#' @param rates Numeric vector. Observed interest rates corresponding to maturities.
#'
#' @return Numeric vector. Optimized parameters in the order:
#'   \itemize{
#'     \item beta0: Long-term level
#'     \item beta1: Short-term component
#'     \item beta2: Medium-term component
#'     \item beta3: Second medium-term component
#'     \item tau1: First decay parameter
#'     \item tau2: Second decay parameter
#'   }
#'
#' @details
#' The function uses the L-BFGS-B optimization method with the following constraints:
#' \itemize{
#'   \item All beta parameters: [-0.5, 0.5]
#'   \item Both tau parameters: [0, 10]
#' }
#' Initial values are set to 0.1 for betas and 1 for taus.
#'
#' @examples
#' maturities <- c(0.5, 1, 2, 5, 10)
#' rates <- c(0.02, 0.025, 0.03, 0.035, 0.04)
#' params <- fit_svensson(maturities, rates)
#'
#' @export
fit_svensson <- function(maturities, rates) {
  objective <- function(params) {
    beta0 <- params[1]
    beta1 <- params[2]
    beta2 <- params[3]
    beta3 <- params[4]
    tau1 <- params[5]
    tau2 <- params[6]

    predicted <- sapply(maturities, function(t) {
      svensson_rate(t, beta0, beta1, beta2, beta3, tau1, tau2)
    })

    sum((rates - predicted)^2)
  }

  # Definindo limites para os parâmetros
  lower <- c(-0.5, -0.5, -0.5, -0.5, 0, 0)
  upper <- c(0.5, 0.5, 0.5, 0.5, 10, 10)

  # Otimização
  result <- optim(
    par = c(0.1, 0.1, 0.1, 0.1, 1, 1),
    fn = objective,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  return(result$par)
}



#' Generate Fixed Maturity Interest Rate Series Using Svensson Model
#'
#' Creates a time series of interest rates for fixed maturities by fitting
#' the Svensson model to treasury bond data for each date.
#'
#' @param dados_tesouro Data frame. Treasury bond data with columns:
#'   \itemize{
#'     \item ref_date: Reference date
#'     \item matur_date: Maturity date
#'     \item yield_bid: Yield to maturity
#'   }
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item data: Reference date
#'     \item titulo_0_25ano: 3-month rate
#'     \item titulo_1ano: 1-year rate
#'     \item titulo_2ano: 2-year rate
#'     \item titulo_3ano: 3-year rate
#'     \item titulo_5ano: 5-year rate
#'   }
#'
#' @details
#' The function performs the following steps:
#' 1. Calculates time to maturity in years for each bond
#' 2. For each unique date with at least 4 observations:
#'    - Fits the Svensson model to the observed yields
#'    - Generates rates for standard maturities
#' 3. Returns a data frame with fixed maturity rates
#'
#' @note
#' Dates with fewer than 4 observations will have NA values for all maturities
#' as there isn't enough data to fit the model reliably.
#'
#' @examples
#' # Assuming dados_tesouro is a data frame with required columns
#' fixed_rates <- generate_fixed_maturity_series(dados_tesouro)
#'
#' @export
generate_fixed_maturity_series <- function(dados_tesouro) {
  # Preparar os dados
  dados <- dados_tesouro |>
    dplyr::mutate(
      maturity = as.numeric(matur_date - ref_date) / 365,
      yield = yield_bid
    ) |>
    dplyr::filter(maturity > 0)

  # Definir as maturidades desejadas (em anos)
  maturities <- c(0.25, 1, 2, 3, 5) # 3 meses, 1 ano, 2 anos, 3 anos, 5 anos

  # Criar data frame de resultado com as novas colunas
  result_df <- data.frame(
    data = unique(dados$ref_date)
  )

  # Adicionar colunas para cada maturidade
  for (mat in maturities) {
    colname <- paste0("titulo_", gsub("\\.", "_", as.character(mat)), "ano")
    result_df[[colname]] <- NA_real_
  }

  # Para cada data única
  for (current_date in unique(dados$ref_date)) {
    # Filtrar dados para a data atual
    current_data <- dados |>
      dplyr::filter(ref_date == current_date)

    # Se houver pelo menos 4 observações
    if (nrow(current_data) >= 4) {
      # Ajustar o modelo
      params <- fit_svensson(current_data$maturity, current_data$yield)

      # Calcular taxas para todas as maturidades
      for (mat in maturities) {
        colname <- paste0("titulo_", gsub("\\.", "_", as.character(mat)), "ano")
        rate <- svensson_rate(
          mat, params[1], params[2], params[3],
          params[4], params[5], params[6]
        )
        result_df[result_df$data == current_date, colname] <- rate
      }
    }
  }

  return(result_df)
}
