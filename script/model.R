rm(list = ls())

source("R/modeling/factor_estimation.R")
source("R/modeling/impulse_responde.R")


data <- readr::read_csv("data/processed/final_data.csv") |>
  dplyr::select(-ref.date) |>
  tidyr::drop_na() |>
  dplyr::select(
    dplyr::contains("juros"),
    dplyr::contains("credito"),
    dplyr::contains("spread"),
    dplyr::contains("veiculos"),
    dplyr::contains("consumo"),
    dplyr::contains("producao"),
    capacidade_instalada_industria,
    dplyr::contains("vendas"),
    dplyr::contains("price"),
    dplyr::contains("trab"),
    dplyr::contains("commodity"),
    dplyr::contains("cambio")
  )


dplyr::glimpse(data)



# Constructing contrains matrix

monetary_vars <- c(
  "juros_selic", "juros_cdi", "juros_cdb_rdb",
  "credito_agro", "credito_industria_total", "credito_construcao",
  "credito_comercio", "credito_transporte", "credito_pessoa_fisica",
  "spread_icc_juridica", "spread_icc_fisica"
)

real_vars <- c(
  "veiculos_automoveis", "veiculos_leves", "veiculos_caminhoes",
  "veiculos_onibus", "consumo_gasolina", "consumo_glp",
  "consumo_oleo_combustivel", "consumo_oleo_diesel",
  "consumo_demais_derivados", "consumo_alcool",
  "consumo_eletricidade_comecial", "consumo_eletricidade_residencial",
  "consumo_eletricidade_industrial", "consumo_eletricidade_outros",
  "producao_bens_consumo", "producao_bens_duraveis",
  "producao_bens_nao_duraveis", "producao_bens_capital",
  "producao_transformacao", "capacidade_instalada_industria",
  "vendas_varejo", "vendas_servicos"
)

price_vars <- c(
  "price_ipca", "price_igp_m", "price_ipc", "price_incc",
  "price_inpc", "price_ipp"
)

employment_vars <- c(
  "trab_caged", "trab_pop_forca_trab", "trab_pop_ocupada",
  "trab_tx_desemprego", "trab_hrs_trabalhadas_industria",
  "trab_razao_vagas_desempregados", "trab_menos_de_1_mes",
  "trab_de_1_mes_a_menos_de_1_ano", "trab_de_1_ano_a_menos_de_2_anos",
  "trab_x2_anos_ou_mais"
)

commodity_vars <- c(
  "indice_commodity_agro", "indice_commodity_metal",
  "indice_commodity_energia"
)

exchange_vars <- c(
  "cambio_usd", "cambio_eur", "cambio_gbp", "cambio_ars",
  "cambio_jpy"
)




list_groups <- list(
  monetary_vars,
  real_vars,
  price_vars
  # employment_vars,
  # commodity_vars,
  # exchange_vars
)

named_factor <- constraint_matrix(list_groups, 6)


# Bancada de teste ----


#' Plot Structural Impulse Response Functions
#' @param sirf Array of structural impulse responses
#' @param variables Variables to plot
#' @param shocks Shocks to plot
#' @param horizon Number of periods to plot
#' @param ci Confidence interval (0-1)
#' @param se Standard errors matrix (optional)
plot_sirf <- function(
    sirf, variables = NULL, shocks = NULL,
    horizon = NULL, ci = 0.95, se = NULL) {
  # Get dimensions
  n_var <- dim(sirf)[1]
  n_shocks <- dim(sirf)[2]
  h <- dim(sirf)[3]

  # Print SIRF info for debugging
  print("SIRF array dimensions:")
  print(dim(sirf))

  print("Sample SIRF values for first few horizons:")
  for (t in 1:min(5, h)) {
    print(paste("Horizon", t - 1))
    print(sirf[, , t])
  }

  # Set defaults
  if (is.null(variables)) variables <- 1:n_var
  if (is.null(shocks)) shocks <- 1:n_shocks
  if (is.null(horizon)) horizon <- h - 1

  time <- 0:horizon

  # Create long format data
  irf_df <- tidyr::crossing(
    variable = variables,
    shock = shocks,
    horizon = time
  ) |>
    dplyr::mutate(
      irf = purrr::map_dbl(
        1:dplyr::n(),
        ~ sirf[variable[.x], shock[.x], horizon[.x] + 1]
      )
    )

  # Add confidence intervals if SE provided
  if (!is.null(se)) {
    irf_df <- irf_df |>
      dplyr::mutate(
        ci_lower = irf - stats::qnorm(1 - (1 - ci) / 2) *
          se[variable, shock, horizon + 1],
        ci_upper = irf + stats::qnorm(1 - (1 - ci) / 2) *
          se[variable, shock, horizon + 1]
      )
  }

  # Create plots
  plots <- list()
  for (s in shocks) {
    p <- ggplot2::ggplot(
      dplyr::filter(irf_df, shock == s),
      ggplot2::aes(x = horizon, y = irf)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "gray50"
      ) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      {
        if (!is.null(se)) {
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
            alpha = 0.2,
            fill = "steelblue"
          )
        }
      } +
      ggplot2::facet_wrap(~variable, scales = "free_y") +
      ggplot2::labs(
        x = "Horizon",
        y = "Response",
        title = paste("Responses to Shock", s)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        strip.background = ggplot2::element_rect(fill = "white")
      )

    plots[[s]] <- p
  }

  return(plots)
}



# Testando o medelo

x <- scale(data)

sdfm <- estimate_sdfm(
  X = x,
  r = 6,
  p_var = 2,
  constraint_matrix = named_factor,
  group_vars = list_groups,
  horizon = 60
)



plots <- plot_sirf(
  sirf = sdfm$sirf,
  variables = 1:6, # Primeiras 6 variáveis
  shocks = 1,
  horizon = 24
)

# Ver resposta das variáveis ao primeiro choque (monetário)
plots[[1]]


print(sdfm$H) # matriz H
