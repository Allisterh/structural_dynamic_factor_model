rm(list = ls())

source("R/modeling/factor_estimation.R")
source("R/modeling/impulse_responde.R")


# Lendo os dados ----
#' Lendo os dados e ordenando para que fique de acordo com a matrix de restrição
#' dos fatores, seguindo a metodologia de stock & watson (2016)
data <- readr::read_csv("data/processed/final_data.csv") |>
  dplyr::select(-ref.date) |>
  tidyr::drop_na() |>
  dplyr::select(
    dplyr::contains("juros"),
    dplyr::contains("spread"),
    dplyr::contains("price"),
    dplyr::contains("cambio"),
    dplyr::contains("credito"),
    dplyr::contains("veiculos"),
    dplyr::contains("consumo"),
    dplyr::contains("producao"),
    capacidade_instalada_industria,
    dplyr::contains("vendas"),
    dplyr::contains("trab"),
    dplyr::contains("commodity")
  )


# Olhando a base
dplyr::glimpse(data)





# determinando o numero de fatores ----
bai_ng <- bai_ng_criteria(data, standardize = TRUE)
bai_ng$r_hat

amengual <- amengual_watson(data, r = 6, p = 2)
amengual$q_hat

r <- 6
q <- 4


# Estimando os fatores estaticos ----

## Normalizando os dados ----
pca <- prcomp(data, center = TRUE, scale. = TRUE)

## Obtendo os fatores estaticos ----
static_factors <- pca$x[, 1:r]
factors_df <- as.data.frame(static_factors)
colnames(factors_df) <- paste0("Factor", 1:r)

## Obtendo os fatores dinamicos ----
var_factors <- vars::VAR(factors_df, p = 1, type = "none")
dynamic_factor <- residuals(var_factors)




# Identificando os fatores ----

## Para identificar os fatores, vamos agrupar as variaveis

grupos_variaveis <- list(
  atividade_real = c(
    "producao_bens_consumo", "producao_transformacao", "vendas_varejo",
    "capacidade_instalada_industria", "producao_bens_duraveis",
    "producao_bens_nao_duraveis", "producao_bens_capital"
  ),
  precos = c(
    "price_ipca", "price_igp_m", "price_ipc", "price_incc",
    "price_inpc", "price_ipp"
  ),
  politica_monetaria = c(
    "juros_selic", "juros_cdi", "juros_cdb_rdb",
    "spread_icc_juridica", "spread_icc_fisica"
  ),
  cambio = c(
    "cambio_usd", "cambio_eur", "cambio_gbp", "cambio_ars",
    "cambio_jpy"
  ),
  trab = data |>
    dplyr::select(dplyr::contains("commodity")) |>
    colnames()
)


## Ajustar os dados para terem o mesmo número de observações

data_adjusted <- data[(nrow(data) - nrow(dynamic_factor) + 1):nrow(data), ]

order_shocks <- function(dynamic_factor, data_adjusted, grupos_variaveis) {
  var_explained <- matrix(NA, nrow = q, ncol = length(grupos_variaveis))
  colnames(var_explained) <- names(grupos_variaveis)
  rownames(var_explained) <- paste0("Shock", 1:q)

  for (i in 1:q) {
    for (j in 1:length(grupos_variaveis)) {
      vars <- grupos_variaveis[[j]]
      correlations <- sapply(vars, function(var) {
        abs(cor(dynamic_factor[, i], data_adjusted[[var]]))
      })
      var_explained[i, j] <- mean(correlations)
    }
  }
  return(var_explained)
}

## Calcular a matriz de correlações para identificar os choques
shock_order <- order_shocks(dynamic_factor, data_adjusted, grupos_variaveis)
print(shock_order)




# Reordenar os choques estruturais conforme a análise
# Ordem: Shock1 (atividade) -> Shock4 (preços) -> Shock3 (pol monetária) -> Shock2 (câmbio)
new_order <- c(1, 4, 3, 2)
# Para os dynamic_shocks ordenados
dynamic_shocks_ordered <- dynamic_factor[, new_order]


# Aplicar Cholesky nos resíduos ordenados
sigma <- cov(dynamic_shocks_ordered)
chol_matrix <- t(chol(sigma))
structural_shocks <- dynamic_shocks_ordered %*% solve(chol_matrix)

colnames(structural_shocks) <- c("atividade", "preço", "pol_monetaria", "cambio")

# Verificar resultado
print("Matriz de correlação dos choques estruturais:")
print(round(cor(structural_shocks), 4))
