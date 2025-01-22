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


# Olhando a base
dplyr::glimpse(data)

# Determinando o numero de fatores ----
## Determinando r ----

static_factor <- bai_ng_criteria(data)
static_factor$r_hat

scree_plot <- scree_analysis(data)

# Bai & ng sugerem usar o IC2 e tanto o IC2 quanto o scree plot, sugere
# estimar 6 fatores estaticos


## determinando q ----
dynamic_factor <- amengual_watson(data, static_factor$r_hat$IC2)
dynamic_factor$q_hat # q = 4 fatores dinamicos




# Constructing contrains matrix ----

# Primeiro vamos separar as varaveis por grupos
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



# Juntar tudo em uma lista para aplicar a funçao que cria a matriz de restriçao
list_groups <- list(
  monetary_vars,
  real_vars,
  price_vars
  # employment_vars,
  # commodity_vars,
  # exchange_vars
)

named_factor <- constraint_matrix(list_groups, 6)




# Testando o medelo ----

# Neste primeiro momento, vamos testar o modelo com r=q=6 apenas para ver se
# tudo roda direitinho

sdfm <- estimate_sdfm(
  X = data,
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
