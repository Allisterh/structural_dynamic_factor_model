# Getting auxiliary functions ----
source("R/data_download/bcb.R")
source("R/data_download/exchange.R")


# Taxa de cambio ----

moedas <- c("USD", "EUR", "GBP", "ARS", "JPY", "CNY", "RUB")

# problema na cotação do rublo(RUB) e yuan(CNY)
cambio <- download_cambio(moedas) |> dplyr::rename(ref.date = date)



# Dados do tesouro direto ----

# treasure <- GetTDData::td_get("LTN", 2010, 2020)

# dplyr::glimpse(treasure)


# Dados do Banco Central do Brasil ----

# Primeiro vamos pegar os dados mensais
## Dados referente a juros ----
vec_juros <- c(
  "selic" = 4189,
  "cdi" = 4392,
  "cdb_rdb" = 28618 # Taxa média acumulada no mês de instituições financeiras
)

juros <- download_bcb_data(vec_juros, parallel = TRUE) |>
  dplyr::rename_with(~ paste0("juros_", .), -ref.date)


## Dados de crédito ----
vec_credito <- c(
  "credito_agro" = 22027,
  "credito_industria_total" = 22043,
  "credito_construcao" = 22030,
  # "credito_siup" = 22034,
  "credito_comercio" = 22036,
  "credito_transporte" = 22037,
  "credito_pessoa_fisica" = 22050,
  "spread_icc_juridica" = 27444,
  "spread_icc_fisica" = 27445
)

credito <- download_bcb_data(vec_credito, parallel = TRUE)


## Dados sobre atividade economica ----

vec_ativ_economica <- c(
  "veiculos_automoveis" = 7384,
  "veiculos_leves" = 7385,
  "veiculos_caminhoes" = 7386,
  "veiculos_onibus" = 7387,
  "consumo_gasolina" = 1393,
  "consumo_glp" = 1394,
  "consumo_oleo_combustivel" = 1395,
  "consumo_oleo_diesel" = 1396,
  "consumo_demais_derivados" = 1397,
  "consumo_alcool" = 1401,
  "consumo_eletricidade_comecial" = 1402,
  "consumo_eletricidade_residencial" = 1403,
  "consumo_eletricidade_industrial" = 1404,
  "consumo_eletricidade_outros" = 1405,
  "producao_bens_consumo" = 21865,
  "producao_bens_duraveis" = 21866,
  "producao_bens_nao_duraveis" = 21867,
  "producao_bens_capital" = 21863,
  "producao_transformacao" = 21862,
  "capacidade_instalada_industria" = 28554,
  "vendas_varejo" = 1455,
  "vendas_servicos" = 23982
)


ativ_economica <- download_bcb_data(vec_ativ_economica, parallel = TRUE)




## Dados sobre emprego ----

vec_emprego <- c(
  "caged" = 28763,
  "pop_forca_trab" = 24378,
  "pop_ocupada" = 28543,
  "tx_desemprego" = 24369,
  "hrs_trabalhadas_industria" = 28556
)


emprego <- download_bcb_data(vec_emprego, parallel = TRUE) |>
  dplyr::mutate(
    razao_vagas_desempregados = caged / ((tx_desemprego / 100) * (pop_forca_trab * 1000))
  ) |>
  tidyr::drop_na() |>
  dplyr::rename_with(~ paste0("trab_", .), -ref.date)





## Dados sobre inflacao ----

vec_inflacao <- c(
  "ipca" = 433,
  # "ipa" = 255,
  "igp_m" = 189,
  "ipc" = 191,
  "incc" = 192,
  "inpc" = 188
)

ipp <- ipeadatar::ipeadata("IPP12_IPPCG12") |>
  dplyr::select(ref.date = date, ipp = value)


inflacao <- download_bcb_data(vec_inflacao, parallel = TRUE) |>
  dplyr::left_join(ipp, by = "ref.date") |>
  tidyr::drop_na() |>
  dplyr::rename_with(~ paste0("price_", .), -ref.date)




## Dados sobre commodities ----

vec_commodity <- c(
  "indice_commodity_agro" = 27575,
  "indice_commodity_metal" = 27576,
  "indice_commodity_energia" = 27577
)

commodity <- download_bcb_data(vec_commodity, parallel = TRUE)



# Dados trimestrais ----
# aqui vamos ter que fazer uma interpolação dos dados para transformar para mensal

tempo_procura_trab <- sidrar::get_sidra(1616, variable = 4092, geo = "Brazil", period = "all") |>
  dplyr::select(Trimestre, `Tempo de procura de trabalho`, Valor) |>
  dplyr::filter(`Tempo de procura de trabalho` != "Total") |>
  tidyr::pivot_wider(
    names_from = `Tempo de procura de trabalho`,
    values_from = Valor
  ) |>
  janitor::clean_names() |>
  dplyr::mutate(
    trimestre = stringr::str_extract(trimestre, "\\d{1}º trimestre \\d{4}"),
    # Primeiro vamos criar o formato ano-trimestre (exemplo: 2012-Q1)
    trimestre = stringr::str_replace(trimestre, "(\\d{1})º trimestre (\\d{4})", "\\2-Q\\1"),
    # Agora convertemos para data usando yearquarter do tsibble
    trimestre = tsibble::yearquarter(trimestre)
  )


# Criar sequência mensal
dados_mensais <- tempo_procura_trab |>
  # Criar datas mensais
  dplyr::mutate(
    data_inicio = as.Date(trimestre),
    meses = purrr::map(data_inicio, ~ seq.Date(
      from = .,
      by = "month",
      length.out = 3
    ))
  ) |>
  tidyr::unnest(meses) |>
  dplyr::select(-trimestre, -data_inicio)



# Aplicar interpolação spline para cada coluna
colunas_interpoladas <- dados_mensais |>
  dplyr::mutate(
    across(
      .cols = -meses,
      .fns = ~ stats::spline(
        x = as.numeric(meses),
        y = .,
        n = length(meses)
      )$y
    )
  )

# Formatar resultado final
resultado_mensal <- colunas_interpoladas |>
  dplyr::mutate(
    mes = lubridate::floor_date(meses, "month")
  ) |>
  dplyr::select(ref.date = mes, everything(), -meses) |>
  dplyr::arrange(ref.date) |>
  dplyr::rename_with(~ paste0("trab_", .), -ref.date)





# Dados da B3 ----

path_b3_csv <- list.files(
  path = "data/raw/indexes/",
  pattern = "*.csv",
  full.names = TRUE
)

# Método 1: Usando map do purrr para criar uma lista de dataframes



indices <- path_b3_csv |>
  purrr::map(~ readr::read_csv(
    .x,
    locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
    show_col_types = FALSE
  )) |>
  purrr::set_names(
    c(
      "asset_ibrx100", "asset_imat", "asset_imat", "asset_ifnc",
      "asset_ifix", "asset_imob", "asset_MLCX", "asset_smll"
    )
  ) |>
  purrr::map(~ .x |>
    janitor::clean_names() |>
    dplyr::select(ref.date = data, ultimo)) |>
  purrr::map(~ .x |>
    dplyr::mutate(ref.date = lubridate::dmy(ref.date))) |>
  purrr::imap(~ .x |>
    dplyr::rename(!!.y := ultimo)) |>
  purrr::reduce(dplyr::left_join, by = "ref.date")



# Dados da AMBIMA ----

path_ambima_csv <- list.files(
  path = "data/raw/mp_index/",
  pattern = "*.csv",
  full.names = TRUE
)



mp_indexes <- path_ambima_csv |>
  purrr::map(~ readr::read_csv(
    .x,
    show_col_types = FALSE
  )) |>
  purrr::set_names(
    c(
      "ida_di", "idk_1a", "idk_2a", "idk_3a",
      "idk_3m", "idk_5a"
    )
  ) |>
  purrr::map(~ .x |>
    janitor::clean_names() |>
    dplyr::select(ref.date = data_de_referencia, numero_indice)) |>
  purrr::map(~ .x |>
    dplyr::mutate(ref.date = lubridate::dmy(ref.date))) |>
  purrr::imap(~ .x |>
    dplyr::rename(!!.y := numero_indice)) |>
  purrr::reduce(dplyr::left_join, by = "ref.date")







# Juntando tudo em apenas um df ----

all_dfs <- list(
  cambio = cambio,
  juros = juros,
  credito = credito,
  ativ_economica = ativ_economica,
  emprego = emprego,
  inflacao = inflacao,
  commodity = commodity,
  tempo_procura = resultado_mensal,
  indices = indices,
  mp_indexes = mp_indexes
)

merged_df <- all_dfs |>
  purrr::reduce(dplyr::left_join, by = "ref.date") |>
  dplyr::arrange(ref.date) |>
  tidyr::drop_na()



# dados vão até jan/2024 por causa do cdb_rdb
# readr::write_csv(merged_df, "data/raw/raw_data.csv")
