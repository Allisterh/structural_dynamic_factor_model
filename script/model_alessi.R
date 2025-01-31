rm(list = ls())

source("R/modeling/factor_estimation.R")
source("R/modeling/impulse_responde.R")


# Lendo os dados ----
data <- readr::read_csv("data/processed/final_data.csv") |>
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


# determinando o numero de fatores ----
bai_ng <- bai_ng_criteria(data, standardize = TRUE)
bai_ng$r_hat


amengual <- amengual_watson(data, r = 7, p = 1, scale = TRUE)
amengual$q_hat


scree_analysis(data)




# Estimando os fatores ----


# Usage example:
r <- 7 # number of static factors
q <- 5 # number of dynamic factors
p <- 1 # VAR lag order

dfm_results <- estimate_dfm(data, r, q, p)




# Computando a IRF ----


# Calcular IRFs
irf_results <- compute_irf_dfm(dfm_results, h = 50, nboot = 800)




# Diagnosticos do modelo ----

## Teste de ariz unitaria em painel dos dados ----


panel_unit_root_tests <- function(data) {
  # Preparar dados em formato de painel
  T <- nrow(data)
  N <- ncol(data)

  panel_data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = as.vector(data)
  )

  # Criar objeto pdata.frame
  panel_data <- plm::pdata.frame(panel_data, index = c("id", "time"))

  # Teste de Maddala e Wu
  mw_test <- plm::purtest(panel_data$y, test = "madwu", exo = "intercept")

  # Teste de Choi
  choi_test <- plm::purtest(panel_data$y, test = "Pm", exo = "intercept")

  # Teste de Levin-Lin-Chu
  llc_test <- plm::purtest(panel_data$y, test = "levinlin", exo = "intercept")

  # Organizar resultados
  results <- list(
    mw = mw_test,
    choi = choi_test,
    llc = llc_test
  )

  # Função de sumário
  summary <- function() {
    cat("\nTestes de Raiz Unitária em Painel:\n")
    cat("================================\n")

    cat("\n1. Teste de Maddala e Wu:\n")
    print(mw_test)

    cat("\n2. Teste de Choi:\n")
    print(choi_test)

    cat("\n3. Teste de Levin-Lin-Chu:\n")
    print(llc_test)
  }

  results$summary <- summary
  return(results)
}

# Usar as funções


# 2. Testes de Raiz Unitária
unit_root_tests <- panel_unit_root_tests(data)
unit_root_tests$summary()




## fatores estaticos (diagnostico) ----
dfm_diagnostics <- function(dfm_results, original_data) {
  # Verificações iniciais
  required_components <- c("static_factors", "static_loadings", "Z")
  if (!all(required_components %in% names(dfm_results))) {
    stop("Componentes necessários não encontrados em dfm_results")
  }

  # Extração dos componentes
  static_factors <- dfm_results$static_factors
  static_loadings <- dfm_results$static_loadings
  Z <- dfm_results$Z

  # 1. Análise de Variância dos Fatores Estáticos
  common_component <- static_factors %*% t(static_loadings)
  total_var <- sum(apply(scale(original_data), 2, var))

  factor_var <- numeric(ncol(static_factors))
  for (i in 1:ncol(static_factors)) {
    common_i <- static_factors[, i, drop = FALSE] %*% t(static_loadings[, i, drop = FALSE])
    factor_var[i] <- sum(apply(common_i, 2, var)) / total_var * 100
  }

  cumulative_var <- cumsum(factor_var)
  total_var_explained <- sum(factor_var)

  # 2. Análise de Componentes Principais
  pca <- prcomp(Z, scale = FALSE)
  pc_var <- pca$sdev^2 / sum(pca$sdev^2) * 100
  pc_cum_var <- cumsum(pc_var)

  # 3. Análise de Correlação Cruzada
  static_cor <- cor(static_factors)
  diag(static_cor) <- 0
  max_cross_cor <- max(abs(static_cor))

  # Estrutura de resultados
  results <- list(
    variance = list(
      individual = factor_var,
      cumulative = cumulative_var,
      total_explained = total_var_explained,
      pc_variance = pc_var,
      pc_cumulative = pc_cum_var
    ),
    correlation = list(
      matrix = static_cor,
      max_correlation = max_cross_cor
    )
  )

  # Função de sumário
  results$summary <- function() {
    cat("\nDiagnósticos do Modelo de Fatores:\n")
    cat("================================\n")

    cat("\n1. Decomposição da Variância (Fatores Estáticos):\n")
    cat(sprintf("Variância Total Explicada: %.2f%%\n", total_var_explained))
    for (i in seq_along(factor_var)) {
      cat(sprintf(
        "- Fator %d: %.2f%% (Acum: %.2f%%)\n",
        i, factor_var[i], cumulative_var[i]
      ))
    }

    cat("\n2. Análise de Componentes Principais:\n")
    for (i in 1:min(5, length(pc_var))) {
      cat(sprintf(
        "- CP %d: %.2f%% (Acum: %.2f%%)\n",
        i, pc_var[i], pc_cum_var[i]
      ))
    }

    cat("\n3. Correlação entre Fatores:\n")
    cat(sprintf("Máxima correlação cruzada: %.3f\n", max_cross_cor))
  }

  return(results)
}




static_diagnostics <- dfm_diagnostics(dfm_results, data)
static_diagnostics$summary()









## Analise de robustez dos fatores estaticos ----
analyze_robustness <- function(data, r_values = c(5:9)) {
  results <- list()

  for (r in r_values) {
    # Converter data para matriz numérica
    data_matrix <- as.matrix(data)
    if (!is.numeric(data_matrix)) {
      stop("Os dados precisam ser numéricos")
    }

    # Ajustar modelo com diferente número de fatores
    dfm_temp <- tryCatch(
      {
        estimate_dfm(data_matrix, r = r, q = min(r - 2, 5), p = 1)
      },
      error = function(e) {
        message("Erro ao estimar modelo com r = ", r)
        message(e)
        return(NULL)
      }
    )

    if (!is.null(dfm_temp)) {
      # Calcular variância explicada
      static_diagnostics <- dfm_diagnostics(dfm_temp, data_matrix)

      results[[paste0("r_", r)]] <- list(
        variance_explained = static_diagnostics$variance$total_explained,
        cumulative_var = static_diagnostics$variance$cumulative,
        model = dfm_temp
      )
    }
  }

  # Verificar se temos resultados
  if (length(results) == 0) {
    stop("Nenhum modelo pôde ser estimado")
  }

  # Criar plot comparativo
  var_explained <- data.frame(
    r = r_values,
    var = sapply(results, function(x) x$variance_explained)
  )

  p <- ggplot2::ggplot(var_explained, ggplot2::aes(x = factor(r), y = var)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Número de Fatores",
      y = "Variância Explicada (%)",
      title = "Variância Explicada por Número de Fatores"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  return(list(
    results = results,
    plot = p
  ))
}



# Garantir que os dados são numéricos
data_numeric <- apply(data, 2, as.numeric)
robustness_analysis <- analyze_robustness(data_numeric)
print(robustness_analysis$plot)


cat("\nAnálise de Robustez - Variância Explicada:\n")
for (r in names(robustness_analysis$results)) {
  cat(sprintf("\n%s: %.2f%%", r, robustness_analysis$results[[r]]$variance_explained))
}



## Analise de correlaçao das cargas fatoriais ----

# 1. Agrupar variáveis por categoria econômica
categorize_variables <- function(var_names) {
  categories <- list(
    consumo = grep("consumo_", var_names, value = TRUE),
    vendas = grep("vendas_", var_names, value = TRUE),
    veiculos = grep("veiculos_", var_names, value = TRUE),
    producao = grep("producao_", var_names, value = TRUE),
    trabalho = grep("trab_", var_names, value = TRUE),
    precos = grep("price_", var_names, value = TRUE),
    cambio = grep("cambio_", var_names, value = TRUE),
    commodities = grep("commodity_", var_names, value = TRUE),
    juros = grep("juros_|titulo_|spread_", var_names, value = TRUE),
    credito = grep("credito_", var_names, value = TRUE),
    ativos = grep("asset_", var_names, value = TRUE),
    outros = grep("capacidade_", var_names, value = TRUE)
  )

  # Criar vetor ordenado com as categorias
  ordered_vars <- unlist(categories)
  categories_vec <- rep(names(categories), sapply(categories, length))
  names(categories_vec) <- ordered_vars

  return(list(
    categories = categories,
    ordered_vars = ordered_vars,
    categories_vec = categories_vec
  ))
}

# 2. Modificar a função de análise dos loadings para incluir categorias
analyze_static_loadings_grouped <- function(dfm_results, data_names) {
  # Categorizar variáveis
  cat_info <- categorize_variables(data_names)

  # Extrair e preparar loadings como antes
  loadings <- as.matrix(dfm_results$static_loadings)
  loadings_df <- data.frame(loadings)
  colnames(loadings_df) <- paste0("Factor", 1:ncol(loadings))
  rownames(loadings_df) <- data_names

  # Preparar dados para o plot com categorias
  loadings_long <- data.frame(
    Variable = rep(rownames(loadings_df), ncol(loadings_df)),
    Factor = rep(colnames(loadings_df), each = nrow(loadings_df)),
    Value = as.vector(as.matrix(loadings_df)),
    Category = rep(cat_info$categories_vec[rownames(loadings_df)], ncol(loadings_df))
  )

  # Criar plot com categorias
  p <- ggplot2::ggplot(loadings_long, ggplot2::aes(
    x = Factor, y = reorder(Variable, -as.numeric(factor(Category))),
    fill = Value
  )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(min(loadings), max(loadings))
    ) +
    ggplot2::facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Fatores Estáticos", y = "", fill = "Loading") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0)
    )

  return(list(
    loadings = loadings_df,
    plot = p,
    categories = cat_info$categories
  ))
}

loadings_grouped <- analyze_static_loadings_grouped(dfm_results, colnames(data))
print(loadings_grouped$plot)










# Bancada de teste ----







p <- plot_irf(irf_results,
  response_vars = response,
  shock = 3, horizon = 20
)
p



# Exemplo de uso
yield <- list(
  c("IDA" = 54),
  c("Yield - 3M" = 49),
  c("Yield - 1A" = 50),
  c("Yield - 2A" = 51),
  c("Yield - 3A" = 52),
  c("Yield - 5A" = 53)
)

finalcial_markets <- list(
  c("IBRx-100" = 64),
  c("Mid Caps" = 63),
  c("Small Caps" = 70),
  c("IFnc" = 66),
  c("IMob" = 68),
  c("IFix" = 65)
)

ex <- list(
  c("selic" = 47),
  c("spread_j" = 55),
  c("spreac_f" = 56),
  c("commodity agro" = 44),
  c("commodity energia" = 46),
  c("ipca" = 33),
  c("ipp" = 38),
  c("USD" = 39)
)

response <- list(
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

p <- plot_irf(irf_results,
  response_vars = response,
  shock = 3, horizon = 20
)
p

# Depois salve usando ggsave
ggplot2::ggsave("img/yield.png", p,
  width = 10, height = 10, # em polegadas
  dpi = 320
)
