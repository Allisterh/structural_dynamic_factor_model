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
    dplyr::contains("spread_"),
    dplyr::contains("credito_"),
    dplyr::contains("asset_")
  ) |>
  as.matrix()


# determinando o numero de fatores ----
bai_ng <- bai_ng_criteria(data, standardize = TRUE)



amengual <- amengual_watson(data, r = 7, p = 1, scale = TRUE)

r <- 7
q <- 5


# Estimando os fatores estaticos ----




estimate_static_factors <- function(data, r) {
  # Dimensões
  T <- nrow(data)

  # Primeira diferença
  y <- diff(data)

  # Desvio padrão das diferenças
  sy <- apply(y, 2, sd)

  # Standardização das diferenças
  yy <- scale(y)

  # Detrending dos dados originais
  regX <- cbind(1, 1:T)
  beta <- solve(crossprod(regX)) %*% crossprod(regX, data)
  X <- data - regX %*% beta

  # Standardização BLL
  Z <- sweep(X, 2, sy, "/")

  # PCA nos dados transformados
  pca_result <- prcomp(Z, scale = FALSE)

  # Extrair fatores e loadings
  factors <- pca_result$x[, 1:r]
  loadings <- pca_result$rotation[, 1:r]

  return(list(
    factors = factors,
    loadings = loadings,
    sy = sy,
    Z = Z
  ))
}



# Função para construir a matriz companion
construct_companion <- function(coef_mat, n_vars, n_lags) {
  companion_size <- n_vars * n_lags
  companion_matrix <- matrix(0, companion_size, companion_size)

  # Preencher com coeficientes VAR (excluir constante)
  companion_matrix[1:n_vars, ] <- t(coef_mat)

  # Preencher com identidade para lags adicionais
  if (n_lags > 1) {
    companion_matrix[(n_vars + 1):companion_size, 1:(companion_size - n_vars)] <-
      diag(1, companion_size - n_vars)
  }

  return(companion_matrix)
}

# Correção de Kilian reescrita
kilian_correction <- function(A, u, n_vars, n_lags, n_obs) {
  # Construir SIGMA (matriz de covariância dos resíduos)
  SIGMA <- crossprod(u) / (n_obs - n_vars * n_lags - 1)

  # Calcular autovalores da matriz companion original
  peigen <- eigen(A)$values

  # Inicializar soma dos autovalores
  companion_size <- n_vars * n_lags
  I <- diag(companion_size)
  B <- t(A)
  sumeig <- matrix(0, companion_size, companion_size)

  # Calcular soma dos autovalores (seguindo MATLAB)
  for (h in 1:companion_size) {
    sumeig <- sumeig + peigen[h] * solve(I - peigen[h] * B)
  }

  # Calcular SIGMA_Y usando kronecker
  vecSIGMAY <- solve(diag(companion_size^2) - kronecker(A, A)) %*%
    as.vector(SIGMA)
  SIGMAY <- matrix(vecSIGMAY, companion_size, companion_size)

  # Calcular viés inicial
  bias <- SIGMA %*% (solve(I - B) + B %*% solve(I - B %*% B) + sumeig) %*%
    solve(SIGMAY)
  bias <- -bias / n_obs

  # Ajuste iterativo para garantir estacionariedade
  delta <- 1
  bcstab <- 9

  while (bcstab >= 1 && delta > 0) {
    bcA <- A - delta * bias

    if (max(abs(eigen(bcA)$values)) >= 1) {
      bcstab <- 1
    } else {
      bcstab <- 0
    }

    delta <- delta - 0.01

    if (delta <= 0) {
      bcstab <- 0
    }
  }

  return(bcA)
}

# Função para estimar VAR com correção
estimate_corrected_var <- function(data, p) {
  # Dimensões
  T <- nrow(data)
  K <- ncol(data)

  # Construir matriz de regressores
  Y <- data[(p + 1):T, ]
  X <- matrix(0, T - p, K * p)

  for (i in 1:p) {
    X[, ((i - 1) * K + 1):(i * K)] <- data[(p - i + 1):(T - i), ]
  }

  # Estimar coeficientes
  beta <- solve(crossprod(X)) %*% crossprod(X, Y)

  # Calcular resíduos
  resid <- Y - X %*% beta

  # Construir matriz companion
  A <- construct_companion(beta, K, p)

  # Aplicar correção de Kilian
  A_corrected <- kilian_correction(A, resid, K, p, T)

  # Extrair coeficientes corrigidos
  beta_corrected <- t(A_corrected[1:K, ])

  return(list(
    coefficients = beta_corrected,
    residuals = resid,
    companion = A_corrected
  ))
}

# Função para computar resíduos
compute_residuals <- function(data, beta, p) {
  T <- nrow(data)
  K <- ncol(data)

  # Construir matriz de regressores
  Y <- data[(p + 1):T, ]
  X <- matrix(0, T - p, K * p)

  for (i in 1:p) {
    X[, ((i - 1) * K + 1):(i * K)] <- data[(p - i + 1):(T - i), ]
  }

  # Calcular resíduos
  residuals <- Y - X %*% beta

  return(residuals)
}


estimate_dynamic_factors <- function(var_residuals, q, r) {
  if (q == r) {
    # Caso especial: q = r
    K <- diag(r)
    M <- diag(r)
    eta <- var_residuals # Não precisamos transformar
  } else {
    # Calcular matriz de covariância dos resíduos
    sigma_u <- cov(var_residuals)

    # Extrair os q maiores autovalores/autovetores
    eig <- eigen(sigma_u)
    idx <- order(abs(eig$values), decreasing = TRUE)[1:q]

    # Ordenar K e M pelos maiores autovalores
    K <- eig$vectors[, idx]
    M <- diag(sqrt(eig$values[idx]))

    # Calcular fatores dinâmicos
    eta <- var_residuals %*% K # Não multiplica por solve(M)
  }

  # Normalizar fatores
  eta_std <- scale(eta, scale = TRUE)

  return(list(
    factors = eta_std, # Fatores normalizados
    factors_raw = eta, # Fatores não normalizados
    K = K,
    M = M
  ))
}



estimate_dfm <- function(data, r, q, p) {
  # 1. Estimação dos fatores estáticos
  static_result <- estimate_static_factors(data, r)

  # 2. Estimação do VAR com correção de Kilian
  var_result <- estimate_corrected_var(static_result$factors, p)

  # 3. Estimação dos fatores dinâmicos
  dynamic_result <- estimate_dynamic_factors(var_result$residuals, q, r)

  # 4. Preparar componentes para IRF
  companion <- var_result$companion

  # Retornar todos os componentes necessários
  return(list(
    static_factors = static_result$factors,
    static_loadings = static_result$loadings,
    var_coefficients = var_result$coefficients,
    var_residuals = var_result$residuals,
    companion_matrix = companion,
    dynamic_factors = dynamic_result$factors,
    factors_raw = dynamic_result$factors_raw, # Corrigido aqui
    dynamic_loadings = dynamic_result$K,
    dynamic_scaling = dynamic_result$M,
    data_sd = static_result$sy,
    Z = static_result$Z
  ))
}

# Usage example:
r <- 7 # number of static factors
q <- 5 # number of dynamic factors
p <- 1 # VAR lag order

dfm_results <- estimate_dfm(data, r, q, p)




# Computando a IRF ----








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

  p <- ggplot(var_explained, aes(x = factor(r), y = var)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(
      x = "Número de Fatores",
      y = "Variância Explicada (%)",
      title = "Variância Explicada por Número de Fatores"
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  return(list(
    results = results,
    plot = p
  ))
}



# Garantir que os dados são numéricos
data_numeric <- apply(data, 2, as.numeric)
robustness_analysis <- analyze_robustness(data_numeric)
print(robustness_analysis$plot)

# Imprimir sumário
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
  p <- ggplot(loadings_long, aes(
    x = Factor, y = reorder(Variable, -as.numeric(factor(Category))),
    fill = Value
  )) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(min(loadings), max(loadings))
    ) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    labs(x = "Fatores Estáticos", y = "", fill = "Loading") +
    theme(
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      strip.text.y = element_text(angle = 0)
    )

  return(list(
    loadings = loadings_df,
    plot = p,
    categories = cat_info$categories
  ))
}

loadings_grouped <- analyze_static_loadings_grouped(dfm_results, colnames(data))
print(loadings_grouped$plot)





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
