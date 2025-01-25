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
    dplyr::contains("credito"),
    dplyr::contains("veiculos"),
    dplyr::contains("consumo"),
    dplyr::contains("producao"),
    capacidade_instalada_industria,
    dplyr::contains("vendas"),
    dplyr::contains("price"),
    dplyr::contains("juros"),
    dplyr::contains("spread"),
    dplyr::contains("cambio"),
    dplyr::contains("trab"),
    dplyr::contains("commodity")
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
  "spread_icc_juridica", "spread_icc_fisica"
)

price_vars <- c(
  "price_ipca", "price_igp_m", "price_ipc", "price_incc",
  "price_inpc", "price_ipp"
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
  price_vars
  # employment_vars,
  # commodity_vars,
  # exchange_vars
)

named_factor <- constraint_matrix(list_groups, 6)



# Testando o medelo ----

# Neste primeiro momento, vamos testar o modelo com r=q=6 apenas para ver se
# tudo roda direitinho




#### teste 1 ----

calculate_sirf <- function(
    factor_est, var_lags = 4, horizon = 20, q = NULL,
    cumulative = FALSE, n_ar_lags = 4, verbose = FALSE) {
  n_series <- nrow(factor_est$loadings)
  r <- ncol(factor_est$loadings)
  if (is.null(q)) q <- r

  if (q > r) stop("q must be less than or equal to r")
  if (horizon < 1) stop("horizon must be positive")

  # AR models for idiosyncratic components
  ar_coefs <- matrix(0, n_series, n_ar_lags)
  ar_sds <- numeric(n_series)

  idio_resids <- factor_est$residuals
  for (i in 1:n_series) {
    ar_model <- stats::ar(idio_resids[, i], aic = FALSE, order.max = n_ar_lags)
    ar_coefs[i, ] <- c(ar_model$ar, rep(0, n_ar_lags - length(ar_model$ar)))
    ar_sds[i] <- sqrt(ar_model$var.pred)
  }

  # VAR and G matrix estimation
  if (r > q) {
    G_est <- estimate_G(factor_est$factors, q, var_lags)
    G <- G_est$G
    eta_t <- G_est$eta
    var_model <- G_est$var_model
  } else {
    var_model <- vars::VAR(factor_est$factors, p = var_lags)
    G <- diag(r)
    eta_t <- residuals(var_model)
  }

  # VAR companion form
  var_coef <- vars::Acoef(var_model)
  M <- matrix(0, r * var_lags, r * var_lags)
  M[1:r, 1:(r * var_lags)] <- t(matrix(unlist(var_coef), nrow = r))
  if (var_lags > 1) {
    M[(r + 1):(r * var_lags), 1:(r * var_lags - r)] <- diag(r * (var_lags - 1))
  }

  # Unit effect normalization via Cholesky
  Sigma_eta <- cov(eta_t)
  H_chol <- t(chol(Sigma_eta))
  D <- diag(diag(H_chol))
  H <- H_chol %*% solve(D)

  if (!all(abs(diag(H) - 1) < 1e-10)) {
    stop("H matrix does not have unit diagonal elements")
  }

  # Structural shocks and orthogonality check
  structural_shocks <- eta_t %*% solve(H)
  shock_cor <- cor(structural_shocks)
  if (!all(abs(shock_cor[lower.tri(shock_cor)]) < 0.1)) {
    warning("Structural shocks may not be sufficiently orthogonal")
    if (verbose) print(shock_cor)
  }

  # Calculate IRFs
  irf_factors <- array(0, c(r, q, horizon))
  irf_variables <- array(0, c(n_series, q, horizon))

  if (cumulative) {
    cum_irf_variables <- array(0, c(n_series, q, horizon))
  }

  for (i_shock in 1:q) {
    x <- matrix(0, r * var_lags, 1)
    x[1:r] <- G %*% H[, i_shock] # impacto inicial do choque i

    for (h in 1:horizon) {
      # Resposta dos fatores: pegar apenas os primeiros r elementos
      irf_factors[, i_shock, h] <- x[1:r]

      # Resposta das variáveis: multiplicar loadings pela resposta dos fatores
      irf_variables[, i_shock, h] <- factor_est$loadings %*% x[1:r]

      # Atualizar para próximo período usando a forma companheira
      x <- M %*% x
    }
  }

  # Check impact effects
  impact_effects <- irf_variables[, , 1]
  if (!all(abs(diag(impact_effects[1:q, 1:q]) - 1) < 1e-10)) {
    warning("Impact effects may not respect unit normalization")
    if (verbose) {
      cat("\nImpact effects matrix (first q x q block):\n")
      print(impact_effects[1:q, 1:q])
    }
  }

  # Variance decomposition
  vd_common <- array(0, c(n_series, q, horizon))
  vd_idio <- matrix(0, n_series, horizon)

  for (i in 1:n_series) {
    M_idio <- matrix(0, n_ar_lags, n_ar_lags)
    M_idio[1, ] <- ar_coefs[i, ]
    if (n_ar_lags > 1) {
      M_idio[2:n_ar_lags, 1:(n_ar_lags - 1)] <- diag(n_ar_lags - 1)
    }

    x_idio <- matrix(0, n_ar_lags, 1)
    x_idio[1] <- ar_sds[i]
    idio_var <- numeric(horizon)

    for (h in 1:horizon) {
      idio_var[h] <- x_idio[1]^2
      x_idio <- M_idio %*% x_idio
    }
    idio_var_cumsum <- cumsum(idio_var)

    for (h in 1:horizon) {
      common_var <- sapply(1:q, function(j) {
        sum(irf_variables[i, j, 1:h]^2)
      })
      total_var <- sum(common_var) + idio_var_cumsum[h]

      vd_common[i, , h] <- common_var / total_var
      vd_idio[i, h] <- idio_var_cumsum[h] / total_var
    }
  }

  validation <- list(
    H_diagonal = all(abs(diag(H) - 1) < 1e-10),
    shock_orthogonality = all(abs(shock_cor[lower.tri(shock_cor)]) < 0.1),
    G_dimensions = if (r > q) all(dim(G) == c(r, q)) else TRUE,
    impact_normalization = all(abs(diag(impact_effects[1:q, 1:q]) - 1) < 1e-10)
  )

  output <- list(
    irf = irf_variables,
    irf_factors = irf_factors,
    variance_decomp_common = vd_common,
    variance_decomp_idio = vd_idio,
    structural_shocks = structural_shocks,
    structural_cov = H %*% t(H),
    state_space = list(M = M, H = H, G = G),
    ar_coefficients = ar_coefs,
    ar_standard_deviations = ar_sds,
    validation = validation,
    impact_effects = impact_effects[1:q, 1:q]
  )

  if (cumulative) {
    output$cum_irf <- cum_irf_variables
  }

  return(output)
}



# Estimação dos fatores
factors <- estimate_factor(
  X = data, n_factors = 6,
  constraint_info = named_factor
)



results <- calculate_sirf(factors,
  var_lags = 2,
  horizon = 40, q = NULL,
  verbose = FALSE
)





# Diagnostico da funçao IRF -----


# # 1. Verificar dimensões e estrutura
# dims_check <- function(irf_results, n_series, n_factors, horizon) {
#   list(
#     irf_dims = dim(irf_results$irf) == c(n_series, n_factors, horizon),
#     vd_dims = dim(irf_results$variance_decomp_common) == c(n_series, n_factors, horizon),
#     vd_sum = all.equal(
#       apply(irf_results$variance_decomp_common, c(1, 3), sum) +
#         irf_results$variance_decomp_idio,
#       matrix(1, n_series, horizon)
#     )
#   )
# }

# # 2. Verificar propriedades dos choques estruturais
# shock_diagnostics <- function(irf_results) {
#   shocks <- irf_results$structural_shocks
#   list(
#     mean_zero = colMeans(shocks),
#     unit_var = apply(shocks, 2, var),
#     normality = apply(shocks, 2, shapiro.test)
#   )
# }

# # 3. Verificar consistência das IRFs
# irf_checks <- function(irf_results) {
#   # Impact effects devem respeitar normalização
#   impact_effects <- irf_results$irf[, , 1]

#   # Persistência das IRFs
#   persistence <- apply(irf_results$irf, c(1, 2), function(x) sum(abs(x)))

#   # Variação cross-section das IRFs
#   cross_variation <- apply(irf_results$irf, c(2, 3), sd)

#   list(
#     impact = impact_effects,
#     persistence = persistence,
#     cross_var = cross_variation
#   )
# }

# # 4. Decomposição de variância
# vd_checks <- function(irf_results) {
#   # Importância relativa dos componentes comum vs idiossincrático
#   relative_importance <- rowMeans(irf_results$variance_decomp_idio)

#   # Contribuição dos choques ao longo do tempo
#   shock_importance <- apply(irf_results$variance_decomp_common, 2, rowMeans)

#   list(
#     idio_share = relative_importance,
#     shock_shares = shock_importance
#   )
# }




# # 1. Verificar normalização unit effect
# all.equal(diag(irf_results$state_space$H), rep(1, ncol(irf_results$state_space$H)))

# # 2. Verificar ortogonalidade dos choques
# cor(irf_results$structural_shocks) # Deve ser próxima da identidade

# # 4. Resposta contemporânea deve refletir restrições
# irf_results$irf[1:3, 1:3, 1] # Verificar diagonal unitária



# dims_check(irf_results, 57, 6, 50)

# shock_diagnostics(irf_results)

# irf_checks(irf_results)

# vd_checks(irf_results)









#### plots ----


# Limpar qualquer configuração gráfica anterior
graphics.off()

# Configuração do layout dos gráficos
par(mfrow = c(2, 3))

# Loop para plotar as respostas
for (i in seq_len(dim(irf_results$irf_factors)[1])) {
  plot(1:dim(irf_results$irf_factors)[3], irf_results$irf_factors[i, 1, ], # Corrigido
    type = "l",
    main = paste("Resposta do fator", i, "ao Choque 1"),
    xlab = "Horizonte",
    ylab = "Resposta"
  )
}




# Bancada de teste ----

#' Calculate Impulse Response Functions for State Space Model
#'
#' @param coef List containing Q, M, G matrices of state space model:
#'             y(t) = Q*z(t)
#'             z(t) = M*z(t-1) + G*u(t)
#'             var(u(t)) = I
#' @param i Vector of shocks to compute IRFs for (0 for all shocks)
#' @param hor Horizon for impulse responses
#' @return Array with dimensions [ndim, nshocks, hor] containing IRFs
#'         where irf[i,j,t] is response of variable i to shock j at horizon t-1
#' @import dplyr
varirf_ss <- function(coef, i, hor) {
  # Extract matrices from coef list
  Q <- coef$Q
  M <- coef$M
  G <- coef$G

  # Get dimensions
  ndim <- nrow(Q)

  # If i=0, compute IRFs for all shocks
  if (i == 0) {
    nshocks <- ncol(G)
    i <- 1:nshocks
  } else {
    nshocks <- length(i)
  }

  # Initialize output array
  irf <- array(NA, dim = c(ndim, nshocks, hor))

  # Loop over shocks
  for (j in 1:nshocks) {
    i_shock <- i[j]
    x <- G[, i_shock, drop = FALSE]

    # Loop over horizons
    for (i_hor in 1:hor) {
      y <- Q %*% x
      irf[, j, i_hor] <- y
      x <- M %*% x
    }
  }

  return(irf)
}


#' Estimate VAR on factors and return companion form matrices
#' @param y Matrix of factors (t x n)
#' @param var_par List with parameters
varest_r <- function(y, var_par) {
  # Extract parameters
  nlag <- var_par$nlag
  icomp <- var_par$icomp
  iconst <- var_par$iconst

  # Estimate VAR
  if (iconst == 1) {
    type <- "const"
  } else {
    type <- "none"
  }

  var_est <- vars::VAR(y, p = nlag, type = type)

  # Check VAR stability
  roots <- abs(roots(var_est))
  if (any(roots >= 1)) {
    warning("VAR may be unstable - some roots are on or outside unit circle")
  }

  # Get coefficients and residuals
  betahat <- t(sapply(coef(var_est), function(x) x[, 1]))
  resid <- residuals(var_est)
  T <- nrow(resid)
  seps <- crossprod(resid) / (T - nrow(betahat))

  if (icomp > 0) {
    ns <- ncol(y)

    # Extract VAR coefficients
    if (iconst == 0) {
      b <- t(betahat)
      const_coef <- matrix(0, ns, 1)
    } else {
      b <- t(betahat[, -1])
      const_coef <- betahat[, 1]
    }

    # Create companion matrix
    comp_size <- ns * nlag
    comp <- matrix(0, comp_size, comp_size)
    comp[1:ns, ] <- b

    if (comp_size > ns) {
      comp[(ns + 1):comp_size, 1:(comp_size - ns)] <- diag(comp_size - ns)
    }

    # Normalize innovations to have unit variance
    seps_chol <- t(chol(seps))
    seps_chol <- seps_chol / sqrt(diag(seps))

    # Create state space matrices
    M <- comp
    Q <- cbind(diag(ns), matrix(0, ns, ncol(M) - ns))
    G <- rbind(seps_chol, matrix(0, nrow(M) - ns, ns))

    if (icomp == 2 && iconst == 1) {
      G <- rbind(G, matrix(0, 1, ncol(G)))
      Q <- cbind(Q, matrix(0, nrow(Q), 1))
      M <- cbind(M, matrix(0, nrow(M), 1))
      M <- rbind(M, matrix(0, 1, ncol(M)))
      M[nrow(M), ncol(M)] <- 1
      M[1:ns, ncol(M)] <- const_coef
    }
  } else {
    Q <- M <- G <- NULL
  }

  return(list(
    betahat = betahat,
    seps = seps,
    resid = resid,
    roots = roots,
    coef = list(
      Q = Q,
      M = M,
      G = G
    )
  ))
}

#' Plot IRFs with stability checks
plot_factor_irfs <- function(irf, n_factors, horizon) {
  # Check for explosive behavior
  max_response <- max(abs(irf))
  if (max_response > 10) {
    warning("Large impulse responses detected - check VAR stability")
  }

  par(
    mfrow = c(ceiling(n_factors / 2), 2),
    mar = c(4, 4, 2, 1),
    oma = c(0, 0, 2, 0)
  )

  periods <- 0:(horizon - 1)

  for (i in 1:n_factors) {
    response <- irf[i, 1, ]

    ylim <- c(min(response), max(response))
    # Ajuste os limites se forem muito grandes
    if (diff(ylim) > 5) {
      ylim <- c(-2, 2)
    }

    plot(periods, response,
      type = "l",
      xlab = "Horizonte",
      ylab = "Resposta",
      main = paste("Fator", i),
      ylim = ylim,
      lwd = 2
    )

    abline(h = 0, lty = 2, col = "gray")
  }

  mtext("Respostas dos Fatores ao Primeiro Choque",
    outer = TRUE,
    cex = 1.2
  )
}

# Parâmetros do VAR
var_par <- list(
  nlag = 2,
  icomp = 1,
  iconst = 0
)

# Estimar VAR e obter matrizes
var_out <- varest_r(factors$factors, var_par)

irfs <- varirf_ss(var_out$coef, 0, 50)

plot_factor_irfs(irfs, n_factors = 6, horizon = 50)













# Estimação dos fatores

factors <- estimate_factor(
  X = data, n_factors = 6,
  constraint_info = named_factor
)


factors$loadings



# Função para ajustar o modelo VAR
ajustar_var <- function(fatores, p = 1) {
  # Ajuste o modelo VAR nos fatores estimados, sempre com type = "none"
  var_model <- vars::VAR(fatores, p = p, type = "none")

  # Obtenha a matriz de covariância dos resíduos do VAR
  residual_cov <- residuals(var_model)
  residual_cov <- cov(residual_cov)

  # Realize a decomposição de Cholesky da matriz de covariância dos resíduos
  cholesky_decomp <- t(chol(residual_cov))

  # Aplique a normalização de efeito unitário
  H <- solve(diag(diag(cholesky_decomp))) %*% cholesky_decomp

  return(list(var_model = var_model, H = H)) # Retorna o modelo ajustado e a matriz H
}

plotar_irfs <- function(var_model, H, impulso, resposta = NULL, horizonte = 15, boot = FALSE, runs = 100) {
  irfs <- vars::irf(var_model,
    impulse = impulso,
    response = resposta,
    n.ahead = horizonte,
    ortho = FALSE,
    B = H,
    cumulative = FALSE,
    boot = boot,
    runs = if (boot) runs else NULL
  )

  plot(irfs)
}

# Exemplo de uso
ajuste <- ajustar_var(fatores = factors$factors, p = 1)
plotar_irfs(
  var_model = ajuste$var_model,
  H = ajuste$H,
  impulso = "y1",
  resposta = "y2",
  horizonte = 15,
  boot = TRUE,
  runs = 100
)

for (i in seq_len(6)) {
  file_name <- paste0("img/grafico_y", i, ".png")
  png(filename = file_name, width = 680, height = 420)
  plotar_irfs(
    var_model = ajuste$var_model,
    H = ajuste$H,
    impulso = "y1",
    resposta = paste0("y", i),
    horizonte = 15,
    boot = TRUE,
    runs = 200
  )
  dev.off() # Fecha o dispositivo gráfico e salva o arquivo
  print(paste("Gráfico salvo como:", file_name))
  Sys.sleep(5) # Pausa de 5 segundos
}











# Teste SVAR() ----


# Função para ajustar o modelo SVAR
ajustar_svar <- function(fatores, p = 1) {
  # Ajuste o modelo VAR nos fatores estimados, sempre com type = "none"
  var_model <- vars::VAR(fatores, p = p, type = "none")

  # Defina a matriz A com todos os elementos como NA
  A <- matrix(NA, nrow = ncol(fatores), ncol = ncol(fatores))
  diag(A) <- 1

  # Ajuste o modelo SVAR usando a matriz A
  svar_model <- vars::SVAR(var_model, Amat = A, Bmat = NULL, max.iter = 100)

  return(svar_model)
}

plotar_irfs <- function(svar_model, impulso, resposta = NULL, horizonte = 15, boot = FALSE, runs = 100) {
  irfs <- vars::irf(svar_model,
    impulse = impulso,
    response = resposta,
    n.ahead = horizonte,
    boot = boot,
    runs = if (boot) runs else NULL
  )
  plot(irfs)
}

# Exemplo de uso
svar_ajuste <- ajustar_svar(fatores = factors$factors, p = 1)
plotar_irfs(
  svar_model = svar_ajuste,
  impulso = "y1",
  resposta = "y3",
  horizonte = 15,
  boot = TRUE,
  runs = 100
)
