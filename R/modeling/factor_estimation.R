#' Calculate Bai and Ng (2002) Information Criteria for Factor Model Selection
#'
#' @description
#' Implements the Bai and Ng (2002) information criteria for determining the optimal
#' number of factors in approximate factor models. The function computes three different
#' penalty versions (IC1, IC2, and IC3) of the criteria.
#'
#' @param X A numeric matrix of dimensions T × N, where T is the number of time periods
#'          and N is the number of variables/series
#' @param max_r An integer specifying the maximum number of factors to consider.
#'              Default is 15
#' @param standardize Logical indicating whether to standardize the data before
#'                   computing the criteria. Default is TRUE
#'
#' @return A list containing three components:
#' \itemize{
#'   \item r_hat: A list with the optimal number of factors according to each criterion
#'                (IC1, IC2, and IC3)
#'   \item criteria: A list containing the values of each information criterion for
#'                  r = 1 to max_r
#'   \item pca: The prcomp object from the principal components analysis
#' }
#'
#' @details
#' The function implements the three panel information criteria proposed by Bai and Ng
#' (2002) for determining the number of factors in approximate factor models. The
#' criteria differ in their penalty terms:
#' \itemize{
#'   \item IC1 uses penalty term r * (N+T)/(NT) * log((NT)/(N+T))
#'   \item IC2 uses penalty term r * (N+T)/(NT) * log(min(N,T))
#'   \item IC3 uses penalty term r * log(min(N,T))/min(N,T)
#' }
#'
#' @references
#' Bai, J., & Ng, S. (2002). Determining the Number of Factors in Approximate Factor
#' Models. Econometrica, 70(1), 191-221.
#'
#' @examples
#' # Generate random data
#' set.seed(123)
#' T <- 100 # time periods
#' N <- 20 # variables
#' X <- matrix(rnorm(T * N), nrow = T)
#'
#' # Calculate Bai-Ng criteria
#' results <- bai_ng_criteria(X, max_r = 10)
#'
#' # View optimal number of factors for each criterion
#' print(results$r_hat)
#'
#' @export
bai_ng_criteria <- function(X, max_r = 15, standardize = TRUE) {
  T <- nrow(X)
  N <- ncol(X)

  if (standardize) {
    X_std <- scale(X, center = TRUE, scale = TRUE)
  } else {
    X_std <- X
  }

  pca <- prcomp(X_std)

  ic1 <- ic2 <- ic3 <- numeric(max_r)

  for (r in 1:max_r) {
    # Extract first r factors and loadings
    F_hat <- pca$x[, 1:r]
    Lambda_hat <- pca$rotation[, 1:r]

    # Compute residuals from r-factor approximation
    resid <- X_std - F_hat %*% t(Lambda_hat)
    V_r <- sum(resid^2) / (N * T)

    # Panel dimension adjustment term
    e <- (N + T) / (N * T)

    # Penalty terms for each IC
    p1 <- r * e * log(1 / e)
    p2 <- r * e * log(min(N, T))
    p3 <- r * log(min(N, T)) / min(N, T)

    ic1[r] <- log(V_r) + p1
    ic2[r] <- log(V_r) + p2
    ic3[r] <- log(V_r) + p3
  }

  return(list(
    r_hat = list(
      IC1 = which.min(ic1),
      IC2 = which.min(ic2),
      IC3 = which.min(ic3)
    ),
    criteria = list(IC1 = ic1, IC2 = ic2, IC3 = ic3),
    pca = pca
  ))
}





#' Test for Dynamic Factors Using Amengual and Watson (2007) Method
#'
#' @description
#' Implements the Amengual and Watson (2007) method for determining the number of
#' dynamic factors in a dynamic factor model. The procedure first estimates static
#' factors via PCA, then applies VAR filtering, and finally uses Bai and Ng (2002)
#' criteria on the residuals to determine the number of dynamic factors.
#'
#' @param X A numeric matrix of dimensions T × N, where T is the number of time periods
#'          and N is the number of variables/series
#' @param r An integer specifying the number of static factors to extract in the first
#'          step
#' @param p An integer specifying the number of lags to include in the VAR filtering.
#'          Default is 4
#' @param max_q An integer specifying the maximum number of dynamic factors to consider.
#'              If NULL (default), it is set equal to r
#'
#' @return A list containing two components:
#' \itemize{
#'   \item aw: A numeric vector containing the Amengual-Watson test statistics (based
#'             on IC2 criterion) for q = 1 to max_q
#'   \item q_hat: An integer indicating the estimated number of dynamic factors based
#'                on the IC2 criterion
#' }
#'
#' @details
#' The function implements the following steps:
#' \itemize{
#'   \item Extracts r static factors using principal components
#'   \item Applies VAR(p) filtering to remove dynamic dependence
#'   \item Uses Bai-Ng IC2 criterion on the residuals to determine q
#' }
#'
#' The procedure allows for the number of dynamic factors (q) to be smaller than
#' the number of static factors (r), which is often the case in empirical applications.
#'
#' @references
#' Amengual, D., & Watson, M. W. (2007). Consistent Estimation of the Number of
#' Dynamic Factors in a Large N and T Panel. Journal of Business & Economic
#' Statistics, 25(1), 91-96.
#'
#' @examples
#' # Generate random data
#' set.seed(123)
#' T <- 100 # time periods
#' N <- 20 # variables
#' X <- matrix(rnorm(T * N), nrow = T)
#'
#' # Test for dynamic factors with 3 static factors
#' results <- amengual_watson(X, r = 3, p = 2)
#'
#' # View estimated number of dynamic factors
#' print(results$q_hat)
#'
#' @export
amengual_watson <- function(X, r, p = 4, max_q = NULL, scale = TRUE) {
  X <- as.matrix(X)

  if (scale) {
    X <- scale(X, center = TRUE, scale = TRUE)
  }

  if (is.null(max_q)) max_q <- r

  pca_X <- prcomp(X, scale. = FALSE, center = FALSE)
  F_hat <- pca_X$x[, 1:r, drop = FALSE]

  T <- nrow(X)
  lagged_factors <- stats::embed(F_hat, p + 1)
  Z <- cbind(1, lagged_factors[, -(1:r), drop = FALSE])

  resid_mat <- matrix(NA, T - p, ncol(X))

  for (i in 1:ncol(X)) {
    y <- X[(p + 1):T, i]
    z <- Z[1:length(y), , drop = FALSE]

    b <- solve(crossprod(z), crossprod(z, y))
    e <- y - z %*% b
    resid_mat[, i] <- e
  }

  bn <- bai_ng_criteria(resid_mat, max_r = max_q, standardize = FALSE)
  q_hat <- bn$r_hat$IC2

  list(
    aw = bn$criteria$IC2,
    q_hat = q_hat
  )
}


#' Generate scree plot analysis following Stock & Watson (2016)
#' @param X Matrix of standardized data
#' @param max_comp Maximum number of components to plot (default = 15)
#' @return List containing variance decomposition and plot
scree_analysis <- function(X, max_comp = 15) {
  # Get dimensions
  T <- nrow(X)
  N <- ncol(X)
  max_comp <- min(max_comp, N) # Ensure max_comp doesn't exceed N

  # Standardize data
  X_std <- scale(X)

  # Get PCA
  pca <- stats::prcomp(X_std, scale. = FALSE) # Already standardized

  # Compute R2 for each number of factors
  r2_vec <- numeric(max_comp)
  for (k in 1:max_comp) {
    F_hat <- pca$x[, 1:k, drop = FALSE]
    ssr <- 0

    # For each variable, compute SSR with k factors
    for (i in 1:N) {
      y <- X_std[, i]
      lambda <- solve(crossprod(F_hat), crossprod(F_hat, y))
      resid <- y - F_hat %*% lambda
      ssr <- ssr + sum(resid^2)
    }

    # Compute R2
    r2_vec[k] <- 1 - ssr / sum(X_std^2)
  }

  # Compute marginal R2
  marg_r2 <- c(r2_vec[1], diff(r2_vec))

  # Create data frame for ggplot
  plot_data <- data.frame(
    k = 1:max_comp,
    marginal_r2 = marg_r2
  )

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = k, y = marginal_r2)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 4) +
    ggplot2::labs(
      x = "Number of Components",
      y = "Proportion of Variance Explained",
      title = "Scree Plot: Marginal Contribution of Each Factor"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  # Return results
  results <- list(
    r2 = r2_vec,
    marginal_r2 = marg_r2,
    eigenvalues = pca$sdev[1:max_comp]^2,
    cumulative_r2 = cumsum(marg_r2),
    plot = p
  )

  return(invisible(results))
}




estimate_static_factors <- function(data, r) {
  # Dimensões
  T <- nrow(data)
  N <- ncol(data)
  
  # ===================================================================
  # PADRONIZAÇÃO BLL (Barigozzi, Lippi, Luciani 2016)
  # Implementação baseada no código MATLAB original DFMest_BLL.m
  # ===================================================================
  
  # 1. Calcular primeira diferença para obter desvio padrão
  y <- diff(data)  # Primeira diferença: Y(t) - Y(t-1)
  sy <- apply(y, 2, sd)  # Desvio padrão das primeiras diferenças
  
  # 2. Padronizar as primeiras diferenças (para cálculo dos autovalores)
  y_centered <- sweep(y, 2, colMeans(y), "-")  # Centrar
  yy <- sweep(y_centered, 2, sy, "/")  # Padronizar: (y - mean(y))/sd(y)
  
  # 3. Detrending dos dados em nível
  # Construir matriz de regressores: [1, t] para remoção de tendência linear
  regX <- cbind(1, 1:T)  # Constante e tendência linear
  
  # Remover tendência de cada série individualmente
  beta <- solve(crossprod(regX)) %*% crossprod(regX, data)  # Coeficientes da regressão
  X <- data - regX %*% beta  # Dados destrendizados
  
  # 4. Padronização BLL final: dados destrendizados / sd das diferenças
  Z <- sweep(X, 2, sy, "/")  # Padronização BLL: X / sy
  
  # ===================================================================
  # EXTRAÇÃO DOS FATORES ESTÁTICOS
  # ===================================================================
  
  # 5. Calcular autovalores e autovetores da matriz de covariância das diferenças padronizadas
  cov_yy <- cov(yy)
  eigen_result <- eigen(cov_yy, symmetric = TRUE)
  
  # 6. Selecionar os r maiores autovalores (em ordem decrescente)
  idx <- order(eigen_result$values, decreasing = TRUE)[1:r]
  lambda <- eigen_result$vectors[, idx]  # Autovetores (loadings)
  
  # 7. Calcular fatores estáticos: F = Z * lambda
  F <- Z %*% lambda
  
  return(list(
    factors = F,                    # Fatores estáticos F = Z * lambda
    loadings = lambda,              # Loadings (autovetores)
    sy = sy,                        # Desvios padrão das diferenças
    Z = Z,                          # Dados padronizados BLL
    eigenvalues = eigen_result$values[idx],  # Autovalores selecionados
    yy = yy,                        # Diferenças padronizadas (para diagnóstico)
    detrended_data = X              # Dados destrendizados (para diagnóstico)
  ))
}



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


kilian_correction <- function(A, SIGMA, t, q, p) {
  # ===================================================================
  # CORREÇÃO DE VIÉS DE KILIAN (1998)
  # Implementação baseada no código MATLAB original kiliancorr.m
  # Fonte: Pope (1990), JTSA; Kilian (1998)
  # ===================================================================
  
  T <- t - p  # Tamanho efetivo da amostra
  I <- diag(q * p)  # Matriz identidade
  B <- t(A)  # Transposta da matriz companion
  
  # Calculate SIGMAY using Lyapunov equation
  A_kron_A <- kronecker(A, A)
  I_full <- diag((q * p)^2)
  
  lyapunov_matrix <- I_full - A_kron_A
  if (det(lyapunov_matrix) == 0) {
    lyapunov_inv <- MASS::ginv(lyapunov_matrix)
  } else {
    lyapunov_inv <- solve(lyapunov_matrix)
  }
  
  vecSIGMAY <- lyapunov_inv %*% c(SIGMA)
  SIGMAY <- matrix(vecSIGMAY, nrow = q * p, ncol = q * p)
  
  # Calculate eigenvalues of companion matrix
  peigen <- eigen(A)$values
  
  # Calculate bias correction terms
  sumeig <- matrix(0, q * p, q * p)
  for (h in 1:(q * p)) {
    if (abs(1 - peigen[h]) > 1e-10) {
      term <- peigen[h] * solve(I - peigen[h] * B)
      sumeig <- sumeig + term
    }
  }
  
  # Calculate analytical bias
  I_minus_B <- I - B
  I_minus_B2 <- I - B %*% B
  
  if (det(I_minus_B) == 0) {
    inv_I_minus_B <- MASS::ginv(I_minus_B)
  } else {
    inv_I_minus_B <- solve(I_minus_B)
  }
  
  if (det(I_minus_B2) == 0) {
    warning("Matriz (I-B^2) é singular. Usando pseudo-inversa.")
    inv_I_minus_B2 <- MASS::ginv(I_minus_B2)
  } else {
    inv_I_minus_B2 <- solve(I_minus_B2)
  }
  
  if (det(SIGMAY) == 0) {
    warning("Matriz SIGMAY é singular. Usando pseudo-inversa.")
    inv_SIGMAY <- MASS::ginv(SIGMAY)
  } else {
    inv_SIGMAY <- solve(SIGMAY)
  }
  
  # Calcular termo do viés
  bias_term <- inv_I_minus_B + B %*% inv_I_minus_B2 + sumeig
  bias <- SIGMA %*% bias_term %*% inv_SIGMAY
  
  # 5. Ajuste do viés (dividido por T conforme Kilian)
  Abias <- -bias / T
  
  # 6. Aplicar correção com fator de ajuste delta
  bcstab <- 9  # Valor arbitrário > 1
  delta <- 1   # Fator de ajuste
  
  cat("Aplicando correção iterativa...\n")
  iter <- 0
  max_iter <- 100
  
  while (bcstab >= 1 && iter < max_iter) {
    iter <- iter + 1
    
    # Ajustar correção proporcionalmente
    bcA <- A - delta * Abias
    
    # Verificar estabilidade (todos os autovalores < 1 em módulo)
    bcmod <- abs(eigen(bcA)$values)
    
    if (any(bcmod >= 1)) {
      bcstab <- 1
    } else {
      bcstab <- 0
      break
    }
    
    delta <- delta - 0.01
    
    if (delta <= 0) {
      bcstab <- 0
      bcA <- A
      break
    }
  }

  return(bcA)
}



estimate_corrected_var <- function(data, p) {
  T <- nrow(data)
  K <- ncol(data)
  
  # Construct regressor matrix
  RHS <- matrix(NA, T - p, K * p + 1)
  
  for (i in 1:p) {
    start_col <- (i - 1) * K + 1
    end_col <- i * K
    RHS[, start_col:end_col] <- data[(p + 1 - i):(T - i), ]
  }
  
  # Add constant
  RHS[, K * p + 1] <- 1
  RHS[, K * p + 1] <- 1
  
  # 2. Variável dependente
  LHS <- data[(p + 1):T, ]
  
  # 3. Estimação por OLS
  bet <- solve(crossprod(RHS)) %*% crossprod(RHS, LHS)
  
  # 4. Calcular resíduos
  u <- LHS - RHS %*% bet
  
  # 5. Construir matriz companion (excluindo constante)
  coeffcompanion <- rbind(
    t(bet[1:(p * K), ]),  # Coeficientes VAR
    cbind(diag((p - 1) * K), matrix(0, (p - 1) * K, K))  # Identidade para lags
  )
  
  # Calculate covariance matrix of residuals
  SIGMA <- cov(u)
  
  # Apply Kilian correction
  eigenvals_orig <- eigen(coeffcompanion)$values
  max_eigen_orig <- max(abs(eigenvals_orig))
  
  coeffcompanion_corrected <- kilian_correction(coeffcompanion, SIGMA, T, K, p)
  
  eigenvals_corr <- eigen(coeffcompanion_corrected)$values
  max_eigen_corr <- max(abs(eigenvals_corr))
  
  # Extract corrected coefficients
  beta_corrected <- t(coeffcompanion_corrected[1:K, ])
  
  # Recalculate residuals with corrected coefficients
  Y_resid <- data[(p + 1):T, ]
  X_resid <- matrix(0, T - p, K * p)
  
  for (i in 1:p) {
    X_resid[, ((i - 1) * K + 1):(i * K)] <- data[(p - i + 1):(T - i), ]
  }
  
  beta_corrected <- as.matrix(Re(beta_corrected))
  X_resid <- as.matrix(Re(X_resid))
  Y_resid <- as.matrix(Re(Y_resid))
  
  u_corrected <- Y_resid - X_resid %*% beta_corrected
  u_corrected <- as.matrix(Re(u_corrected))

  return(list(
    coefficients = beta_corrected,
    residuals = u_corrected,
    companion = coeffcompanion_corrected,
    residuals_original = u,
    coefficients_original = t(bet[1:(p * K), ]),
    companion_original = coeffcompanion,
    covariance_matrix = SIGMA
  ))
}



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
  # ===================================================================
  # ESTIMAÇÃO DOS FATORES DINÂMICOS
  # Baseado no código MATLAB DFMest_BLL.m (linhas 52-57)
  # ===================================================================
  
  cat("=== ESTIMAÇÃO DOS FATORES DINÂMICOS ===\n")
  cat("Número de fatores dinâmicos q:", q, "\n")
  cat("Número de fatores estáticos r:", r, "\n")
  
  if (q == r) {
    # Caso especial: q = r (fatores dinâmicos = fatores estáticos)
    cat("Caso especial: q = r. Fatores dinâmicos = fatores estáticos.\n")
    K <- diag(r)
    M <- diag(r)
    eta <- var_residuals
  } else {
    # Caso geral: q < r
    cat("Extraindo", q, "fatores dinâmicos de", r, "fatores estáticos.\n")
    
    # Calcular matriz de covariância dos resíduos VAR
    sigma_u <- cov(var_residuals)
    
    # Calcular autovalores e autovetores (usando os q maiores autovalores)
    eigen_result <- eigen(sigma_u, symmetric = TRUE)
    
    # Selecionar os q maiores autovalores
    idx <- order(eigen_result$values, decreasing = TRUE)[1:q]
    eigenvals <- eigen_result$values[idx]
    eigenvects <- eigen_result$vectors[, idx]
    
    # Construir matrizes K e M conforme MATLAB
    K <- eigenvects  # Autovetores
    M <- diag(sqrt(eigenvals))  # Raiz quadrada dos autovalores
    
    # Calcular fatores dinâmicos: eta = u * K / M
    eta <- var_residuals %*% K %*% solve(M)
    
    cat("Autovalores selecionados:", round(eigenvals, 4), "\n")
    cat("Proporção da variância explicada pelos fatores dinâmicos:", 
        round(sum(eigenvals) / sum(eigen_result$values) * 100, 2), "%\n")
  }
  
  cat("Dimensões dos fatores dinâmicos eta:", nrow(eta), "x", ncol(eta), "\n")
  cat("========================================\n")

  return(list(
    factors = eta,    # Fatores dinâmicos eta
    K = K,           # Matriz de autovetores
    M = M,           # Matriz diagonal com raiz dos autovalores
    eigenvalues = if(q == r) rep(1, r) else eigenvals  # Autovalores para diagnóstico
  ))
}



estimate_dfm <- function(data, r, q, p) {
  # ===================================================================
  # ESTIMAÇÃO COMPLETA DO MODELO DE FATORES DINÂMICOS ESTRUTURAIS (SDFM)
  # Implementação baseada em Alessi & Kerssenfischer com metodologia BLL
  # ===================================================================
  
  cat("\n")
  cat("=======================================================\n")
  cat("    ESTIMAÇÃO DO MODELO DE FATORES DINÂMICOS (SDFM)   \n")
  cat("=======================================================\n")
  cat("Parâmetros do modelo:\n")
  cat("- Número de fatores estáticos (r):", r, "\n")
  cat("- Número de fatores dinâmicos (q):", q, "\n")
  cat("- Ordem do VAR (p):", p, "\n")
  cat("- Dimensões dos dados:", nrow(data), "x", ncol(data), "\n")
  cat("=======================================================\n")

  # 1. ESTIMAÇÃO DOS FATORES ESTÁTICOS (com padronização BLL corrigida)
  cat("\nETAPA 1: ESTIMAÇÃO DOS FATORES ESTÁTICOS\n")
  static_result <- estimate_static_factors(data, r)

  static_result <- estimate_static_factors(data, r)

  var_result <- estimate_corrected_var(static_result$factors, p)

  dynamic_result <- estimate_dynamic_factors(var_result$residuals, q, r)

  # Diagnostics
  max_eigenval <- max(abs(eigen(var_result$companion)$values))
  is_stable <- max_eigenval < 1
  
  eta_corr <- cor(dynamic_result$factors)
  max_corr <- max(abs(eta_corr[upper.tri(eta_corr)]))
  
  total_var_static <- sum(static_result$eigenvalues) / sum(eigen(cov(static_result$yy))$values) * 100
  
  if (q < r) {
    total_var_dynamic <- sum(dynamic_result$eigenvalues) / sum(eigen(cov(var_result$residuals))$values) * 100
  } else {
    total_var_dynamic <- NA
  }

  return(list(
    # Static factors
    static_factors = static_result$factors,
    static_loadings = static_result$loadings,
    static_eigenvalues = static_result$eigenvalues,
    
    # VAR on static factors
    var_coefficients = var_result$coefficients,
    var_residuals = var_result$residuals,
    companion_matrix = var_result$companion,
    var_covariance = var_result$covariance_matrix,
    
    # Dynamic factors
    dynamic_factors = dynamic_result$factors,
    dynamic_loadings = dynamic_result$K,
    dynamic_scaling = dynamic_result$M,
    dynamic_eigenvalues = dynamic_result$eigenvalues,
    
    # Dados transformados e auxiliares
    data_sd = static_result$sy,
    Z = static_result$Z,
    yy = static_result$yy,
    detrended_data = static_result$detrended_data,
    
    # Componentes para IRF
    loadings_lambda = static_result$loadings,  # λ no código MATLAB
    scaling_matrix_K = dynamic_result$K,       # K no código MATLAB  
    scaling_matrix_M = dynamic_result$M,       # M no código MATLAB
    
    # Diagnósticos
    diagnostics = list(
      max_eigenvalue = max_eigenval,
      is_stable = max_eigenval < 1,
      max_dynamic_correlation = max_corr,
      static_variance_explained = total_var_static,
      dynamic_variance_explained = if(q < r) total_var_dynamic else 100
    )
  ))
}
