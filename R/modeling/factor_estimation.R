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

  # Print cumulative R2
  cat("Cumulative proportion of variance explained:\n")
  print(round(results$cumulative_r2, 4))

  # Display plot
  print(p)

  return(invisible(results))
}




estimate_static_factors <- function(data, r) {
  # Dimensões
  T <- nrow(data)

  # Primeira diferença apenas para calcular o desvio padrão
  y <- diff(data)
  sy <- apply(y, 2, sd)

  # Standardização BLL: dados em nível divididos pelo sd das diferenças
  Z <- sweep(data, 2, sy, "/")

  # PCA nos dados transformados
  pca_result <- prcomp(Z, scale = FALSE)

  return(list(
    factors = pca_result$x[, 1:r],
    loadings = pca_result$rotation[, 1:r],
    sy = sy,
    Z = Z
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

  return(Re(bcA))
}



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
    K <- diag(r)
    M <- diag(r)
    eta <- var_residuals
  } else {
    sigma_u <- cov(var_residuals)
    eig <- eigen(sigma_u)
    idx <- order(abs(eig$values), decreasing = TRUE)[1:q]
    K <- eig$vectors[, idx]
    M <- diag(sqrt(eig$values[idx]))
    eta <- var_residuals %*% K %*% solve(M)
  }

  return(list(
    factors = eta,
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
    dynamic_loadings = dynamic_result$K,
    dynamic_scaling = dynamic_result$M,
    data_sd = static_result$sy,
    Z = static_result$Z
  ))
}
