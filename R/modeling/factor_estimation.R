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
amengual_watson <- function(X, r, p = 4, max_q = NULL) {
  X <- as.matrix(X)

  if (is.null(max_q)) max_q <- r

  pca_X <- prcomp(X, scale. = FALSE)
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

  bn <- bai_ng_criteria(resid_mat, max_r = max_q, standardize = TRUE)
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




#' Create a Constraint Matrix from Group List
#'
#' @description
#' Generates a constraint matrix where each group is represented by a column
#' of ones in a sparse matrix, useful for various statistical and machine
#' learning applications.
#'
#' @param list_groups A list of groups, where each element contains the
#' members of a specific group.
#' @param r An integer representing the number of groups (default is 6).
#' This determines the number of columns in the resulting matrix.
#'
#' @return A matrix with binary values, where each group has a column
#' of ones corresponding to its members, and zeros elsewhere.
#'
#' @details
#' The function creates a constraint matrix by:
#' 1. Initializing a matrix of zeros for each group
#' 2. Setting the corresponding group column to ones
#' 3. Combining these matrices row-wise
#'
#' @examples
#' # Create groups
#' groups <- list(
#'   group1 = c("A", "B", "C"),
#'   group2 = c("D", "E"),
#'   group3 = c("F", "G", "H")
#' )
#'
#' # Generate constraint matrix
#' result <- constraint_matrix(groups)
#' print(result)
#'
#' @export
constraint_matrix <- function(list_groups, r = 6) {
  # Initialize an empty list to store matrices
  constraint_list <- list()

  # Loop through each group in the list
  for (i in seq_along(list_groups)) {
    # Create a matrix of zeros
    temp <- matrix(0, nrow = length(list_groups[[i]]), ncol = r)
    # Set the i-th column to 1
    temp[, i] <- 1
    # Add the matrix to the constraint_list
    constraint_list[[i]] <- temp
  }

  # Combine all matrices row-wise into a single matrix
  constraint_matrix <- do.call(rbind, constraint_list)

  return(constraint_matrix)
}






#' Estimate Factors with Constrained Loadings
#'
#' @description
#' Performs factor estimation using a constrained approach, allowing
#' specification of group-based loading constraints.
#'
#' @param X A matrix or data frame of variables to be factor analyzed
#'
#' @param n_factors Number of factors to extract
#'
#' @param constraint_info Matrix of predefined constraints for factor loadings
#'
#' @param group_vars List of variable groups with predefined constraints
#'
#' @param tol Convergence tolerance (default: 1e-8)
#'
#' @param max_iter Maximum number of iterations (default: 1000)
#'
#' @return A list containing:
#' - factors: Estimated factor scores
#' - loadings: Factor loadings matrix
#' - R2: Overall model R-squared
#' - R2_series: R-squared for individual variables
#' - iterations: Number of iterations performed
#' - objective: Final objective function value
#' - residuals: Residual matrix
#'
#' @details
#' The function uses an iterative algorithm to estimate factors with
#' constrained loadings. It begins with standard PCA initialization
#' and updates factor scores and loadings while respecting specified
#' group constraints.
#'
#' @examples
#' # Assuming X is your data matrix and you have predefined constraints
#' result <- estimate_factor(
#'   X = data_matrix,
#'   n_factors = 3,
#'   constraint_info = constraint_matrix,
#'   group_vars = list_of_groups
#' )
#'
#' @export
estimate_factor <- function(
    X, n_factors,
    constraint_info,
    group_vars,
    tol = 1e-8,
    max_iter = 1000) {
  X <- X |>
    as.matrix() |>
    scale()

  # Initialize dimensions
  T <- nrow(X)
  N <- ncol(X)

  # Get indices for constrained variables
  constrained_indices <- unlist(lapply(group_vars, function(x) {
    match(x, colnames(X))
  }))
  unconstrained_indices <- setdiff(1:N, constrained_indices)

  # Initialize factors using standard PCA
  svd_result <- svd(X)
  F_hat <- sqrt(T) * svd_result$u[, 1:n_factors]
  Lambda_hat <- matrix(0, N, n_factors)

  # Set initial constrained loadings
  Lambda_hat[constrained_indices, ] <- constraint_info

  # Initialize objective function value
  V_old <- sum((X - F_hat %*% t(Lambda_hat))^2) / (N * T)

  # Iterate until convergence
  for (iter in 1:max_iter) {
    # Update unconstrained lambdas all at once
    Lambda_hat[unconstrained_indices, ] <- t(X[, unconstrained_indices]) %*% F_hat / T

    # Update factors all at once using matrix operations
    F_hat <- X %*% Lambda_hat %*% solve(t(Lambda_hat) %*% Lambda_hat)

    # Calculate new objective function value
    V_new <- sum((X - F_hat %*% t(Lambda_hat))^2) / (N * T)

    # Check convergence
    if (abs(V_new - V_old) < tol) {
      break
    }

    V_old <- V_new
  }

  # Calculate fit statistics
  res <- X - F_hat %*% t(Lambda_hat)
  TSS <- sum(X^2)
  RSS <- sum(res^2)
  R2 <- 1 - RSS / TSS
  R2_series <- 1 - colSums(res^2) / colSums(X^2)

  return(list(
    factors = F_hat,
    loadings = Lambda_hat,
    R2 = R2,
    R2_series = R2_series,
    iterations = iter,
    objective = V_new,
    residuals = res
  ))
}




#' Estimate G matrix when r > q
#' @param F_t Matrix of static factors
#' @param q Number of dynamic factors
#' @param p_var VAR lag order
#' @return List containing G matrix and factor innovations
estimate_G <- function(F_t, q, p_var = 4) {
  r <- ncol(F_t)

  # Estimate VAR on factors to get innovations
  var_model <- vars::VAR(F_t, p = p_var)
  factor_innovations <- residuals(var_model)

  # Compute sample covariance of innovations
  Sigma_a <- cov(factor_innovations)

  # Partition Sigma_a
  Sigma_a11 <- Sigma_a[1:q, 1:q]
  Sigma_a21 <- Sigma_a[(q + 1):r, 1:q]

  # Construct G
  G <- matrix(0, r, q)
  G[1:q, ] <- diag(q) # Upper block is identity matrix
  G[(q + 1):r, ] <- Sigma_a21 %*% solve(Sigma_a11) # Lower block from regression

  # Get eta_t (q-dimensional innovations)
  eta_t <- factor_innovations[, 1:q]

  return(list(
    G = G,
    eta = eta_t,
    var_model = var_model
  ))
}



#' Estimate Structural Dynamic Factor Model (SDFM)
#'
#' @description
#' Estimates a structural dynamic factor model with constrained loadings
#' and dynamic factor extraction.
#'
#' @param X Input data matrix
#' @param r Number of static factors (default: 8)
#' @param q Number of dynamic factors (default: equal to r)
#' @param p_var Lag order for VAR model (default: 4)
#' @param constraint_matrix Matrix of loading constraints
#' @param group_vars List of variable groups
#' @param horizon Forecast horizon for structural impulse response functions
#'
#' @return A list containing:
#' - Lambda: Factor loadings
#' - F: Factor scores
#' - G: Transformation matrix (if r > q)
#' - H: Orthogonalization matrix
#' - eta: Innovations
#' - Sigma_eta: Normalized covariance of innovations
#' - Sigma_e: Error covariance
#' - var_model: VAR model
#' - sirf: Structural impulse response functions
#' - convergence: Model convergence information
#'
#' @details
#' Handles two scenarios:
#' 1. Standard VAR approach when r = q
#' 2. Reduced-rank dynamic factor model when r > q
#'
#' @export
estimate_sdfm <- function(
    X, r = 8, q = NULL, p_var = 4,
    constraint_matrix,
    group_vars,
    horizon) {
  if (is.null(q)) q <- r # Default case r=q

  if (q > r) stop("Number of dynamic factors cannot exceed number of static factors")

  # Estimate named factors
  named_factors <- estimate_factor(X,
    n_factors = r,
    constraint_info = constraint_matrix,
    group_vars = group_vars
  )


  if (r == q) {
    # Case r=q: Standard VAR approach
    var_model <- vars::VAR(named_factors$factors, p = p_var)
    resid_var <- residuals(var_model)
    Sigma_eta <- cov(resid_var)

    # Unit Effect Normalization via Cholesky
    H_chol <- t(chol(Sigma_eta))
    D <- diag(diag(H_chol))
    H <- H_chol %*% solve(D)
    Sigma_eta_normalized <- H %*% t(H)

    G <- NULL
    eta <- resid_var
  } else {
    # Case r>q: Estimate G and get q-dimensional innovations
    G_estimates <- estimate_G(named_factors$factors, q, p_var)
    G <- G_estimates$G
    eta <- G_estimates$eta
    var_model <- G_estimates$var_model

    # Compute Sigma_eta and H for q-dimensional innovations
    Sigma_eta <- cov(eta)
    H_chol <- t(chol(Sigma_eta))
    D <- diag(diag(H_chol))
    H <- H_chol %*% solve(D)
    Sigma_eta_normalized <- H %*% t(H)
  }

  # Compute SIRFs
  sirf <- compute_SIRF(named_factors$loadings, var_model, G, H, horizon = horizon)

  return(list(
    Lambda = named_factors$loadings,
    F = named_factors$factors,
    G = G,
    H = H,
    eta = eta,
    Sigma_eta = Sigma_eta_normalized,
    Sigma_e = named_factors$Sigma_e,
    var_model = var_model,
    sirf = sirf,
    convergence = named_factors$convergence
  ))
}
