compute_SIRF <- function(Lambda, var_model, G = NULL, H, horizon = 40) {
  A_coef <- vars::Acoef(var_model)
  r <- nrow(A_coef[[1]])
  p <- length(A_coef)
  n <- nrow(Lambda)

  C <- array(0, dim = c(r, r, horizon + 1))
  C[, , 1] <- diag(r)

  # Recursão correta dos coeficientes MA
  for (h in 2:(horizon + 1)) {
    temp <- matrix(0, r, r)
    for (j in 1:min(h - 1, p)) {
      temp <- temp + C[, , h - j] %*% A_coef[[j]]
    }
    C[, , h] <- temp
  }

  # SIRF = Lambda * C * H
  SIRF <- array(0, dim = c(n, ncol(H), horizon + 1))
  if (is.null(G)) {
    for (h in 1:(horizon + 1)) {
      SIRF[, , h] <- Lambda %*% C[, , h] %*% H
    }
  }

  return(SIRF)
}
