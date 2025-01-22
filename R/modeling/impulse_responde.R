#' Compute Structural Impulse Response Functions (SIRF)
#'
#' @description
#' Calculates structural impulse response functions for a VAR model
#' with optional factor transformation.
#'
#' @param Lambda Factor loading matrix
#' @param var_model VAR model object
#' @param G Optional transformation matrix (default: NULL)
#' @param H Orthogonalization matrix
#' @param horizon Forecast horizon (default: 40)
#'
#' @return A 3D array of structural impulse response functions
#'
#' @details
#' Computes impulse response functions using VAR model coefficients
#' and orthogonalization matrices. Handles both standard and
#' reduced-rank factor models.
#'
#' @export
compute_SIRF <- function(Lambda, var_model, G = NULL, H, horizon = 40) {
  A_coef <- vars::Acoef(var_model)
  r <- nrow(A_coef[[1]])
  p <- length(A_coef)
  n <- nrow(Lambda)

  # Inicializar SIRF e C
  SIRF <- array(0, dim = c(n, ncol(H), horizon + 1))
  C <- array(0, dim = c(r, r, horizon + 1))
  C[, , 1] <- diag(r)

  # Recursão corrigida para coeficientes MA
  for (h in 1:horizon) {
    temp <- matrix(0, r, r)
    for (j in 1:min(h, p)) {
      temp <- temp + A_coef[[j]] %*% C[, , h - j + 1]
    }
    C[, , h + 1] <- temp
  }

  # Computar SIRF
  if (is.null(G)) {
    for (h in 1:(horizon + 1)) {
      SIRF[, , h] <- Lambda %*% C[, , h] %*% H
    }
  } else {
    for (h in 1:(horizon + 1)) {
      SIRF[, , h] <- Lambda %*% C[, , h] %*% G %*% H
    }
  }

  return(SIRF)
}



#' Plot Structural Impulse Response Functions (SIRF)
#'
#' @description
#' Generates plots of structural impulse response functions for specified
#' variables and shocks.
#'
#' @param sirf 3D array of structural impulse response functions
#' @param variables Indices of variables to plot (default: all variables)
#' @param shocks Indices of shocks to plot (default: all shocks)
#' @param horizon Maximum time horizon to plot (default: full horizon)
#'
#' @return A list of ggplot2 plot objects, one for each shock
#'
#' @details
#' Creates line plots showing the response of each variable to specific
#' structural shocks across different time horizons.
#'
#' @export
plot_sirf <- function(
    sirf, variables = NULL, shocks = NULL,
    horizon = NULL) {
  # Get dimensions
  n_var <- dim(sirf)[1]
  n_shocks <- dim(sirf)[2]
  h <- dim(sirf)[3]

  # Set defaults
  if (is.null(variables)) variables <- 1:n_var
  if (is.null(shocks)) shocks <- 1:n_shocks
  if (is.null(horizon)) horizon <- h - 1

  time <- 0:horizon

  # Create long format data
  irf_df <- tidyr::crossing(
    variable = variables,
    shock = shocks,
    horizon = time
  ) |>
    dplyr::mutate(
      irf = purrr::map_dbl(
        1:dplyr::n(),
        ~ sirf[variable[.x], shock[.x], horizon[.x] + 1]
      )
    )

  # Create plots
  plots <- list()
  for (s in shocks) {
    p <- ggplot2::ggplot(
      dplyr::filter(irf_df, shock == s),
      ggplot2::aes(x = horizon, y = irf)
    ) +
      ggplot2::geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "gray50"
      ) +
      ggplot2::geom_line(linewidth = 1, color = "steelblue") +
      ggplot2::facet_wrap(~variable, scales = "free_y") +
      ggplot2::labs(
        x = "Horizon",
        y = "Response",
        title = paste("Responses to Shock", s)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        strip.background = ggplot2::element_rect(fill = "white")
      )

    plots[[s]] <- p
  }

  return(plots)
}
