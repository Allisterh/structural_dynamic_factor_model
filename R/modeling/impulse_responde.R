compute_irf_dfm <- function(dfm_results, h = 24, nboot = 300) {
  # Extract components
  Lambda <- dfm_results$static_loadings
  A <- dfm_results$companion_matrix
  K <- dfm_results$dynamic_loadings
  M <- dfm_results$dynamic_scaling
  sy <- dfm_results$data_sd

  # Dimensions
  n_vars <- nrow(Lambda)
  r <- ncol(Lambda)
  q <- ncol(K)

  # Normalize standard deviations and shock size as in MATLAB
  shock_size <- 5 # 50 basis points
  sy_normalized <- sy * 100 # Convert to percentage changes

  # Point IRF
  irf_point <- array(0, dim = c(n_vars, h + 1, q))
  B <- array(0, dim = c(r, r, h + 1))
  B[, , 1] <- diag(r)
  B[, , 2] <- A[1:r, 1:r]

  for (i in 3:(h + 1)) {
    B[, , i] <- A[1:r, 1:r] %*% B[, , i - 1]
  }

  # Calculate point IRF with correct normalization
  for (s in 1:q) {
    for (i in 1:(h + 1)) {
      temp <- Re((Lambda %*% B[, , i] %*% K %*% M))
      # Scale to get percentage changes
      irf_point[, i, s] <- temp[, s] * sy_normalized * shock_size
    }
  }

  # Bootstrap with same normalization
  irf_boot <- array(0, dim = c(n_vars, h + 1, nboot, q))
  rr <- matrix(1 - 2 * (runif(nrow(dfm_results$var_residuals) * nboot) > 0.5),
    nrow = nrow(dfm_results$var_residuals),
    ncol = nboot
  )

  for (b in 1:nboot) {
    resid_boot <- dfm_results$var_residuals * rr[, b]

    F_boot <- matrix(0, nrow = nrow(dfm_results$static_factors), ncol = r)
    F_boot[1, ] <- dfm_results$static_factors[1, ]

    for (t in 2:nrow(F_boot)) {
      F_boot[t, ] <- A[1:r, 1:r] %*% F_boot[t - 1, ] + resid_boot[t - 1, ]
    }

    var_boot <- estimate_corrected_var(F_boot, 1)
    A_boot <- var_boot$companion

    B_boot <- array(0, dim = c(r, r, h + 1))
    B_boot[, , 1] <- diag(r)
    B_boot[, , 2] <- A_boot[1:r, 1:r]

    for (i in 3:(h + 1)) {
      B_boot[, , i] <- A_boot[1:r, 1:r] %*% B_boot[, , i - 1]
    }

    for (s in 1:q) {
      for (i in 1:(h + 1)) {
        temp <- Re((Lambda %*% B_boot[, , i] %*% K %*% M))
        irf_boot[, i, b, s] <- temp[, s] * sy_normalized * shock_size
      }
    }
  }

  # Confidence intervals
  irf_ci <- array(0, dim = c(n_vars, h + 1, 5, q))
  for (s in 1:q) {
    irf_ci[, , 1, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.05)
    irf_ci[, , 2, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.10)
    irf_ci[, , 3, s] <- irf_point[, , s]
    irf_ci[, , 4, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.90)
    irf_ci[, , 5, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.95)
  }

  return(irf_ci)
}

plot_irf <- function(irf_results, response_vars, shock = 3, horizon = 20) {
  # Criar lista para armazenar os plots individuais
  plot_list <- list()

  # Loop através das variáveis de resposta
  for (i in seq_along(response_vars)) {
    var_index <- as.numeric(response_vars[[i]])
    var_name <- names(response_vars[[i]])

    df_plot <- data.frame(
      tempo = seq_len(horizon + 1) - 1,
      irf = irf_results[var_index, seq_len(horizon + 1), 3, shock],
      ic_05 = irf_results[var_index, seq_len(horizon + 1), 1, shock],
      ic_10 = irf_results[var_index, seq_len(horizon + 1), 2, shock],
      ic_90 = irf_results[var_index, seq_len(horizon + 1), 4, shock],
      ic_95 = irf_results[var_index, seq_len(horizon + 1), 5, shock]
    )

    plot_list[[i]] <- ggplot2::ggplot(df_plot, ggplot2::aes(x = tempo)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ic_05, ymax = ic_95),
        fill = "grey90", alpha = 0.5
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ic_10, ymax = ic_90),
        fill = "grey70", alpha = 0.5
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = irf),
        color = "black",
        linewidth = 1
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::theme_classic() +
      ggplot2::labs(y = var_name) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.4)),
        axis.text = ggplot2::element_text(size = ggplot2::rel(1.2)),
        plot.title = ggplot2::element_blank()
      ) +
      ggplot2::annotate(
        "text",
        x = Inf,
        y = mean(range(df_plot$ic_95, df_plot$ic_05)),
        label = var_name,
        hjust = 0,
        vjust = 0.5
      )
  }

  # Combinar os plots em um grid vertical
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)

  return(combined_plot)
}
