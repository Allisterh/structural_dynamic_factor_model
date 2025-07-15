compute_irf_dfm <- function(dfm_results, h = 24, nboot = 300) {
  # ===================================================================
  # CÁLCULO DE FUNÇÕES DE IMPULSO-RESPOSTA (IRF) PARA MODELO DFM
  # ===================================================================
  
  # Extract components
  Lambda <- dfm_results$static_loadings
  A <- dfm_results$companion_matrix
  K <- dfm_results$dynamic_loadings
  M <- dfm_results$dynamic_scaling
  sy <- dfm_results$data_sd
  p <- dfm_results$p  # Get VAR order from DFM results

  # Dimensions
  n_vars <- nrow(Lambda)
  r <- ncol(Lambda)
  q <- dfm_results$q  # Use stored q value directly

  # Point IRF calculation
  
  # Point IRF
  irf_point <- array(0, dim = c(n_vars, h + 1, q))
  B <- array(0, dim = c(r, r, h + 1))
  B[, , 1] <- diag(r)
  B[, , 2] <- A[1:r, 1:r]

  for (i in 3:(h + 1)) {
    B[, , i] <- A[1:r, 1:r] %*% B[, , i - 1]
  }

  # Calculate point IRF using formula
  # C(:,:,ii)=(lambda*B(:,:,ii)*K*M).*repmat(sy',1,q)
  for (i in 1:(h + 1)) {
    if (!is.matrix(K) && !is.matrix(M)) {
      # Case q=r with scalar K,M (K=1, M=1)
      # In this case, q factors but K=1, M=1 so we replicate across q
      temp <- Re(Lambda %*% B[, , i] * K * M)
      for (s in 1:q) {
        # Each dynamic factor gets the same response (since K=1, M=1)
        irf_point[, i, s] <- temp[, s] * sy  # Extract s-th static factor response
      }
    } else {
      # General case with matrix K,M
      temp <- Re(Lambda %*% B[, , i] %*% K %*% M)
      for (s in 1:q) {
        irf_point[, i, s] <- temp[, s] * sy  # Element-wise multiplication with sy
      }
    }
  }

  # Bootstrap with same normalization
  if (nboot > 0) {
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

      var_boot <- estimate_corrected_var(F_boot, p)
      A_boot <- var_boot$companion

      B_boot <- array(0, dim = c(r, r, h + 1))
      B_boot[, , 1] <- diag(r)
      B_boot[, , 2] <- A_boot[1:r, 1:r]

      for (i in 3:(h + 1)) {
        B_boot[, , i] <- A_boot[1:r, 1:r] %*% B_boot[, , i - 1]
      }

      for (i in 1:(h + 1)) {
        if (!is.matrix(K) && !is.matrix(M)) {
          # Case q=r with scalar K,M
          temp <- Re(Lambda %*% B_boot[, , i] * K * M)
          for (s in 1:q) {
            irf_boot[, i, b, s] <- temp[, s] * sy
          }
        } else {
          # General case with matrix K,M
          temp <- Re(Lambda %*% B_boot[, , i] %*% K %*% M)
          for (s in 1:q) {
            irf_boot[, i, b, s] <- temp[, s] * sy
          }
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
  } else {
    # No bootstrap case - return only point estimates
    irf_ci <- array(0, dim = c(n_vars, h + 1, 5, q))
    for (s in 1:q) {
      irf_ci[, , 1, s] <- irf_point[, , s]  # Same as point estimate
      irf_ci[, , 2, s] <- irf_point[, , s]  # Same as point estimate
      irf_ci[, , 3, s] <- irf_point[, , s]  # Point estimate
      irf_ci[, , 4, s] <- irf_point[, , s]  # Same as point estimate
      irf_ci[, , 5, s] <- irf_point[, , s]  # Same as point estimate
    }
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

# ===================================================================
# FUNÇÕES AUXILIARES PARA LIDAR COM CASOS q=r (K,M ESCALARES)
# ===================================================================

#' Helper function to get number of dynamic factors
#' Handles case when K is scalar (q=r) or matrix (q<r)
get_q_from_K <- function(K) {
  if (is.matrix(K)) {
    return(ncol(K))
  } else {
    return(1)  # K is scalar, so q=1
  }
}

#' Validate DFM results
#' @param dfm_results Results from estimate_dfm function
#' @return List of validation checks
validate_dfm_results <- function(dfm_results) {
  checks <- list()
  
  # Check if all required components exist
  required_components <- c("static_loadings", "companion_matrix", 
                          "dynamic_loadings", "dynamic_scaling", 
                          "data_sd", "p", "r", "q")
  
  checks$missing_components <- setdiff(required_components, names(dfm_results))
  
  # Check q=r case handling
  q <- dfm_results$q
  r <- dfm_results$r
  K <- dfm_results$dynamic_loadings
  M <- dfm_results$dynamic_scaling
  
  if (q == r) {
    checks$qr_case_K_scalar <- !is.matrix(K) && length(K) == 1
    checks$qr_case_M_scalar <- !is.matrix(M) && length(M) == 1
    checks$qr_case_K_equals_1 <- K == 1
    checks$qr_case_M_equals_1 <- M == 1
  } else {
    checks$qr_case_K_matrix <- is.matrix(K) && ncol(K) == q
    checks$qr_case_M_matrix <- is.matrix(M) && ncol(M) == q
  }
  
  # Check companion matrix stability
  eigenvals <- eigen(dfm_results$companion_matrix)$values
  checks$companion_stable <- all(abs(eigenvals) < 1)
  checks$max_eigenvalue <- max(abs(eigenvals))
  
  return(checks)
}
