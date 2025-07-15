compute_irf_dfm <- function(dfm_results, h = 24, nboot = 300, bootstrap_seed = NULL) {
  # ===================================================================
  # CÁLCULO DE FUNÇÕES DE IMPULSO-RESPOSTA (IRF) PARA MODELO DFM
  # COM WILD BOOTSTRAP SEGUINDO GERTLER & KARADI (2015)
  # ===================================================================
  
  # Set seed for reproducibility
  set_bootstrap_seed(bootstrap_seed)
  
  # Extract components
  Lambda <- dfm_results$static_loadings
  A <- dfm_results$companion_matrix
  K <- dfm_results$dynamic_loadings
  M <- dfm_results$dynamic_scaling
  sy <- dfm_results$data_sd
  p <- dfm_results$p  # Get VAR order from DFM results
  SIGMA <- dfm_results$var_covariance # Covariance matrix of VAR residuals

  # Dimensions
  n_vars <- nrow(Lambda)
  r <- ncol(Lambda)
  q <- dfm_results$q  # Use stored q value directly

  # Point IRF calculation
  
  # Identification via Cholesky decomposition (lower triangular)
  # This ensures consistency with MATLAB's chol, which returns upper by default
  # We use the transpose to get the lower triangular factor
  S <- t(chol(SIGMA))

  # Point IRF
  irf_point <- array(0, dim = c(n_vars, h + 1, q))
  B <- array(0, dim = c(r, r, h + 1))
  B[, , 1] <- diag(r)
  B[, , 2] <- A[1:r, 1:r]

  for (i in 3:(h + 1)) {
    B[, , i] <- A[1:r, 1:r] %*% B[, , i - 1]
  }

  # Calculate point IRF using formula
  # C(:,:,ii)=(lambda*B(:,:,ii)*S*K*M).*repmat(sy',1,q)
  for (i in 1:(h + 1)) {
    if (!is.matrix(K) && !is.matrix(M)) {
      # Case q=r with scalar K,M (K=1, M=1)
      # In this case, q factors but K=1, M=1 so we replicate across q
      temp <- Re(Lambda %*% B[, , i] %*% S * K * M)
      for (s in 1:q) {
        # Each dynamic factor gets the same response (since K=1, M=1)
        irf_point[, i, s] <- temp[, s] * sy  # Extract s-th static factor response
      }
    } else {
      # General case with matrix K,M
      temp <- Re(Lambda %*% B[, , i] %*% S %*% K %*% M)
      for (s in 1:q) {
        irf_point[, i, s] <- temp[, s] * sy  # Element-wise multiplication with sy
      }
    }
  }

  # Wild Bootstrap following Gertler & Karadi (2015) and Alessi & Kerssenfischer methodology
  if (nboot > 0) {
    cat("Iniciando wild bootstrap com", nboot, "replicações...\n")
    irf_boot <- array(0, dim = c(n_vars, h + 1, nboot, q))
    
    # Progress tracking
    pb <- utils::txtProgressBar(min = 0, max = nboot, style = 3)

    for (b in 1:nboot) {
      utils::setTxtProgressBar(pb, b)
      
      tryCatch({
        # Generate wild bootstrap draw (Rademacher: ±1 with prob 0.5 each)
        rr <- 1 - 2 * (runif(nrow(dfm_results$var_residuals)) > 0.5)
        resid_boot <- dfm_results$var_residuals * rr
        
        # Reconstruct factors with bootstrapped residuals
        F_boot <- matrix(0, nrow = nrow(dfm_results$static_factors), ncol = r)
        F_boot[1:p, ] <- dfm_results$static_factors[1:p, ]
        
        for (t in (p + 1):nrow(F_boot)) {
          lagged_vars <- as.vector(t(F_boot[(t-1):(t-p), ]))
          F_boot[t, ] <- lagged_vars %*% dfm_results$var_coefficients + resid_boot[t - p, ]
        }
        
        # Re-estimate only the VAR on bootstrapped factors
        var_boot <- estimate_corrected_var(F_boot, p)
        
        # Re-estimate dynamic components from bootstrapped VAR residuals
        dynamic_boot <- estimate_dynamic_factors(var_boot$residuals, q, r)
        
        # Calculate bootstrapped IRF matrices
        A_boot <- var_boot$companion
        SIGMA_boot <- var_boot$covariance_matrix
        S_boot <- t(chol(SIGMA_boot)) # Cholesky on bootstrapped covariance
        K_boot <- dynamic_boot$K
        M_boot <- dynamic_boot$M
        
        B_boot <- array(0, dim = c(r, r, h + 1))
        B_boot[, , 1] <- diag(r)
        if (h >= 1) B_boot[, , 2] <- A_boot[1:r, 1:r]
        
        for (i in 3:(h + 1)) {
          B_boot[, , i] <- A_boot[1:r, 1:r] %*% B_boot[, , i - 1]
        }

        # Calculate bootstrapped IRFs using original Lambda and sy
        for (i in 1:(h + 1)) {
          if (!is.matrix(K_boot) && !is.matrix(M_boot)) {
            temp <- Re(Lambda %*% B_boot[, , i] %*% S_boot * K_boot * M_boot)
            for (s in 1:q) {
              irf_boot[, i, b, s] <- temp[, s] * sy
            }
          } else {
            temp <- Re(Lambda %*% B_boot[, , i] %*% S_boot %*% K_boot %*% M_boot)
            for (s in 1:q) {
              irf_boot[, i, b, s] <- temp[, s] * sy
            }
          }
        }
      }, error = function(e) {
        warning("Bootstrap iteration ", b, " failed: ", e$message)
        # Use point estimates as fallback for the failed iteration
        for (s in 1:q) {
          irf_boot[, , b, s] <- irf_point[, , s]
        }
      })
    }
    
    close(pb)
    cat("\nWild bootstrap concluído!\n")
    
    # Validate bootstrap results
    bootstrap_validation <- validate_bootstrap_results(irf_boot, irf_point)
    cat("Taxa de sucesso do bootstrap:", round(bootstrap_validation$success_rate * 100, 1), "%\n")

    # Confidence intervals
    irf_ci <- array(0, dim = c(n_vars, h + 1, 5, q))
    for (s in 1:q) {
      irf_ci[, , 1, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.05, na.rm = TRUE)
      irf_ci[, , 2, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.10, na.rm = TRUE)
      irf_ci[, , 3, s] <- irf_point[, , s]
      irf_ci[, , 4, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.90, na.rm = TRUE)
      irf_ci[, , 5, s] <- apply(irf_boot[, , , s], c(1, 2), quantile, probs = 0.95, na.rm = TRUE)
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

#' Wild Bootstrap for Dynamic Factor Models
#'
#' @description
#' Implements wild bootstrap following Gertler & Karadi (2015) methodology
#' as used in Alessi & Kerssenfischer. The procedure reconstructs the full
#' dataset by combining bootstrapped common components with original 
#' idiosyncratic components.
#'
#' @param dfm_results Results from estimate_dfm function
#' @param nboot Number of bootstrap replications
#' @param r Number of static factors  
#' @param q Number of dynamic factors
#' @param p VAR lag order
#'
#' @return Array of bootstrapped DFM results
#' 
#' @details
#' The wild bootstrap procedure:
#' 1. Generates Rademacher draws (±1 with prob 0.5)
#' 2. Applies draws to VAR residuals 
#' 3. Reconstructs factors using bias-corrected VAR
#' 4. Combines with idiosyncratic components
#' 5. Re-estimates complete DFM
#'
#' @references
#' Gertler, M., & Karadi, P. (2015). Monetary policy surprises, credit costs, 
#' and economic activity. American Economic Journal: Macroeconomics, 7(1), 44-76.
wild_bootstrap_dfm <- function(dfm_results, nboot, r, q, p) {
  
  # Pre-calculate components that don't change across bootstrap iterations
  T_data <- nrow(dfm_results$Z)
  regX <- cbind(1, 1:T_data)
  
  # Reconstruct original data from DFM components
  original_data <- dfm_results$Z
  for (i in 1:ncol(original_data)) {
    original_data[, i] <- original_data[, i] * dfm_results$data_sd[i]
  }
  
  # Add back linear trend
  trend_data <- matrix(0, T_data, ncol(original_data))
  for (i in 1:ncol(original_data)) {
    beta_trend <- solve(crossprod(regX)) %*% crossprod(regX, original_data[, i])
    trend_data[, i] <- regX %*% beta_trend
  }
  original_data <- original_data + trend_data
  
  # Calculate common and idiosyncratic components
  common_components <- dfm_results$static_factors %*% t(dfm_results$static_loadings)
  common_components_scaled <- sweep(common_components, 2, dfm_results$data_sd, "*")
  common_components_full <- common_components_scaled + trend_data
  idiosyncratic <- original_data - common_components_full
  
  # Bootstrap array to store results
  bootstrap_results <- list()
  
  for (b in 1:nboot) {
    # Generate wild bootstrap draw
    rr <- 1 - 2 * (runif(nrow(dfm_results$var_residuals)) > 0.5)
    resid_boot <- dfm_results$var_residuals * rr
    
    # Reconstruct factors with bootstrapped residuals
    F_boot <- matrix(0, nrow = nrow(dfm_results$static_factors), ncol = r)
    F_boot[1:p, ] <- dfm_results$static_factors[1:p, ]
    
    for (t in (p + 1):nrow(F_boot)) {
      lagged_vars <- as.vector(t(F_boot[(t-1):(t-p), ]))
      F_boot[t, ] <- lagged_vars %*% dfm_results$var_coefficients + resid_boot[t - p, ]
    }
    
    # Reconstruct data
    common_boot <- F_boot %*% t(dfm_results$static_loadings)
    common_boot_scaled <- sweep(common_boot, 2, dfm_results$data_sd, "*")
    common_boot_full <- common_boot_scaled + trend_data
    X_boot <- common_boot_full + idiosyncratic
    
    # Store bootstrapped dataset
    bootstrap_results[[b]] <- X_boot
  }
  
  return(bootstrap_results)
}

#' Set seed for reproducible wild bootstrap
#' @param seed Integer seed value
set_bootstrap_seed <- function(seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
    cat("Bootstrap seed definido para:", seed, "\n")
  }
}

#' Validate wild bootstrap results
#' @param irf_boot Array of bootstrap IRF results
#' @param irf_point Point estimate IRFs
#' @return List of validation metrics
validate_bootstrap_results <- function(irf_boot, irf_point) {
  n_vars <- dim(irf_boot)[1]
  h <- dim(irf_boot)[2] - 1
  nboot <- dim(irf_boot)[3]
  q <- dim(irf_boot)[4]
  
  # Check for failed bootstrap iterations
  failed_iterations <- 0
  for (b in 1:nboot) {
    for (s in 1:q) {
      if (all(irf_boot[, , b, s] == irf_point[, , s])) {
        failed_iterations <- failed_iterations + 1
        break
      }
    }
  }
  
  # Calculate bootstrap statistics
  bootstrap_stats <- list(
    total_iterations = nboot,
    failed_iterations = failed_iterations,
    success_rate = (nboot - failed_iterations) / nboot,
    mean_abs_irf = mean(abs(irf_boot), na.rm = TRUE),
    bootstrap_variance = var(as.vector(irf_boot), na.rm = TRUE)
  )
  
  if (failed_iterations > nboot * 0.1) {
    warning("Mais de 10% das iterações do bootstrap falharam. Considere ajustar os parâmetros.")
  }
  
  return(bootstrap_stats)
}

#' Compare Wild Bootstrap vs Simple Bootstrap
#' 
#' @description
#' Diagnostic function to compare wild bootstrap results with simple bootstrap
#' for validation purposes. Helps assess the impact of the bootstrap method choice.
#'
#' @param dfm_results Results from estimate_dfm
#' @param h IRF horizon
#' @param nboot Number of bootstrap replications (small number for comparison)
#' @return List comparing both methods
compare_bootstrap_methods <- function(dfm_results, h = 12, nboot = 50) {
  cat("Comparando métodos de bootstrap (", nboot, " replicações)...\n")
  
  # Wild bootstrap (current implementation)
  irf_wild <- compute_irf_dfm(dfm_results, h = h, nboot = nboot, bootstrap_seed = 123)
  
  # Simple bootstrap (original method) - just for comparison
  r <- dfm_results$r
  q <- dfm_results$q
  n_vars <- nrow(dfm_results$static_loadings)
  
  # Point IRF for comparison
  Lambda <- dfm_results$static_loadings
  A <- dfm_results$companion_matrix
  K <- dfm_results$dynamic_loadings
  M <- dfm_results$dynamic_scaling
  sy <- dfm_results$data_sd
  
  irf_point <- array(0, dim = c(n_vars, h + 1, q))
  B <- array(0, dim = c(r, r, h + 1))
  B[, , 1] <- diag(r)
  B[, , 2] <- A[1:r, 1:r]
  
  for (i in 3:(h + 1)) {
    B[, , i] <- A[1:r, 1:r] %*% B[, , i - 1]
  }
  
  for (i in 1:(h + 1)) {
    if (!is.matrix(K) && !is.matrix(M)) {
      temp <- Re(Lambda %*% B[, , i] * K * M)
      for (s in 1:q) {
        irf_point[, i, s] <- temp[, s] * sy
      }
    } else {
      temp <- Re(Lambda %*% B[, , i] %*% K %*% M)
      for (s in 1:q) {
        irf_point[, i, s] <- temp[, s] * sy
      }
    }
  }
  
  # Simple bootstrap for comparison
  irf_simple <- array(0, dim = c(n_vars, h + 1, nboot, q))
  set.seed(123)
  
  for (b in 1:nboot) {
    # Simple residual resampling
    resid_indices <- sample(1:nrow(dfm_results$var_residuals), replace = TRUE)
    resid_boot <- dfm_results$var_residuals[resid_indices, ]
    
    F_boot <- matrix(0, nrow = nrow(dfm_results$static_factors), ncol = r)
    F_boot[1, ] <- dfm_results$static_factors[1, ]
    
    for (t in 2:nrow(F_boot)) {
      F_boot[t, ] <- A[1:r, 1:r] %*% F_boot[t - 1, ] + resid_boot[t - 1, ]
    }
    
    var_boot <- estimate_corrected_var(F_boot, dfm_results$p)
    A_boot <- var_boot$companion
    
    B_boot <- array(0, dim = c(r, r, h + 1))
    B_boot[, , 1] <- diag(r)
    B_boot[, , 2] <- A_boot[1:r, 1:r]
    
    for (i in 3:(h + 1)) {
      B_boot[, , i] <- A_boot[1:r, 1:r] %*% B_boot[, , i - 1]
    }
    
    for (i in 1:(h + 1)) {
      if (!is.matrix(K) && !is.matrix(M)) {
        temp <- Re(Lambda %*% B_boot[, , i] * K * M)
        for (s in 1:q) {
          irf_simple[, i, b, s] <- temp[, s] * sy
        }
      } else {
        temp <- Re(Lambda %*% B_boot[, , i] %*% K %*% M)
        for (s in 1:q) {
          irf_simple[, i, b, s] <- temp[, s] * sy
        }
      }
    }
  }
  
  # Calculate comparison metrics
  wild_variance <- apply(irf_wild[, , 3, ], c(1, 2), var, na.rm = TRUE)
  simple_variance <- apply(irf_simple[, , , ], c(1, 2, 4), var, na.rm = TRUE)
  
  cat("Comparação concluída!\n")
  cat("Variância média - Wild bootstrap:", mean(wild_variance, na.rm = TRUE), "\n")
  cat("Variância média - Simple bootstrap:", mean(simple_variance, na.rm = TRUE), "\n")
  
  return(list(
    wild_bootstrap = irf_wild,
    simple_bootstrap = irf_simple,
    point_estimate = irf_point,
    variance_comparison = list(
      wild = wild_variance,
      simple = simple_variance
    )
  ))
}
