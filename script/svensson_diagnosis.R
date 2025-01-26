rm(list = ls())

source("R/modeling/svensson_model.R")

# Definir os títulos de interesse
titulos <- c("LTN")

# Baixar os dados dos títulos
dados_tesouro <- GetTDData::td_get(titulos) |> suppressMessages()


# Aplicar o modelo aos dados
resultado <- generate_fixed_maturity_series(dados_tesouro)






library(ggplot2)
library(dplyr)
library(tidyr)
library(tseries) # For additional time series tests
library(forecast) # For time series decomposition
library(corrplot) # For correlation visualization
library(moments) # For distribution analysis
library(zoo) # For rolling calculations

run_enhanced_diagnostics <- function(resultado) {
  # Logging function for diagnostic messages
  log_msg <- function(msg) {
    cat(paste0("\n", msg, "\n"))
  }

  # Safe wrapper for potentially failing operations
  safe_compute <- function(expr) {
    tryCatch(
      expr,
      error = function(e) {
        log_msg(paste("Warning: Operation failed -", e$message))
        return(NULL)
      }
    )
  }

  # 1. Basic Time Series Visualization ----
  log_msg("Computing time series visualization...")
  p1 <- safe_compute({
    ggplot(resultado, aes(x = data)) +
      geom_line(aes(y = titulo_0_25ano, color = "3 meses")) +
      geom_line(aes(y = titulo_1ano, color = "1 ano")) +
      geom_line(aes(y = titulo_2ano, color = "2 anos")) +
      geom_line(aes(y = titulo_3ano, color = "3 anos")) +
      geom_line(aes(y = titulo_5ano, color = "5 anos")) +
      labs(
        title = "Evolução das Taxas de Juros",
        y = "Taxa",
        x = "Data",
        color = "Maturidade"
      ) +
      theme_minimal() +
      scale_color_manual(values = c(
        "3 meses" = "purple",
        "1 ano" = "blue",
        "2 anos" = "red",
        "3 anos" = "green",
        "5 anos" = "orange"
      ))
  })

  # 2. Term Structure Analysis ----
  log_msg("Computing term structure statistics...")
  term_structure <- safe_compute({
    resultado %>%
      select(starts_with("titulo")) %>%
      summarise(across(
        everything(),
        list(
          mean = ~ mean(., na.rm = TRUE),
          sd = ~ sd(., na.rm = TRUE),
          q25 = ~ quantile(., 0.25, na.rm = TRUE),
          q75 = ~ quantile(., 0.75, na.rm = TRUE)
        )
      ))
  })

  # 3. Enhanced Spread Analysis ----
  log_msg("Computing spread analysis...")
  spreads <- safe_compute({
    resultado %>%
      mutate(
        spread_2y1y = titulo_2ano - titulo_1ano,
        spread_3y1y = titulo_3ano - titulo_1ano,
        spread_5y1y = titulo_5ano - titulo_1ano,
        spread_curve = titulo_5ano - titulo_0_25ano
      )
  })

  spread_stats <- safe_compute({
    if (!is.null(spreads)) {
      spreads %>%
        select(starts_with("spread")) %>%
        summary()
    }
  })

  # 4. Distribution Analysis ----
  log_msg("Computing distribution characteristics...")
  rate_distributions <- safe_compute({
    resultado %>%
      select(starts_with("titulo")) %>%
      summarise(across(
        everything(),
        list(
          skewness = ~ skewness(., na.rm = TRUE),
          kurtosis = ~ kurtosis(., na.rm = TRUE)
        )
      ))
  })

  # Add Jarque-Bera test separately to handle errors
  jb_tests <- safe_compute({
    sapply(resultado %>% select(starts_with("titulo")), function(x) {
      tryCatch(
        jarque.bera.test(na.omit(x))$statistic,
        error = function(e) NA
      )
    })
  })

  # 5. Enhanced Correlation Analysis ----
  log_msg("Computing correlation analysis...")
  cor_matrix <- safe_compute({
    cor(resultado %>% select(starts_with("titulo")),
      use = "pairwise.complete.obs"
    )
  })

  p_corr <- safe_compute({
    if (!is.null(cor_matrix)) {
      corrplot(cor_matrix,
        method = "color",
        type = "upper",
        addCoef.col = "black"
      )
    }
  })

  # 6. Volatility Analysis ----
  log_msg("Computing volatility analysis...")
  vol_analysis <- safe_compute({
    resultado %>%
      mutate(across(starts_with("titulo"),
        list(vol = ~ rollapply(.,
          width = 21,
          FUN = function(x) sd(x, na.rm = TRUE) * sqrt(252),
          fill = NA, align = "right",
          partial = TRUE # Allow partial windows
        )),
        .names = "{col}_vol"
      ))
  })

  p_vol <- safe_compute({
    if (!is.null(vol_analysis)) {
      vol_analysis %>%
        select(data, ends_with("vol")) %>%
        pivot_longer(cols = -data) %>%
        ggplot(aes(x = data, y = value, color = name)) +
        geom_line() +
        labs(
          title = "Volatilidade Anualizada (21 dias)",
          y = "Volatilidade",
          x = "Data"
        ) +
        theme_minimal()
    }
  })

  # 7. Seasonal Analysis ----
  log_msg("Computing seasonal patterns...")
  seasonal_analysis <- safe_compute({
    resultado %>%
      mutate(
        month = format(data, "%m"),
        year = format(data, "%Y")
      ) %>%
      group_by(month) %>%
      summarise(across(
        starts_with("titulo"),
        list(mean = ~ mean(., na.rm = TRUE))
      ))
  })

  # 8. Principal Component Analysis ----
  log_msg("Computing PCA...")
  pca_results <- safe_compute({
    complete_data <- na.omit(resultado %>% select(starts_with("titulo")))
    if (nrow(complete_data) > 0) {
      prcomp(complete_data, scale. = TRUE)
    }
  })

  variance_explained <- safe_compute({
    if (!is.null(pca_results)) {
      summary(pca_results)$importance[2, ] * 100
    }
  })

  # 9. Mean Reversion Tests ----
  log_msg("Computing stationarity tests...")
  adf_tests <- safe_compute({
    lapply(
      resultado %>% select(starts_with("titulo")),
      function(x) {
        tryCatch(
          adf.test(na.omit(x)),
          error = function(e) NULL
        )
      }
    )
  })

  # 10. Model Fit Quality ----
  log_msg("Computing model fit quality metrics...")

  # Improved R-squared calculation using all available data
  daily_r2 <- safe_compute({
    resultado %>%
      group_by(data) %>%
      summarise(
        r2 = if (sum(!is.na(c_across(starts_with("titulo")))) >= 2) {
          # Calculate R² only if we have at least 2 valid observations
          cor_matrix <- cor(c_across(starts_with("titulo")),
            use = "pairwise.complete.obs"
          )
          mean(cor_matrix[upper.tri(cor_matrix)]^2, na.rm = TRUE)
        } else {
          NA_real_
        },
        .groups = "drop"
      )
  })

  # 11. Generate Comprehensive Report ----
  log_msg("Generating comprehensive report...")

  cat("\n=== Enhanced Svensson Model Diagnostics Report ===\n")

  if (!is.null(term_structure)) {
    cat("\n1. Term Structure Statistics:\n")
    print(term_structure)
  }

  if (!is.null(spread_stats)) {
    cat("\n2. Spread Analysis:\n")
    print(spread_stats)
  }

  if (!is.null(rate_distributions)) {
    cat("\n3. Distribution Characteristics:\n")
    print(rate_distributions)
  }

  if (!is.null(variance_explained)) {
    cat("\n4. Principal Components Analysis:\n")
    cat("Variance explained by components:\n")
    print(variance_explained)
  }

  if (!is.null(adf_tests)) {
    cat("\n5. Mean Reversion Tests (p-values):\n")
    print(sapply(adf_tests, function(x) if (!is.null(x)) x$p.value else NA))
  }

  if (!is.null(daily_r2)) {
    cat("\n6. Model Fit Quality:\n")
    cat("Average R-squared:", mean(daily_r2$r2, na.rm = TRUE))
  }

  # 12. Data Quality Metrics ----
  log_msg("Computing data quality metrics...")
  data_quality <- safe_compute({
    resultado %>%
      summarise(across(
        starts_with("titulo"),
        list(
          missing = ~ sum(is.na(.)),
          pct_missing = ~ mean(is.na(.)) * 100,
          consecutive_missing = ~ max(rle(is.na(.))$lengths)
        )
      ))
  })

  # Return plots and analysis objects
  return(list(
    plots = list(
      rates_evolution = p1,
      correlation = p_corr,
      volatility = p_vol
    ),
    analysis = list(
      term_structure = term_structure,
      spreads = spread_stats,
      distributions = rate_distributions,
      pca = pca_results,
      adf_tests = adf_tests,
      fit_quality = daily_r2,
      data_quality = data_quality
    )
  ))
}

# Example usage:
results <- run_enhanced_diagnostics(resultado)

# Function to plot residual analysis
plot_residuals <- function(resultado) {
  residuals_long <- resultado %>%
    select(data, starts_with("titulo")) %>%
    pivot_longer(
      cols = -data,
      names_to = "maturity",
      values_to = "rate"
    ) %>%
    group_by(maturity) %>%
    mutate(residual = rate - mean(rate, na.rm = TRUE))

  ggplot(residuals_long, aes(x = data, y = residual, color = maturity)) +
    geom_line() +
    facet_wrap(~maturity, scales = "free_y") +
    labs(
      title = "Residuals Analysis by Maturity",
      x = "Date",
      y = "Residual"
    ) +
    theme_minimal()
}

plot_residuals(resultado)

# Function to analyze curve inversions
analyze_inversions <- function(resultado) {
  inversions <- resultado %>%
    mutate(
      inversion_2y1y = titulo_2ano < titulo_1ano,
      inversion_3y1y = titulo_3ano < titulo_1ano,
      inversion_5y1y = titulo_5ano < titulo_1ano
    ) %>%
    summarise(
      across(
        starts_with("inversion"),
        list(
          total_days = ~ sum(., na.rm = TRUE),
          pct_time = ~ mean(., na.rm = TRUE) * 100
        )
      )
    )

  return(inversions)
}


analyze_inversions(resultado)
