rm(list = ls())

source("R/preprocessing/seasonality.R")
source("R/preprocessing/stationarity.R")


# loading data ----
raw_data <- readr::read_csv("data/raw/raw_data.csv")



# Apply log transformation in nominal variables ----

data <- raw_data |>
  dplyr::mutate(
    dplyr::across(
      .cols = c(
        dplyr::contains("credito"),
        dplyr::contains("consumo"),
        dplyr::contains("veiculos"),
        dplyr::contains("caged"),
        dplyr::contains("pop"),
        dplyr::contains("trab"),
        -c(trab_tx_desemprego, trab_razao_vagas_desempregados)
      ),
      .fns = ~ log(.x)
    )
  )





# Treating data seasonality ----

season_result <- check_seasonality(data)


data_no_season <- data |>
  dplyr::select(dplyr::matches(season_result$season_vars)) |>
  purrr::map(
    ~ ts(
      .x,
      start = c(
        lubridate::year(min(data$ref.date)),
        1
      ),
      frequency = 12
    ) |>
      seasonal::seas(
        x11 = "",
        transform.function = "auto",
        regression.aictest = NULL,
        outlier = NULL
      ) |>
      seasonal::final()
  ) |>
  purrr::map_dfc(~ as.numeric(.x))


# Create new dataframe using tidyverse approach
data <- data |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(season_result$season_vars),
      ~ data_no_season[[dplyr::cur_column()]]
    )
  )

# Verify the seasonality was removed
check_if_final_was_removed <- check_seasonality(data)
print(check_if_final_was_removed)
# Yes, it was






# We now gotta check for unity roots ----

unity_root_test <- adf_test(data)
print(unity_root_test)



final_data <- remove_unit_root(data, max_diff = 5)

# readr::write_csv(final_data$data, "data/processed/final_data.csv")

final_data$control |> print(n = 100)
