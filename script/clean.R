source("R/preprocessing/seasonality.R")
source("R/preprocessing/stationarity.R")


# loading data ----
raw_data <- readr::read_csv("data/raw/raw_data.csv")



# Treating data seasonality ----

season_result <- check_seasonality(raw_data)


data_no_season <- raw_data |>
  dplyr::select(dplyr::matches(season_result$season_vars)) |>
  purrr::map(
    ~ ts(
      .x,
      start = c(
        lubridate::year(min(raw_data$ref.date)),
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
data <- raw_data |>
  dplyr::mutate(
    dplyr::across(
      # Select columns that match the seasonal variables
      dplyr::all_of(season_result$season_vars),
      # Replace with seasonally adjusted values
      ~ data_no_season[[dplyr::cur_column()]]
    )
  )

# Verify the seasonality was removed
check_if_final_was_removed <- check_seasonality(data)
# Yes, it was




# Apply log transformation in nominal variables ----

data <- data |>
  dplyr::mutate(
    dplyr::across(
      .cols = c(
        dplyr::contains("credito"),
        dplyr::contains("consumo"),
        dplyr::contains("veiculos"),
        dplyr::contains("caged"),
        dplyr::contains("pop"),
        dplyr::contains("trab")
      ),
      .fns = ~ log(.x)
    )
  )




# We now gotta check for unity roots ----

unity_root_test <- adf_test(data)
print(unity_root_test)




final_data <- remove_unit_root(data, max_diff = 5)

# readr::write_csv(final_data$data, "data/processed/final_data.csv")
