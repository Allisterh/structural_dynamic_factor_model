#' Test for Seasonality in Numeric Variables
#'
#' @description
#' Analyzes numeric columns in a data frame for seasonal patterns using the
#' `seastests::isSeasonal()` function. The test is performed with a monthly
#' frequency (freq = 12).
#'
#' @param data A data frame containing numeric columns to be tested for seasonality
#'
#' @return A list containing two elements:
#'   \itemize{
#'     \item season_test_result: A data frame with test results for each numeric variable
#'     \item season_vars: A character vector containing names of variables that
#'           exhibit seasonal patterns
#'   }
#'
#' @details
#' The function performs the following steps:
#' 1. Selects numeric columns from the input data frame
#' 2. Applies seasonality test with monthly frequency
#' 3. Returns both detailed test results and a simplified list of seasonal variables
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'     date = seq.Date(
#'         from = as.Date("2020-01-01"),
#'         by = "month", length.out = 24
#'     ),
#'     sales = rnorm(24),
#'     temperature = sin(seq(0, 4 * pi, length.out = 24))
#' )
#' results <- check_seasonality(data)
#' print(results$season_vars) # Variables with seasonal patterns
#' }
#'
#' @importFrom dplyr select_if filter pull
#' @importFrom purrr map_dfr
#' @importFrom tidyr pivot_longer
#' @importFrom seastests isSeasonal
#'
#' @export
check_seasonality <- function(data) {
    df <- data |>
        dplyr::select_if(is.numeric) |>
        purrr::map_dfr(~ seastests::isSeasonal(., freq = 12),
            .id = "variable"
        ) |>
        tidyr::pivot_longer(
            cols = dplyr::everything(),
            names_to = "variavel",
            values_to = "is_seasonal"
        )


    sasonality <- df |>
        dplyr::filter(is_seasonal == TRUE) |>
        dplyr::pull(variavel)

    return(list("season_test_result" = df, "season_vars" = sasonality))
}
