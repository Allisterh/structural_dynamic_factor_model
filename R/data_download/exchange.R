#' Download historical exchange rate data
#'
#' @description
#' Function that downloads historical exchange rate data for different currencies
#' relative to a base currency, using Yahoo Finance data through the tidyquant package.
#'
#' @param currency_tickers Character vector with the desired currency codes (e.g., c("USD", "EUR"))
#' @param base_currency Base currency for conversion. Default is "BRL" (Brazilian Real)
#' @param start_date Start date in "YYYY-MM-DD" format. Default is "2010-01-01"
#' @param end_date End date in "YYYY-MM-DD" format. Default is "2024-12-01"
#'
#' @return
#' Returns a dataframe with monthly data containing:
#' - date: Date (first day of the month)
#' - cambio_[currency]: Columns with monthly average exchange rate for each currency
#'
#' @details
#' The function:
#' 1. Formats currency tickers to Yahoo Finance standard
#' 2. Downloads daily exchange rate data using tidyquant
#' 3. Selects only date, symbol and closing price
#' 4. Pivots data to wide format
#' 5. Renames columns adding "cambio_" prefix
#' 6. Calculates monthly averages
#'
#' @examples
#' # Download Dollar and Euro exchange rates
#' df <- download_cambio(c("USD", "EUR"))
#'
#' # Download Dollar exchange rate relative to Euro
#' df_eur <- download_cambio("USD", base_currency = "EUR")
#'
#' @importFrom tidyquant tq_get
#' @importFrom dplyr select rename_with mutate group_by summarise arrange across
#' @importFrom tidyr pivot_wider
#' @importFrom lubridate floor_date
download_cambio <- function(
    currency_tickers,
    base_currency = "BRL",
    start_date = "2010-01-01",
    end_date = "2024-12-01") {
    currency_tickers <- paste0(currency_tickers, base_currency, "=X") |> tolower()

    df <- tidyquant::tq_get(
        currency_tickers,
        from = start_date,
        to = end_date,
        get = "stock.prices"
    ) |>
        dplyr::select(date, symbol, close) |>
        tidyr::pivot_wider(
            names_from = symbol,
            values_from = close
        ) |>
        dplyr::rename_with(
            ~ paste0("cambio_", gsub("brl=x", "", .x)),
            -date
        )

    df_mensal <- df |>
        dplyr::mutate(date = lubridate::floor_date(date, "month")) |>
        dplyr::group_by(date) |>
        dplyr::summarise(
            across(
                starts_with("cambio_"),
                ~ mean(.x, na.rm = TRUE)
            )
        ) |>
        dplyr::arrange(date)

    return(df_mensal)
}
