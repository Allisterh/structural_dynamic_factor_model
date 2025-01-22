#' Teste de Raiz Unitária ADF (Augmented Dickey-Fuller)
#'
#' @description
#' Realiza o teste ADF para verificar estacionariedade em séries temporais.
#' Se encontradas séries não estacionárias, automaticamente realiza o teste na primeira diferença.
#'
#' @param data DataFrame contendo as séries temporais a serem testadas
#' @param type Tipo de teste ADF a ser realizado:
#'   - "none": sem constante nem tendência
#'   - "drift": com constante (padrão)
#'   - "trend": com constante e tendência
#'
#' @return
#' Se todas as séries forem estacionárias em nível:
#'   - DataFrame com resultados do teste contendo: variável, estatística tau,
#'     valor crítico (5%) e resultado (estacionário/não estacionário)
#'
#' Se houver séries não estacionárias:
#'   - Lista com dois DataFrames:
#'     1. "Teste em nivel": resultados do teste em nível
#'     2. "Teste com uma diferenca": resultados do teste na primeira diferença
#'        para séries não estacionárias em nível
#'
#' @details
#' - Aplica o teste apenas em colunas numéricas do DataFrame
#' - Usa valores críticos de 5% de significância
#' - Utiliza os pacotes urca para o teste ADF e dplyr/purrr para manipulação dos dados
#' - Para cada tipo de teste usa diferentes estatísticas:
#'   * none: tau1
#'   * drift: tau2
#'   * trend: tau3
#'
#' @examples
#' # Teste com constante (padrão)
#' resultados <- adf_test(dados)
#'
#' # Teste com constante e tendência
#' resultados <- adf_test(dados, type = "trend")
#'
#' # Acessando resultados quando há séries não estacionárias
#' resultados$`Teste em nivel`
#' resultados$`Teste com uma diferenca`
adf_test <- function(data, type = "drift") {
    get_test_values <- function(teste, type) {
        if (type == "none") {
            # Para type="none", usa tau1
            stat <- teste@teststat[1] # Estatística tau1
            crit <- teste@cval[2] # Valor crítico de 5% para tau1
        } else if (type == "drift") {
            # Para type="drift", usa tau2
            stat <- teste@teststat[1] # Estatística tau2
            crit <- teste@cval[1, 2] # Valor crítico de 5% para tau2
        } else if (type == "trend") {
            # Para type="trend", usa tau3
            stat <- teste@teststat[1] # Estatística tau3
            crit <- teste@cval[1, 2] # Valor crítico de 5% para tau3
        }
        return(list(stat = stat, crit = crit))
    }

    # Aplica o teste ADF para todas as colunas numéricas
    adf_test <- data |>
        dplyr::select(dplyr::where(is.numeric)) |>
        purrr::map(~ urca::ur.df(.x, type = type))

    # Inicializa o dataframe para armazenar os resultados
    df <- data.frame(
        variavel = as.character(),
        tau = as.numeric(),
        valor_critico = as.numeric()
    )

    # Itera pelos resultados e extrai as estatísticas e valores críticos
    for (i in seq_along(adf_test)) {
        test_values <- get_test_values(adf_test[[i]], type)
        df <- rbind(df, data.frame(
            variavel = names(adf_test)[i],
            tau = round(test_values$stat, digits = 3),
            valor_critico = round(test_values$crit, digits = 3)
        ))
    }

    # Classifica as variáveis como estacionárias ou não
    df$resultado <- dplyr::case_when(
        df$tau > df$valor_critico ~ "nao estacionario",
        df$tau <= df$valor_critico ~ "estacionario"
    )

    # Se houver variáveis não estacionárias, realiza teste na 1ª diferença
    if (any(df$resultado == "nao estacionario")) {
        variaveis_ne <- df |>
            dplyr::filter(resultado == "nao estacionario") |>
            dplyr::pull(var = "variavel") # Extrai as variáveis não estacionárias

        # Aplica o teste ADF na 1ª diferença
        adf_test2 <- data |>
            dplyr::select(dplyr::any_of(variaveis_ne)) |>
            purrr::map(~ diff(.x) |>
                urca::ur.df(type = type))

        # Inicializa outro dataframe para armazenar os novos resultados
        df2 <- data.frame(
            variavel = as.character(),
            tau = as.numeric(),
            valor_critico = as.numeric()
        )

        # Itera pelos novos resultados e extrai estatísticas e valores críticos
        for (i in seq_along(adf_test2)) {
            test_values <- get_test_values(adf_test2[[i]], type)
            df2 <- rbind(df2, data.frame(
                variavel = names(adf_test2)[i],
                tau = round(test_values$stat, digits = 3),
                valor_critico = round(test_values$crit, digits = 3)
            ))
        }

        # Classifica os resultados da 1ª diferença
        df2$resultado <- dplyr::case_when(
            df2$tau > df2$valor_critico ~ "nao estacionario",
            df2$tau <= df2$valor_critico ~ "estacionario"
        )

        # Retorna os resultados para nível e 1ª diferença
        return(list(
            "Teste em nivel" = df,
            "Teste com uma diferenca" = df2
        ))
    }

    # Retorna os resultados do teste em nível
    return(df)
}



#' Remove Unit Root Through Sequential Differencing
#'
#' @description
#' Applies sequential differencing to non-stationary variables until stationarity
#' is achieved or maximum differences are reached. Uses Augmented Dickey-Fuller test
#' to check for stationarity.
#'
#' @param data A data frame containing time series variables to be tested and
#'   potentially differenced
#' @param max_diff Maximum number of differences to apply. Default is 3
#'
#' @return A list containing two elements:
#'   \itemize{
#'     \item data: The transformed data frame with differenced variables
#'     \item control: A tibble tracking which variables were differenced and how many times
#'   }
#'
#' @details
#' The function implements the following algorithm:
#' 1. Tests all variables for stationarity using ADF test
#' 2. Identifies non-stationary variables
#' 3. Applies first difference to non-stationary variables
#' 4. Repeats process until either:
#'    - All variables become stationary
#'    - Maximum number of differences (max_diff) is reached
#'
#' First observation for each difference is set to NA due to the differencing process.
#'
#' @examples
#' \dontrun{
#' # Create sample non-stationary data
#' data <- data.frame(
#'     date = seq.Date(
#'         from = as.Date("2020-01-01"),
#'         by = "month", length.out = 24
#'     ),
#'     random_walk = cumsum(rnorm(24)),
#'     stationary = rnorm(24)
#' )
#'
#' # Apply unit root removal
#' result <- remove_unit_root(data, max_diff = 2)
#'
#' # Check which variables were differenced
#' print(result$control)
#' }
#'
#' @importFrom dplyr filter pull mutate across all_of bind_rows
#' @importFrom tidyr drop_na
#' @importFrom tibble tibble
#'
#' @export
remove_unit_root <- function(data, max_diff = 3) {
    df_temp <- data
    df_control <- tibble::tibble(
        variable = character(),
        times_diff = numeric()
    )

    for (i in 1:max_diff) {
        test_data <- df_temp |>
            tidyr::drop_na()

        test_result <- adf_test(test_data)

        if (!is.list(test_result) || !("Teste em nivel" %in% names(test_result))) {
            test_result <- list("Teste em nivel" = test_result)
        }

        names_unit_root <- test_result[["Teste em nivel"]] |>
            dplyr::filter(resultado == "nao estacionario") |>
            dplyr::pull(variavel)

        if (length(names_unit_root) == 0) break

        df_temp <- df_temp |>
            dplyr::mutate(
                dplyr::across(
                    .cols = dplyr::all_of(names_unit_root),
                    .fns = ~ c(NA, diff(.x))
                )
            )

        new_controls <- tibble::tibble(
            variable = names_unit_root,
            times_diff = i
        )

        df_control <- dplyr::bind_rows(
            df_control,
            new_controls
        )
    }

    return(list("data" = df_temp, "control" = df_control))
}
