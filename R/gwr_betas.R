#' estimate `gwr` model coefficients
#'
#' @param formula A formula of `gwr` model.
#' @param data An sf object of observation data
#' @param bw (optional) The bandwidth used in selecting models. The optimal bandwidth can be
#' selected using one of two methods: `AIC`, and `CV`. Default will use `AIC`.
#' @param adaptive (optional) Whether the bandwidth value is adaptive or not. Default is `TRUE`.
#' @param kernel (optional) Kernel function. Default is `gaussian`.
#' @param intercept (optional) Whether to include the intercept term in the coefficient `tibble.
#' Default is `FALSE`.
#'
#' @return A `tibble`
#' @export
#'
#' @examples
#' depression = system.file('extdata/Depression.csv',package = 'gdverse') |>
#'   readr::read_csv() |>
#'   sf::st_as_sf(coords = c('X','Y'), crs = 4326)
#' gwr_betas(Depression_prevelence ~ ., data = depression)
#' gwr_betas(Depression_prevelence ~ ., data = depression, intercept = TRUE)
#'
gwr_betas = \(formula, data, bw = "AIC", adaptive = TRUE,
              kernel = "gaussian", intercept = FALSE){
  suppressWarnings({g = GWmodel3::gwr_basic(
    formula, data, bw = bw, adaptive = adaptive, kernel = kernel
  )})
  betas = tibble::as_tibble(stats::coef(g))
  if(!intercept){betas = dplyr::select(betas,-1)}
  return(betas)
}
