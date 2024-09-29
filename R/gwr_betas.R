#' Estimated Model Coefficients By `GWR`
#'
#' @param formula A formula of `GWR` model.
#' @param data An sf object of observation data
#' @param bw (optional) The bandwidth used in selecting models. The optimal bandwidth can be selected
#' using one of two methods: `AIC`, and `CV`. Default will use `AIC`.
#' @param adaptive (optional) Whether the bandwidth value is adaptive or not. Default is `TRUE`.
#' @param kernel (optional) Kernel function. Default is `gaussian`.
#'
#' @return A `tibble`
#' @export
#'
#' @examples
#' depression = system.file('extdata/Depression.csv',package = 'gdverse') |>
#'   readr::read_csv() |>
#'   sf::st_as_sf(coords = c('X','Y'), crs = 4326)
#' gwr_betas(Depression_prevelence ~ ., data = depression)
#'
gwr_betas = \(formula, data, bw = "AIC",
              adaptive = TRUE, kernel = "gaussian"){
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)
  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }
  yname = formula.vars[1]
  gname = sdsfun::sf_geometry_name(data)
  xname = colnames(data)[-which(colnames(data) %in% c(yname,gname))]

  suppressWarnings({g = GWmodel3::gwr_basic(
    formula, data, bw = bw, adaptive = adaptive, kernel = kernel
  )})

  betas = coef(g) |>
    tibble::as_tibble() |>
    dplyr::select(dplyr::all_of(xname))
  return(betas)
}
