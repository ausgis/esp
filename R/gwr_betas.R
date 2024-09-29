#' Title
#'
#' @param formula
#' @param data
#' @param adaptive
#' @param kernel
#' @param multiscale
#'
#' @return A `tibble`
#' @export
#'
#' @examples
#' depression = system.file('extdata/Depression.csv',package = 'gdverse') |>
#'   readr::read_csv() |>
#'   sf::st_as_sf(coords = c('X','Y'), crs = 4326)
#' gwr_betas(Depression_prevelence ~ ., data = depression)
#' gwr_betas(Depression_prevelence ~ ., data = depression, multiscale = TRUE)
#'
gwr_betas = \(formula, data,
              adaptive = TRUE,
              kernel = "gaussian",
              multiscale = FALSE){
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)
  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }
  yname = formula.vars[1]
  gname = sdsfun::sf_geometry_name(data)
  xname = colnames(data)[-which(colnames(data) %in% c(yname,gname))]

  longlat = dplyr::if_else(sf::st_is_longlat(data),TRUE,FALSE,FALSE)
  if (multiscale) {
    suppressWarnings({g = GWmodel3::gwr_multiscale(
      formula, data,
      config = list(GWmodel3::mgwr_config(adaptive = adaptive,
                                          kernel = kernel,
                                          longlat = longlat))
    )})
  } else {
    suppressWarnings({g = GWmodel3::gwr_basic(
      formula, data, adaptive = adaptive, kernel = kernel
    )})
  }

  betas = coef(g) %>%
    tibble::as_tibble() %>%
    dplyr::select(dplyr::all_of(xname))
  return(betas)
}
