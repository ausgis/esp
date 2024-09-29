gwr_betas = \(formula, data,
              adaptive = TRUE,
              kernel = "gaussian",
              multiscale = FALSE){
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
  return(coef(m))
}
