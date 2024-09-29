esp = \(formula, data, discvar = NULL, discnum = 3:8,
        alpha = 0.75, bw = "AIC", adaptive = TRUE,
        kernel = "gaussian", overlay = 'and', ...) {
  formula = stats::as.formula(formula)
  formula.vars = all.vars(formula)
  if (formula.vars[2] != "."){
    data = dplyr::select(data,dplyr::all_of(formula.vars))
  }
  yname = formula.vars[1]
  geom = sf::st_geometry(data)
  gdist = sdsfun::sf_distance_matrix(data)
  gname = sdsfun::sf_geometry_name(data)
  xname = colnames(data)[-which(colnames(data) %in% c(yname,gname))]

  if (is.null(discvar)) {
    xdiscname = xname
    xundiscname = NULL
  } else {
    xdiscname = discvar
    xundiscname = xname[-which(xname %in% discvar)]
  }
  discdf = dplyr::select(data,dplyr::all_of(c(yname,xdiscname)))

  gwrcoefs = esp::gwr_betas(paste0(yname," ~ ."),discdf,bw,adaptive,kernel)
  discdf = purrr::map_dfc(gwrcoefs, \(.coef) {
    D0 = stats::dist(.coef)
    res = ClustGeo::hclustgeo(D0,as.dist(gdist),alpha)
    return(as.factor(stats::cutree(res,.k)))
  })
  return(discdf)
}

